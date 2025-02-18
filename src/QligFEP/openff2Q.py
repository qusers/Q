"""Module containing the OpenFF2Q class to process ligands and generate OpenFF parameter files for QligFEP."""

from io import StringIO
from pathlib import Path
from typing import Optional, TextIO

import numpy as np
from joblib import Parallel, delayed, parallel_config
from openff.toolkit import ForceField, Molecule, Topology
from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
from tqdm import tqdm

from .chemIO import MoleculeIO
from .IO import get_force_field_paths, parse_prm
from .logger import logger
from .pdb_utils import (
    append_pdb_to_another,
    pdb_parse_out,
    read_pdb_to_dataframe,
    write_dataframe_to_pdb,
)
from .settings.settings import FF_DIR


class OpenFF2Q(MoleculeIO):
    """Class to process ligands and generate OpenFF parameter files for QligFEP. Dictionary
    variables use ligand names as keys, as setup in the MoleculeIO class.

    Attributes:
        mapping: Dictionary with enumerated list of atoms containing [index, name, type, charge, Xcoord, Ycoord, Zcoord]
        forcefield: OpenFF ForceField object.
        topologies: Dict with mol names as keys and `openff.toolkit.topology.Topology` objects as values.
        parameters: Dict with mol names as keys and `openff.toolkit.typing.engines.smirnoff.parameters.ParameterList`
            objects as values.
        charges_list_magnitude: Dictionary to store the partial charge magnitude for each
            atom in the ligand for each ligand.
        total_charges: Dictionary to store the total charges for each ligand.
    """

    def __init__(
        self, lig, pattern: str = "*.sdf", reindex_hydrogens: bool = True, nagl: bool = False, n_jobs: int = 1
    ):
        """Initializes a new instance of OpenFF2Q to process the `.sdf` input as lig.

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to None.
            reindex_hydrogens: If True, loading molecules will assert that hydrogen atoms are at the end
                of the atom list and reindex them if they are not (needed by restraint setting algorithm).
                If False, the molecules will be loaded as is. Defaults to True.
            nagl: if True, the partial charges will be calculated with nagl, which does so with
                deep learning in the backend, and is faster than the default AM1-BCC method.
            n_jobs: number of jobs to calculate the partial charges in parallel.
        """
        super().__init__(lig, pattern=pattern, reindex_hydrogens=reindex_hydrogens)
        self.n_jobs = n_jobs
        self.out_dir = Path(self.lig).parent
        self.mapping = {lname: {} for lname in self.lig_names}
        self.forcefield = self._set_forcefield(None)
        self.topologies, self.parameters = self.set_topologies_and_parameters()
        self.charges_list_magnitude = {}  # store charge magnitude for each ligand
        self.total_charges = {}  # store the total charges
        self._set_nagl(nagl=nagl)

    def _set_forcefield(self, ffstring: Optional[str]) -> ForceField:
        if ffstring is None:
            # why not the constrained: https://docs.openforcefield.org/projects/toolkit/en/stable/faq.html
            ffstring = "openff-2.2.1.offxml"
            # for a list of the forcefields: https://github.com/openforcefield/openff-forcefields/tree/main/openforcefields/offxml
        logger.debug(f"Forcefield for the ligand parameters: {ffstring}")
        return ForceField(ffstring)

    def _set_nagl(self, nagl: bool):
        """Set the forcefield to be used to calculate the molecules' partial charges."""
        if nagl:
            self.nagl = NAGLToolkitWrapper()
            # Implementation following the OMSF demo @ Naturalis, Leiden:
            # https://github.com/openforcefield/symposium_2024_demo/tree/main
            logger.debug("Warming up the NAGL toolkit")
            self.nagl_partial_charg_str = "openff-gnn-am1bcc-0.1.0-rc.2.pt"
            self.nagl.assign_partial_charges(Molecule.from_smiles("C"), self.nagl_partial_charg_str)
        else:
            self.nagl_partial_charg_str = None
            self.nagl = None

    def _assign_charge(self, molecule: Molecule) -> np.ndarray:
        """Private method that assigns partial charges to an input Molecule and returns their magnitudes"""
        try:
            if self.nagl is not None:
                # Nagl is much faster than the default amber method, but not all ligands are supported
                self.nagl.assign_partial_charges(molecule, partial_charge_method=self.nagl_partial_charg_str)
            else:
                # Seems like the fastest way of doing this with Amber (if you're not OpenEye licensed);
                # for details on this, see: https://github.com/openforcefield/openff-toolkit/issues/1853
                molecule.assign_partial_charges(partial_charge_method="am1bcc")
        except Exception as e:
            print(f"Failed to assign charges for a molecule: {e}")
        charges_magnitudes = np.array([c._magnitude for c in molecule.partial_charges])
        return charges_magnitudes

    def set_topologies_and_parameters(self):
        topologies = {}
        parameters = {}
        for lname, mol in zip(self.lig_names, self.molecules):
            topology = Topology.from_molecules(mol)
            topologies.update({lname: topology})
            parameters.update({lname: self.forcefield.label_molecules(topology)[0]})
        return topologies, parameters

    def process_ligands(self) -> None:
        """Assigns partial charges and writes the .lib, .prm and .pdb files for each ligand."""
        logger.info("Calculating charges")
        backend = "threading" if self.nagl else "multiprocessing"
        with parallel_config(n_jobs=self.n_jobs, backend=backend):
            molecules = tqdm(self.molecules)
            charges_magnitudes = Parallel()(delayed(self._assign_charge)(molecule) for molecule in molecules)
        logger.info("Done! Writing .lib, .prm and .pdb files for each ligand")
        logger.debug(f"Output path: {self.out_dir}")
        for lname, charges in zip(self.lig_names, charges_magnitudes):
            self.charges_list_magnitude.update({lname: charges})
            formatted_sum = f"{round(charges.sum(), 10):.3f}"
            if formatted_sum == "-0.000":
                formatted_sum = "0.000"
            self.total_charges.update({lname: formatted_sum})
            self.create_atom_prm_mapping(lname)
        all_formal_charges = [self.total_charges[n] for n in self.lig_names]
        if np.unique(all_formal_charges).size > 1:
            logger.warning(f"Formal charges of ligands in .sdf are not unique: {self.total_charges}")
        else:
            logger.info(f"Output files written for {len(self.lig_names)} ligands")

    def write_ligand_files(self, prefix=None, residue_name: str = "LIG") -> None:
        """Writes the .lib, .prm and .pdb files for each ligand.

        Args:
            prefix: a custom prefix to be added to the atom names on the .lib file. Element symbols will be
                lowercase not to be confused with protein atoms. Prefix should be used if more than one
                ligand will be present in the system. E.g.: "X", used as standard for ligand 2 in QligFEP.
                Defaults to None.
            residue_name: name of the residue in the .lib file, matching the residue on the pdb file.
                Defaults to "LIG".
        """
        for name, _ in self:
            self.write_lib_Q(name, prefix=prefix, residue_name=residue_name)
            self.write_prm_Q(name, prefix=prefix)
            self.write_PDB(name, residue_name=residue_name)

    def create_atom_prm_mapping(self, lname):
        """Get the mapping of the ligand atoms to the forcefield parameters.

        Args:
            lname: name of the ligand for the mapping.
        """
        # Splitting the SDF content into individual entries
        sdf_content = self.sdf_contents[lname]
        if len(sdf_content) < 4:
            logger.error("Entry is too short to have atom data!")

        # Safely extract atom count
        try:
            atom_count = int(sdf_content[3].split()[0])
        except ValueError:
            # Skip this ligand if counts line is not properly formatted
            logger.error("Could not extract atom count from SDF file!")

        cnt = -1

        for line in sdf_content[4 : 4 + atom_count]:  # Process only atom lines
            cnt += 1
            atom_data = line.split()

            if len(atom_data) < 4:
                continue  # Skip if line does not have enough data

            atom_index = cnt + 1
            charge = round(self.charges_list_magnitude[lname][cnt], 3)  # round for Q
            self.mapping[lname][cnt] = [
                str(atom_index),  # atom index
                atom_data[3] + str(atom_index),  # atom name
                atom_data[3],  # atom type
                str(charge),  # charge
                atom_data[0],  # X coordinate
                atom_data[1],  # Y coordinate
                atom_data[2],  # Z coordinate
            ]

    def write_lib_Q(
        self,
        lname: str,
        outfile: Optional[TextIO] = None,
        prefix: Optional[str] = None,
        residue_name: str = "LIG",
    ):
        """Writes Q's .lib file for a given ligand.

        Args:
            lname: name of the ligand for the .lib file.
            outfile: file to write the .lib file to. Defaults to None, writing self.out_dir / lname.lib.
            prefix: a custom prefix to be added to the atom names on the .lib file. Element symbols will be
                lowercase not to be confused with protein atoms. Prefix should be used if more than one
                ligand will be present in the system. E.g.: "X", used as standard for ligand 2 in QligFEP.
                Defaults to None.
            residue_name: name of the residue in the .lib file, matching the residue on the pdb file.
                Defaults to "LIG".
        """
        parameters = self.parameters[lname]
        mapping = self.mapping[lname]
        total_charge = self.total_charges[lname]

        if outfile is None:
            libfile_out = self.out_dir / f"{lname}.lib"
            outfile = open(libfile_out, "w")  # noqa: SIM115
            should_close = True
        else:
            should_close = False

        try:
            label = "{" + residue_name + "}"
            outfile.write(f"{label}    ! atoms no {len(mapping)}   total charge {total_charge} \n\n")

            outfile.write("[info] \n SYBYLtype RESIDUE \n\n")

            # atom and charge block:
            outfile.write("[atoms] \n")
            for at in mapping:
                at_name = prefix + mapping[at][1].lower() if prefix is not None else mapping[at][1].lower()
                outfile.write(
                    f"{mapping[at][0]:>4s}   {mapping[at][1]:10}{at_name:11}{mapping[at][3]:>10s}\n"
                )
            # bonded block
            outfile.write("\n[bonds]\n")
            for bond in parameters["Bonds"]:
                ai = mapping[bond[0]][1]
                aj = mapping[bond[1]][1]
                outfile.write(f"{ai:10s}{aj}\n")

            # improper block
            outfile.write("\n[impropers]\n")
            for torsion in parameters["ImproperTorsions"]:
                ai = mapping[torsion[0]][1]
                aj = mapping[torsion[1]][1]
                ak = mapping[torsion[2]][1]
                al = mapping[torsion[3]][1]
                outfile.write(f"{ai:10}{aj:10}{ak:10}{al}\n")
        finally:
            if should_close:
                outfile.close()

    def write_prm_Q(self, lname: str, outfile: Optional[TextIO] = None, prefix: Optional[str] = None):
        """Writes Q's .prm file for a given ligand.

        Args:
            lname: name of the ligand for the .prm file.
            outfile: file to write the .prm file to. Defaults to None, writing self.out_dir / lname.prm.
            prefix: a custom prefix to be added to the atom names on the .prm file. Element symbols will be
                lowercase not to be confused with protein atoms. Prefix should be used if more than one
                ligand will be present in the system. E.g.: "X", used as standard for ligand 2 in QligFEP.
                Defaults to None.
            residue_name: name of the residue in the .prm file, matching the residue on the pdb file.
                Defaults to "LIG".
        """

        def insert_prefix(at_name, prefix: Optional[str]):
            if prefix is not None:
                return prefix + at_name
            return at_name

        prm_template = str(FF_DIR / "NOMERGE.prm")
        parameters = self.parameters[lname]
        mapping = self.mapping[lname]

        if outfile is None:
            prmfile_out = f"{self.out_dir / lname}.prm"
            outfile = open(prmfile_out, "w")  # noqa: SIM115
            should_close = True
        else:
            should_close = False

        try:
            mol = self[lname]
            with open(prm_template) as infile:
                for line in infile:
                    block = 0
                    outfile.write(line)
                    if len(line) > 1:
                        if line == "! Ligand vdW parameters\n":
                            block = 1
                        if line == "! Ligand bond parameters\n":
                            block = 2
                        if line == "! Ligand angle parameters\n":
                            block = 3
                        if line == "! Ligand torsion parameters\n":
                            block = 4
                        if line == "! Ligand improper parameters\n":
                            block = 5

                    if block == 1:
                        for atom_indices, parameter in parameters["vdW"].items():
                            ai = atom_indices[0]
                            ai_name = insert_prefix(mapping[ai][1].lower(), prefix)
                            # This is a bit hacky, check how to get the float out directly
                            epsilon = float(f"{parameter.epsilon}".split()[0])
                            epsilon23 = epsilon / 2
                            # TO DO: CHECK IF THIS IS CORRECT!
                            Rmin = f"{parameter.rmin_half}"
                            Rmin = Rmin.split()[0]
                            Rmin = float(Rmin)
                            assert len(atom_indices) == 1, f"More than 1 atom indices present: {atom_indices}"
                            mass = f"{round(mol.atoms[atom_indices[0]].mass.magnitude, 2):.2f}"
                            outfile.write(
                                f"{ai_name:<6}{Rmin:>16.3f}{0.00: 12.1f}{epsilon:>10.3f}{Rmin:>12.3f}{epsilon23:>11.5f}{mass:>11s}\n"
                            )

                    if block == 2:
                        for atom_indices, parameter in parameters["Bonds"].items():
                            ai = atom_indices[0]
                            ai_name = insert_prefix(mapping[ai][1].lower(), prefix)
                            aj = atom_indices[1]
                            aj_name = insert_prefix(mapping[aj][1].lower(), prefix)
                            fc = float(f"{parameter.k}".split()[0])
                            leng = float(f"{parameter.length}".split()[0])
                            outfile.write(f"{ai_name:13}{aj_name:13}{fc:10.1f}{leng:>11.4f}\n")

                    if block == 3:
                        for atom_indices, parameter in parameters["Angles"].items():
                            ai = atom_indices[0]
                            ai_name = insert_prefix(mapping[ai][1].lower(), prefix)
                            aj = atom_indices[1]
                            aj_name = insert_prefix(mapping[aj][1].lower(), prefix)
                            ak = atom_indices[2]
                            ak_name = insert_prefix(mapping[ak][1].lower(), prefix)
                            fc = float(f"{parameter.k}".split()[0])
                            angle = float(f"{parameter.angle}".split()[0])

                            outfile.write(f"{ai_name:13}{aj_name:13}{ak_name:13}{fc:>10.2f}{angle:>11.3f}\n")

                    if block == 4:
                        for atom_indices, parameter in parameters["ProperTorsions"].items():
                            forces = []
                            ai = atom_indices[0]
                            ai_name = insert_prefix(mapping[ai][1].lower(), prefix)
                            aj = atom_indices[1]
                            aj_name = insert_prefix(mapping[aj][1].lower(), prefix)
                            ak = atom_indices[2]
                            ak_name = insert_prefix(mapping[ak][1].lower(), prefix)
                            al = atom_indices[3]
                            al_name = insert_prefix(mapping[al][1].lower(), prefix)
                            max_phase = len(parameter.phase)

                            # Now check if there are multiple minima
                            for i in range(0, max_phase):
                                fc = float(f"{parameter.k[i]}".split()[0])
                                phase = float(f"{parameter.phase[i]}".split()[0])
                                paths = int(parameter.idivf[i])

                                if i != max_phase - 1 and max_phase > 1:
                                    minimum = float(parameter.periodicity[i]) * -1

                                else:
                                    minimum = float(parameter.periodicity[i])

                                force = (fc, minimum, phase, paths)
                                forces.append(force)

                            for force in forces:
                                outfile.write(
                                    f"""{ai_name:13}{aj_name:13}{ak_name:13}{al_name:13}{force[0]:>10.4f}{force[1]:>6.1f}{force[2]:>11.1f}{force[3]:>6.1f}\n"""
                                )

                    if block == 5:
                        for atom_indices, parameter in parameters["ImproperTorsions"].items():
                            ai = atom_indices[0]
                            ai_name = insert_prefix(mapping[ai][1].lower(), prefix)
                            aj = atom_indices[1]
                            aj_name = insert_prefix(mapping[aj][1].lower(), prefix)
                            ak = atom_indices[2]
                            ak_name = insert_prefix(mapping[ak][1].lower(), prefix)
                            al = atom_indices[3]
                            al_name = insert_prefix(mapping[al][1].lower(), prefix)
                            fc = float(f"{parameter.k[0]}".split()[0])
                            phase = float(f"{parameter.phase[0]}".split()[0])
                            outfile.write(
                                f"""{ai_name:13}{aj_name:13}{ak_name:13}{al_name:13}{fc:10.1f}{phase:11.1f}\n"""
                            )
        finally:
            if should_close:
                outfile.close()

    def write_PDB(self, lname: str, outfile: Optional[TextIO] = None, residue_name: str = "LIG"):
        """Writes pdb file for a given ligand.

        Args:
            lname: name of the ligand to be written.
            outfile: file to write the .pdb file to. Defaults to None, writing self.out_dir / lname.pdb.
            prefix: a custom prefix to be added to the atom names in Q's .prm file. Element symbols will be
                lowercase not to be confused with protein atoms. Prefix should be used if more than one
                ligand will be present in the system. E.g.: "X", used as standard for ligand 2 in QligFEP.
                Defaults to None.
            residue_name: name of the residue in the .prm file, matching the residue on the pdb file.
                Defaults to "LIG".
        """
        mapping = self.mapping[lname]

        if outfile is None:
            pdbfile_out = self.out_dir / f"{lname}.pdb"
            outfile = open(pdbfile_out, "w")  # noqa: SIM115
            should_close = True
        else:
            should_close = False

        try:
            for atom in mapping:
                ai = atom + 1
                ai_name = mapping[atom][1]
                a_el = mapping[atom][2]
                ax = float(mapping[atom][4])
                ay = float(mapping[atom][5])
                az = float(mapping[atom][6])
                at_entry = [
                    "HETATM",  #  0 ATOM/HETATM
                    ai,  #  1 ATOM serial number
                    ai_name,  #  2 ATOM name
                    "",  #  3 Alternate location indicator
                    residue_name,  #  4 Residue name
                    "",  #  5 Chain identifier
                    1,  #  6 Residue sequence number
                    "",  #  7 Code for insertion of residue
                    ax,  #  8 Orthogonal coordinates for X
                    ay,  #  9 Orthogonal coordinates for Y
                    az,  # 10 Orthogonal coordinates for Z
                    0.0,  # 11 Occupancy
                    0.0,  # 12 Temperature factor
                    a_el,  # 13 Element symbol
                    "",  # 14 Charge on atom
                ]
                outfile.write(pdb_parse_out(at_entry) + "\n")
        finally:
            if should_close:
                outfile.close()

    def write_cofactor_plus_ff_files(self, ff: str):
        """Method for writing the .lib, .prm and .pdb files for the cofactor and the forcefield."""
        ff_lib, ff_prm = get_force_field_paths(ff)
        protein_lib_lines = Path(ff_lib).read_text().splitlines()
        protein_prm_lines = Path(ff_prm).read_text().splitlines()
        protein_prm_sections = parse_prm(protein_prm_lines)
        protein_prm_header = [line for line in protein_prm_lines[:15] if line.startswith("*")]

        prefixes = list("QWERTY")
        residues = ["CFA", "CFB", "CFC", "CFD", "CFE", "CFF"]  # cofactor a, b, ... f

        # raise an error if there are more cofactors than prefixes available
        if len(self.lig_names) > len(prefixes):
            raise ValueError(
                "More cofactors than prefixes available. Are you sure you have that many "
                "cofactors? If so, add extra prefixes & to the residues list."
            )

        lig_prm_contents = {}
        for name, prefix, res in zip(self.lig_names, prefixes, residues):
            lib_out = StringIO()
            self.write_lib_Q(name, outfile=lib_out, prefix=prefix, residue_name=res)
            lib_lines = lib_out.getvalue().split("\n")
            protein_lib_lines.extend(lib_lines)
            protein_lib_lines.append("*" + ("-" * 80))
            lib_out.close()

            prm_out = StringIO()
            self.write_prm_Q(name, outfile=prm_out, prefix=prefix)
            clean_output = []
            for line in prm_out.getvalue().split("\n"):
                if not line.strip():
                    continue
                if line.startswith("!"):
                    continue
                clean_output.append(line)
            sections = parse_prm(clean_output)
            lig_prm_contents[name] = sections
            prm_out.close()

            pdb_out = StringIO()  # get the pdb file to output as cofactors.pdb
            self.write_PDB(name, outfile=pdb_out, residue_name=res)
            if prefix == prefixes[0]:
                cofactor = read_pdb_to_dataframe(pdb_out.getvalue().split("\n"))
            else:
                cofactor = append_pdb_to_another(
                    cofactor, read_pdb_to_dataframe(pdb_out.getvalue().split("\n"))
                )
            pdb_out.close()

        final_prm_contents = [*protein_prm_header]
        for header, lines in protein_prm_sections.items():
            final_prm_contents.append(f"[{header}]")
            final_prm_contents.append(protein_prm_sections[header])

            if header != "options":
                final_prm_contents.append(f"! Cofactor {header} parameters")

            for name in self.lig_names:
                try:
                    cofactor_contents = lig_prm_contents[name][header]
                except KeyError:
                    logger.info(f'No "{header}" parameters found for cofactor {name}')
                    continue
                if cofactor_contents:
                    final_prm_contents.append(cofactor_contents)

            if header != "options":
                final_prm_contents.append(f"! End cofactor {header} parameters")
            final_prm_contents.append("")

        # Write the final files
        lib_out = self.out_dir / f"{ff}_plus_cofactor.lib"
        prm_out = self.out_dir / f"{ff}_plus_cofactor.prm"

        with open(lib_out, "w") as f:
            f.write("\n".join(protein_lib_lines))
        with open(prm_out, "w") as f:
            f.write("\n".join(final_prm_contents))
        write_dataframe_to_pdb(cofactor, self.out_dir / "all_cofactors.pdb")
