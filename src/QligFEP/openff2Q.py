"""Module containing the OpenFF2Q class to process ligands and generate OpenFF parameter files for QligFEP."""

import numpy as np
from joblib import Parallel, delayed, parallel_config
from openff.toolkit import ForceField, Molecule, Topology
from tqdm import tqdm

from .chemIO import MoleculeIO
from .logger import logger
from .pdb_utils import pdb_parse_out
from .settings.settings import FF_DIR


class OpenFF2Q(MoleculeIO):
    """Class to process ligands and generate OpenFF parameter files for QligFEP. Dictionary
    variables use ligand names as keys, as setup in the MoleculeIO class.

    Attributes:
        mapping: Dictionary to store the mapping of the ligand atoms to the forcefield parameters.
        forcefield: OpenFF ForceField object.
        topologies: Dictionary with `openff.toolkit.topology.Topology` objects as values.
        parameters: Dictionary with `openff.toolkit.typing.engines.smirnoff.parameters.ParameterList`
            objects as values.
        charges_list_magnitude: Dictionary to store the partial charge magnitude for each
            atom in the ligand for each ligand.
        total_charges: Dictionary to store the total charges for each ligand.
    """

    def __init__(self, lig, pattern="*.sdf", n_jobs=1):
        """Initializes a new instance of OpenFF2Q to process the `.sdf` input as lig.

        Args:
            lig: sdf file containing several molecules or directory containing the sdf files.
            pattern: If desired, a pattern can be used to search for sdf files within a directory with
                `glob`. If lig is a sdf file, this argument will be ignored. Defaults to None.
        """
        super().__init__(lig, pattern=pattern)
        self.n_jobs = n_jobs
        self.mapping = {lname: {} for lname in self.lig_names}
        self.forcefield = ForceField("openff-2.2.0.offxml")
        self.topologies, self.parameters = self.set_topologies_and_parameters()
        self.charges_list_magnitude = {}  # store charge magnitude for each ligand
        self.total_charges = {}  # store the total charges

    @staticmethod
    def _assign_charge(molecule: Molecule) -> np.ndarray:
        """Private method that assigns partial charges to an input Molecule and returns their magnitudes"""
        try:
            # Seems like the fastest way of doing this for now (if you're not OpenEye licensed);
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

    def process_ligands(self):
        """Assigns partial charges and writes the .lib, .prm and .pdb files for each ligand."""
        logger.info("Calculating charges")
        with parallel_config(
            n_jobs=self.n_jobs, backend="multiprocessing"
        ):  # backend="threading"
            molecules = tqdm(self.molecules)
            charges_magnitudes = Parallel()(
                delayed(self._assign_charge)(molecule) for molecule in molecules
            )
        logger.info("Done! Writing .lib, .prm and .pdb files for each ligand")
        for lname, charges in zip(self.lig_names, charges_magnitudes):
            self.charges_list_magnitude.update({lname: charges})
            formatted_sum = f'{round(charges_magnitudes.sum(), 10):.3f}'
            if formatted_sum == '-0.000':
                formatted_sum = '0.000'
            self.total_charges.update({lname: formatted_sum})
            self.get_mapping(lname)
            self.write_lib_Q(lname)
            self.write_prm_Q(lname)
            self.write_PDB(lname)
        all_formal_charges = [self.total_charges[n] for n in self.lig_names]
        if np.unique(all_formal_charges).size > 1:
            logger.warning(
                f"Formal charges of ligands in .sdf are not unique: {self.total_charges}"
            )

    def get_mapping(self, lname):
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

    def write_lib_Q(self, lname: str):
        """Writes Q's .lib file for a given ligand.

        Args:
            lname: name of the ligand for the .lib file.
        """
        parameters = self.parameters[lname]
        mapping = self.mapping[lname]
        total_charge = self.total_charges[lname]
        with open(lname + ".lib", "w") as outfile:
            outfile.write(
                "{}    ! atoms no {}   total charge {} \n\n".format(
                    "{LIG}", len(mapping), total_charge
                )
            )

            outfile.write("[info] \n SYBYLtype RESIDUE \n\n")

            # atom and charge block:
            outfile.write("[atoms] \n")
            for i, at in enumerate(mapping):
                outfile.write(
                    "{:>4s}   {:10}{:11}{:>10s}\n".format(
                        mapping[at][0],
                        mapping[at][1],
                        mapping[at][1].lower(),
                        mapping[at][3],
                    )
                )
            # bonded block
            outfile.write("\n[bonds]\n")
            for i, bond in enumerate(parameters["Bonds"]):
                ai = mapping[bond[0]][1]
                aj = mapping[bond[1]][1]
                outfile.write("{:10s}{:}\n".format(ai, aj))

            # improper block
            outfile.write("\n[impropers]\n")
            for i, torsion in enumerate(parameters["ImproperTorsions"]):
                ai = mapping[torsion[0]][1]
                aj = mapping[torsion[1]][1]
                ak = mapping[torsion[2]][1]
                al = mapping[torsion[3]][1]
                outfile.write("{:10}{:10}{:10}{}\n".format(ai, aj, ak, al))

    def write_prm_Q(self, lname: str):
        """Writes Q's .prm file for a given ligand.

        Args:
            lname: name of the ligand for the .prm file.
        """
        prm_file = str(FF_DIR / "NOMERGE.prm")
        parameters = self.parameters[lname]
        mapping = self.mapping[lname]
        prm_file_out = f"{lname}.prm"
        mol = self.molecules[self.lig_names.index(lname)]
        with open(prm_file) as infile, open(prm_file_out, "w") as outfile:
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
                        ai_name = mapping[ai][1].lower()
                        # This is a bit hacky, check how to get the float out directly
                        epsilon = float("{}".format(parameter.epsilon).split()[0])
                        epsilon23 = epsilon / 2
                        # TO DO: CHECK IF THIS IS CORRECT!
                        Rmin = "{}".format(parameter.rmin_half)
                        Rmin = Rmin.split()[0]
                        Rmin = float(Rmin)
                        assert (
                            len(atom_indices) == 1
                        ), f"More than 1 atom indices present: {atom_indices}"
                        mass = str(round(mol.atoms[atom_indices[0]].mass.magnitude, 4))
                        outfile.write(
                            """{:6}{: 8.3f}{: 10.3f}{: 10.3f}{: 10.3f}{: 10.3f}{:>10s}\n""".format(
                                ai_name, Rmin, 0.00, epsilon, Rmin, epsilon23, mass
                            )
                        )

                if block == 2:
                    for atom_indices, parameter in parameters["Bonds"].items():
                        ai = atom_indices[0]
                        ai_name = mapping[ai][1].lower()
                        aj = atom_indices[1]
                        aj_name = mapping[aj][1].lower()
                        fc = float("{}".format(parameter.k).split()[0])
                        l = float("{}".format(parameter.length).split()[0])
                        outfile.write(
                            "{:10}{:10}{:10.1f}{:>10.3f}\n".format(
                                ai_name, aj_name, fc, l
                            )
                        )

                if block == 3:
                    for atom_indices, parameter in parameters["Angles"].items():
                        ai = atom_indices[0]
                        ai_name = mapping[ai][1].lower()
                        aj = atom_indices[1]
                        aj_name = mapping[aj][1].lower()
                        ak = atom_indices[2]
                        ak_name = mapping[ak][1].lower()
                        fc = float("{}".format(parameter.k).split()[0])
                        angle = float("{}".format(parameter.angle).split()[0])

                        outfile.write(
                            """{:10}{:10}{:10}{: 8.2f}{:>12.3f}\n""".format(
                                ai_name, aj_name, ak_name, fc, angle
                            )
                        )

                if block == 4:
                    for atom_indices, parameter in parameters["ProperTorsions"].items():
                        forces = []
                        ai = atom_indices[0]
                        ai_name = mapping[ai][1].lower()
                        aj = atom_indices[1]
                        aj_name = mapping[aj][1].lower()
                        ak = atom_indices[2]
                        ak_name = mapping[ak][1].lower()
                        al = atom_indices[3]
                        al_name = mapping[al][1].lower()
                        max_phase = len(parameter.phase)

                        # Now check if there are multiple minima
                        for i in range(0, max_phase):
                            fc = float("{}".format(parameter.k[i]).split()[0])
                            phase = float("{}".format(parameter.phase[i]).split()[0])
                            paths = int(parameter.idivf[i])

                            if i != max_phase - 1 and max_phase > 1:
                                minimum = float(parameter.periodicity[i]) * -1

                            else:
                                minimum = float(parameter.periodicity[i])

                            force = (fc, minimum, phase, paths)
                            forces.append(force)

                        for force in forces:
                            outfile.write(
                                """{:10}{:10}{:10}{:10}{:>10.3f}{:>10.3f}{:>10.3f}{:>5d}\n""".format(
                                    ai_name,
                                    aj_name,
                                    ak_name,
                                    al_name,
                                    force[0],
                                    force[1],
                                    force[2],
                                    force[3],
                                )
                            )

                if block == 5:
                    for atom_indices, parameter in parameters["ImproperTorsions"].items():
                        ai = atom_indices[0]
                        ai_name = mapping[ai][1].lower()
                        aj = atom_indices[1]
                        aj_name = mapping[aj][1].lower()
                        ak = atom_indices[2]
                        ak_name = mapping[ak][1].lower()
                        al = atom_indices[3]
                        al_name = mapping[al][1].lower()
                        fc = float("{}".format(parameter.k[0]).split()[0])
                        phase = float("{}".format(parameter.phase[0]).split()[0])
                        outfile.write(
                            """{:10}{:10}{:10}{:10}{:10.3f}{:10.3f}\n""".format(
                                ai_name, aj_name, ak_name, al_name, fc, phase
                            )
                        )

    def write_PDB(self, lname: str):
        """Writes a PDB file for a given ligand.

        Args:
            lname: name of the ligand for the .pdb file.
        """
        mapping = self.mapping[lname]
        with open(lname + ".pdb", "w") as outfile:
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
                    "LIG",  #  4 Residue name
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

    # def report_missing_parameters(self): # TODO: This doesn't work yet...
    #     """
    #     Analyze a molecule using a provided ForceField, generating a report of any
    #     chemical groups in the molecule that are lacking parameters.

    #     Parameters
    #     ----------
    #     molecule : an openforcefield.topology.FrozenMolecule
    #         The molecule to analyze
    #     forcefield : an openforcefield.typing.engine.smirnoff.ForceField
    #         The ForceField object to use

    #     Returns
    #     -------
    #     missing_parameters : dict[tagname: list[dict[tagged_smiles:string, image:PIL.Image, atom indices:list[int]]]]
    #         A hierarchical dictionary, with first level keys indicating ForceField tag
    #         names (eg. "Bonds"), and first-level values which are lists of dictionaries.
    #         Each dictionary in this list reflects one missing parameter, and contains the
    #         following key:value pairs :
    #         * "image": PIL.Image
    #             * shows a 2D drawing, highlighting the feature that could not be parametrized
    #         * "tagged_smiles": string
    #             * SMILES of the whole molecule, tagging the atom indices which could not be
    #               parametrized
    #         * "atom_indices": tuple(int)
    #             * The indices of atoms which could not be parametrized

    #     """
    #     highlight_color = (0.75, 0.75, 0.75)

    #     # Make deepcopies of both inputs, since we may modify them in this function
    #     forcefield = deepcopy(self.forcefield)

    #     ###### lig.sdf file ######
    #     molecule = deepcopy(self.molecule)

    #     # Set partial charges to placeholder values so that we can skip AM1-BCC
    #     # during parameterization

    #     # the partial charges have to be added to the molecule somehow with an array like thing multiplied with the elementary charge
    #     # They have to be the same format. It works for the molecule.partial_charges = (np.zeros(molecule.n_atoms) + 0.1) * unit.elementary_charge
    #     # molecule.n_atoms just gives 45 (number of atoms). So that should not be
    #     molecule.partial_charges = np.array(self.charges_list_magnitude) * unit.elementary_charge

    #     # Prepare dictionary to catch parameterization failure info
    #     success = False
    #     missing_params = {}

    #     while not success:
    #         # Try to parameterize the system, catching the exception if there is one.
    #         try:

    #             ###### error occurs, Molecule Brc1cccc(Nc2nc(OCC3CCCCC3)c3nc[nH]c3n2)c1 has a net charge of 4.5 ######
    #                 #### why do you calculate charges when you have them already? ####
    #                 #### why molecule.to_topology() when we have that already in the openff function? ####
    #                 #### charge_from_molecules (List[openff.toolkit.molecule.Molecule], optional) – If specified, partial #
    #                 #    charges will be taken from the given molecules instead of being determined by the force field. #####

    #             forcefield.create_openmm_system(molecule.to_topology(), charge_from_molecules=[molecule])
    #             success = True

    #         ###### NameError: name 'UnassignedValenceParameterException' is not defined ######
    #         except UnassignedValenceParameterException as e:
    #             success = False

    #             # Ensure that there is a list initialized for missing parameters
    #             # under this tagname
    #             handler_tagname = e.handler_class._TAGNAME
    #             if handler_tagname not in missing_params:
    #                 missing_params[handler_tagname] = []

    #             # Create a shortcut to the topology atom tuples attached to
    #             # the parametrization error
    #             top_atom_tuples =  e.unassigned_topology_atom_tuples

    #             # Make a summary of the missing parameters from this attempt and add it to
    #             # the missing_params dict
    #             rdmol = molecule.to_rdkit()
    #             for top_atom_tuple in top_atom_tuples:
    #                 orig_atom_indices = [i.topology_atom_index for i in top_atom_tuple]
    #                 # Make a copy of the input RDMol so that we don't modify the original
    #                 this_rdmol = deepcopy(rdmol)

    #                 # Attach tags to relevant atoms so that a tagged SMILES can be written
    #                 orig_rdatoms = []
    #                 for tag_idx, atom_idx in enumerate(orig_atom_indices):
    #                     rdatom = this_rdmol.GetAtomWithIdx(atom_idx)
    #                     rdatom.SetAtomMapNum(tag_idx + 1)
    #                     orig_rdatoms.append(rdatom)

    #                 tagged_smiles = Chem.MolToSmiles(this_rdmol)

    #                 # Make tagged hydrogens into deuteriums so that RemoveHs doesn't get rid of them
    #                 for rdatom in orig_rdatoms:
    #                     if rdatom.GetAtomicNum() == 1:
    #                         rdatom.SetIsotope(2)

    #                 # Remove hydrogens, since they clutter up the 2D drawing
    #                 # (tagged Hs are not removed, since they were converted to deuterium)
    #                 h_less_rdmol = Chem.RemoveHs(this_rdmol)

    #                 # Generate 2D coords, since drawing from 3D can look really weird
    #                 Draw.rdDepictor.Compute2DCoords(h_less_rdmol)

    #                 # Search over the molecule to find the indices of the tagged atoms
    #                 # after hydrogen removal
    #                 h_less_atom_indices = [None for i in orig_atom_indices]
    #                 for rdatom in h_less_rdmol.GetAtoms():
    #                     # Convert deuteriums back into hydrogens
    #                     if rdatom.GetAtomicNum() == 1:
    #                         rdatom.SetIsotope(1)

    #                     atom_map_num = rdatom.GetAtomMapNum()
    #                     if atom_map_num == 0:
    #                         continue
    #                     h_less_atom_indices[atom_map_num-1] = rdatom.GetIdx()

    #                 # Once the new atom indices are found, use them to find the H-less
    #                 # bond indices
    #                 h_less_rdbonds = []
    #                 for i in range(len(h_less_atom_indices)-1):
    #                     rdbond = h_less_rdmol.GetBondBetweenAtoms(
    #                                             h_less_atom_indices[i],
    #                                             h_less_atom_indices[i+1])
    #                     h_less_rdbonds.append(rdbond)
    #                 h_less_bond_indices = [bd.GetIdx() for bd in h_less_rdbonds]

    #                 # Create a 2D drawing of the molecule, highlighting the
    #                 # parameterization failure
    #                 highlight_atom_colors = {idx:highlight_color for idx in h_less_atom_indices}
    #                 highlight_bond_colors = {idx:highlight_color for idx in h_less_bond_indices}
    #                 image = Draw.MolsToGridImage([h_less_rdmol],
    #                                              highlightAtomLists=[h_less_atom_indices],
    #                                              highlightBondLists=[h_less_bond_indices],
    #                                              molsPerRow=1,
    #                                              highlightAtomColors=[highlight_atom_colors],
    #                                              highlightBondColors=[highlight_bond_colors],
    #                                              subImgSize=(600,600)
    #                                             )

    #                 # Structure and append the relevant info to the missing_params dictionary
    #                 param_description = {'atom_indices': orig_atom_indices,
    #                                      'image': image,
    #                                      'tagged_smiles': tagged_smiles
    #                                     }
    #                 missing_params[handler_tagname].append(param_description)

    #             # Add a "super generic" parameter to the top of this handler's ParameterList,
    #             # which will make it always find parameters for each term. This will prevent the same
    #             # parameterization exception from being raised in the next attempt.
    #             param_list = forcefield.get_parameter_handler(handler_tagname).parameters
    #             param_list.insert(0, super_generics[handler_tagname])

    #     if success is not True:
    #         print(missing_params)
    #     else:
    #         print('Parameters succesfully assigned')
