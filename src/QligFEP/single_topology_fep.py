"""Single-topology single-Hamiltonian FEP builder.

A single-topology hybrid represents the shared core of two congeneric ligands once
(those atoms morph A->B) and only grows/shrinks the unique atoms, so it inherits the
small-difference convergence of a multistate edge while keeping the soft-core acting
on the appearing/vanishing atoms. `SingleTopologyFEP` subclasses `FEP` and overrides
only the topology-specific build steps; the shared machinery (lambda schedule, prm
merge, qprep run, counter-ion handling, ...) is reused from the base.

The hybrid = all of ligand A (shared core + vanishing atoms) plus ligand B's unique
(appearing) atoms. Shared atoms morph A->B (charge interpolates; LJ type morphs to the
X-prefixed ligand-B type); vanishing atoms go real->DUM; appearing atoms go DUM->real.
Ligand B's atom types clash with ligand A's (both use c1,c2,...), so every ligand-B
type is X-prefixed exactly as the dual-topology builder does -- only here the shared
atoms are collapsed to one copy.
"""
from pathlib import Path

from rdkit import Chem

from .pdb_utils import read_pdb_to_dataframe
from .qligfep import FEP
from .single_topology import map_single_topology

#: Residue name of the single hybrid ligand (QligFEP convention, kept literal).
HYBRID_RESNAME = "{LIG}"


def _parse_lib(path):
    """Return atoms=[(name, type, charge)], bonds=[(n1, n2)], impropers=[(n1..n4)]."""
    atoms, bonds, impropers, block = [], [], [], None
    for line in Path(path).read_text().splitlines():
        s = line.split()
        if not s:
            continue
        if s[0].startswith("["):
            block = s[0]
            continue
        if block == "[atoms]" and s[0].isdigit():
            atoms.append((s[1], s[2], float(s[3])))
        elif block == "[bonds]":
            bonds.append((s[0], s[1]))
        elif block == "[impropers]":
            impropers.append(tuple(s[:4]))
    return atoms, bonds, impropers


def _parse_pdb_coords(path):
    """name -> (x, y, z) from ATOM/HETATM records."""
    coords = {}
    for line in Path(path).read_text().splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            coords[line[12:16].strip()] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            )
    return coords


def _element(name):
    """Element symbol = the alphabetic part of an atom name (C1 -> C, Cl11 -> Cl)."""
    return "".join(ch for ch in name if ch.isalpha())


class SingleTopologyFEP(FEP):
    """Build a single-topology hybrid for single-Hamiltonian FEP of two ligands."""

    def read_files(self):
        """Map the two ligands and build the hybrid atom table.

        Returns the same nested shape the orchestration expects:
        ([_, type_replacements], [change_charges, change_atoms], [n_hybrid, 0]).
        The X-prefix `type_replacements` feeds the inherited change_prm; the
        per-atom change_charges/change_atoms and the hybrid table feed
        write_FEP_file. The whole hybrid counts as one restraint group, so the
        second "ligand size" is 0.
        """
        atoms_a, bonds_a, impropers_a = _parse_lib(self.lig1 + ".lib")
        atoms_b, bonds_b, impropers_b = _parse_lib(self.lig2 + ".lib")
        coords_a = _parse_pdb_coords(self.lig1 + ".pdb")
        coords_b = _parse_pdb_coords(self.lig2 + ".pdb")

        mol_a = Chem.MolFromMolFile(self.lig1 + ".sdf", removeHs=False)
        mol_b = Chem.MolFromMolFile(self.lig2 + ".sdf", removeHs=False)
        if len(atoms_a) != mol_a.GetNumAtoms() or len(atoms_b) != mol_b.GetNumAtoms():
            raise ValueError("lib atom count does not match sdf; lib order must equal sdf order")

        self.st_map = map_single_topology(mol_a, mol_b, overlap_thresh=0.6)
        a_to_b = self.st_map.a_to_b
        shared_a = set(a_to_b)
        shared_b = set(a_to_b.values())

        # ligand-B atom index -> hybrid name: shared keep their A-partner's name,
        # appearing atoms are X-prefixed.
        name_b = {j: "X" + atoms_b[j][0] for j in self.st_map.appear}
        for i_a, j_b in a_to_b.items():
            name_b[j_b] = atoms_a[i_a][0]

        hyb = []
        for i, (name, ljtype, charge) in enumerate(atoms_a):
            if i in shared_a:
                j = a_to_b[i]
                state_b = ("X" + atoms_b[j][1], atoms_b[j][2])
                role = "shared"
            else:
                state_b = ("DUM", 0.0)
                role = "vanish"
            hyb.append({
                "name": name, "type": ljtype, "charge": charge, "xyz": coords_a[name],
                "role": role, "sA": (ljtype, charge), "sB": state_b, "elem": _element(name),
            })
        for j in self.st_map.appear:
            state_b = ("X" + atoms_b[j][1], atoms_b[j][2])
            hyb.append({
                "name": name_b[j], "type": state_b[0], "charge": atoms_b[j][2],
                "xyz": coords_b[atoms_b[j][0]], "role": "appear",
                "sA": ("DUM", 0.0), "sB": state_b, "elem": _element(atoms_b[j][0]),
            })

        self.hyb = hyb
        self.hyb_index = {h["name"]: k + 1 for k, h in enumerate(hyb)}
        self.bonds = self._merge_terms(bonds_a, bonds_b, atoms_b, name_b, shared_b)
        self.impropers = self._merge_terms(impropers_a, impropers_b, atoms_b, name_b, shared_b)

        type_replacements = {ljtype: "X" + ljtype for (_, ljtype, _) in atoms_b}
        change_charges = [[k + 1, h["sA"][1], h["sB"][1]] for k, h in enumerate(hyb)]
        change_atoms = [[k + 1, h["sA"][0], h["sB"][0]] for k, h in enumerate(hyb)]

        self.charge_lig1 = round(sum(h["sA"][1] for h in hyb))
        self.charge_lig2 = round(sum(h["sB"][1] for h in hyb))
        self.same_charge = self.charge_lig1 == self.charge_lig2

        return ([type_replacements, type_replacements],
                [change_charges, change_atoms],
                [len(hyb), 0])

    @staticmethod
    def _merge_terms(terms_a, terms_b, atoms_b, name_b, shared_b):
        """Bonds/impropers of ligand A as-is, plus ligand-B terms that touch at
        least one appearing atom (mapped to hybrid names)."""
        b_names = [a[0] for a in atoms_b]
        b_name_to_hyb = {atoms_b[j][0]: name_b[j] for j in range(len(atoms_b))}
        merged = list(terms_a)
        for term in terms_b:
            indices = [b_names.index(n) for n in term]
            if any(j not in shared_b for j in indices):
                merged.append(tuple(b_name_to_hyb[n] for n in term))
        return merged

    def change_lib(self, replacements, writedir):
        """Write the single hybrid library (replaces the dual two-library scheme)."""
        with open(writedir + "/singletop.lib", "w") as f:
            f.write(HYBRID_RESNAME + "\n\n[info]\n SYBYLtype RESIDUE\n\n[atoms]\n")
            for k, h in enumerate(self.hyb, 1):
                f.write(f"{k:5d}   {h['name']:<9s} {h['type']:<12s} {h['charge']:8.4f}\n")
            f.write("\n[bonds]\n")
            for x, y in self.bonds:
                f.write(f"{x:<9s} {y:<9s}\n")
            f.write("\n[impropers]\n")
            for improper in self.impropers:
                f.write("  ".join(f"{n:<8s}" for n in improper) + "\n")

    def _fep_atom_offset(self, writedir):
        """FEP-file topology offset. The water leg's only solute is the hybrid, so
        topology index == FEP index (offset 0). In the protein leg the hybrid is the
        last residue, so the offset is every preceding atom (protein + solvent)."""
        if self.system == "water":
            return 0
        total = read_pdb_to_dataframe(Path(writedir) / "top_p.pdb").shape[0]
        return total - len(self.hyb)

    def write_FEP_file(
        self, change_charges, change_atoms, FEP_vdw, writedir, lig_size1, lig_size2,
        softcore_method="gapsys", single_hamiltonian=True,
    ):
        """Write the single-Hamiltonian FEP file: per-atom charge/type morph with
        DUM-conditional soft-core (0 on real ends, 20 on a ghost end)."""
        offset = self._fep_atom_offset(writedir)
        with open(writedir + "/FEP1.fep", "w") as f:
            f.write(f"!info: {self.lig1} --> {self.lig2} (single topology)\n")
            f.write("[FEP]\nstates 2\nsoftcore_use_max_potential on\n")
            f.write(f"softcore_method {softcore_method}\n")
            if single_hamiltonian:
                f.write("single_hamiltonian on\n")

            f.write("\n[atoms]\n")
            for k in range(1, len(self.hyb) + 1):
                f.write(f"{k:<5d}{k + offset}\n")

            f.write("\n[change_charges]\n")
            for idx, charge_a, charge_b in change_charges:
                f.write(f"{idx:<5d}{charge_a:9.4f}{charge_b:9.4f}\n")

            f.write("\n[atom_types]\n")
            f.write("DUM       0.0000    0.0000    0    0    0.0000    0.0000    1.0080\n")
            for line in FEP_vdw:
                f.write(line + "\n")

            f.write("\n[change_atoms]\n")
            for idx, type_a, type_b in change_atoms:
                f.write(f"{idx:<5d}{type_a:<9s}{type_b:<9s}\n")

            f.write("\n[softcore]\n")
            for k, h in enumerate(self.hyb, 1):
                sc_a = 20 if h["sA"][0] == "DUM" else 0
                sc_b = 20 if h["sB"][0] == "DUM" else 0
                f.write(f"{k:<5d}{sc_a:<6d}{sc_b}\n")
