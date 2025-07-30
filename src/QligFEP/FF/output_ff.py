import math
import re

import requests
import simtk
import xmltodict
from openmm.app import element


def array_calc(atom, sig, eps, mass=None):
    # calculate LJ types
    sig = float(sig)
    eps = float(eps)
    if mass is None:
        mass = 0.0
    return f"{atom:<13}{sig:>10}{0.0:>11}" f"{eps:>11}{sig:>11}{eps/2:>11.6f}{mass:>11.2f}\n"


def atom_type_forces_from_openmm(sig, eps):
    # convert eps from kJ/mol to kcal/mol
    sig = round(sig * 5.612, 4)
    # convert sig to A
    eps = round(eps * 23.9005736 / 100, 4)
    return sig, eps


def convert_kj_mol_nm2_to_kcal_mol_a2(value_kj_mol_nm2):
    # 100 joules = 23.9005736 cal
    conversion_factor = 23.9005736
    value_kcal_mol_a2 = round(value_kj_mol_nm2 * conversion_factor / 10000, 1)
    return value_kcal_mol_a2


def radians_to_degrees(radians):
    return radians * 180 / math.pi


def extract_atom_types(atom_types: dict):
    extracted = []
    for atom in atom_types:
        atom_name = atom["@type"].replace("protein-", "")
        sig, eps = atom_type_forces_from_openmm(float(atom["@sigma"]), float(atom["@epsilon"]))
        mass = float(atom_type_to_mass[atom_name])
        string = array_calc(atom_name, sig, eps, mass)
        extracted.append(string)
    return extracted


def extract_bonds(bonds: dict):
    extracted = []
    for bond in bonds:
        _type1 = bond["@type1"].replace("protein-", "")
        _type2 = bond["@type2"].replace("protein-", "")
        k = float(bond["@k"])
        converted_k = convert_kj_mol_nm2_to_kcal_mol_a2(k)
        converted_len = round(float(bond["@length"]) * 10, 4)
        string = f"{_type1:<13}{_type2:<13}{converted_k:>10}{converted_len:>11}\n"
        all_lines.append(string)
    return extracted


def extract_angles(angles: dict):
    extracted = []
    for angle in angles:
        _type1 = angle["@type1"].replace("protein-", "")
        _type2 = angle["@type2"].replace("protein-", "")
        _type3 = angle["@type3"].replace("protein-", "")
        k = float(angle["@k"])
        converted_k = round(convert_kj_mol_nm2_to_kcal_mol_a2(k) * 100, 1)
        converted_theta = round(radians_to_degrees(float(angle["@angle"])), 1)
        string = f"{_type1:<13}{_type2:<13}{_type3:<13}{converted_k:>10}{converted_theta:>11}\n"
        extracted.append(string)
    return extracted


def extract_torsions(torsions: dict):
    extracted = []
    for torsion in torsions:
        types = [
            torsion[f"@type{idx}"].replace("protein-", "") if torsion[f"@type{idx}"] != "" else "?"
            for idx in range(1, 5)
        ]
        periodicities = []
        phases = []
        k_values = []
        for idx in range(1, 5):
            try:
                periodicities.append(int(torsion[f"@periodicity{idx}"]))
                phases.append(float(torsion[f"@phase{idx}"]))
                k_values.append(float(torsion[f"@k{idx}"]))
            except KeyError:
                break
        # FIXME: the number of paths is currently wrong... See definition from Q manual:
        # The number of paths is the number of ways that a two-atom torsion can
        # be defined, i:e: the product of the number of atoms bonded to the two middle atoms
        # I opened an issue on that on Q's repository: https://github.com/qusers/Q6/issues/20
        # TODO: try to figure it out; this piece of code might help with answering it:
        # https://github.com/mpurg/qtools/blob/master/qscripts-cli/q_amber2q.py#L179-L183
        num_paths = len(periodicities)

        for idx, period in enumerate(periodicities):
            phase = round(radians_to_degrees(phases[idx]), 1)
            k = round(k_values[idx] / 4.184, 1)
            # check if more components to follow:
            if period != periodicities[-1]:
                period = -period
            fmt_string = (
                f"{types[0]:<13}{types[1]:<13}{types[2]:<13}{types[3]:<13}"
                f"{k:>10.2f}{period:>6.1f}{phase:>11.1f}{num_paths:>6.1f}\n"
            )
            extracted.append(fmt_string)
    return extracted


phosaa14SB_path = "https://raw.githubusercontent.com/openmm/openmmforcefields/main/openmmforcefields/ffxml/amber/phosaa14SB.xml"
amber14SB_path = "https://raw.githubusercontent.com/openmm/openmm/8.1.1/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml"

phoosa_content = requests.get(phosaa14SB_path).text
amber_content = requests.get(amber14SB_path).text

phoosa_dict = xmltodict.parse(phoosa_content)
amber_dict = xmltodict.parse(amber_content)

atom_type_to_mass = {}
p_atomtypes = phoosa_dict["ForceField"]["AtomTypes"]["Type"]
for atomtype in p_atomtypes:
    atom_type_to_mass[atomtype["@class"]] = atomtype["@mass"]
a_atomtypes = amber_dict["ForceField"]["AtomTypes"]["Type"]
for atomtype in a_atomtypes:
    atom_type_to_mass[atomtype["@class"]] = atomtype["@mass"]

aa_names = re.compile(
    "ALA|ARG|ASH|ASN|ASP|CYM|CYS|CYX|GLH|GLN|GLU|GLY|HID|HIE|HIP|HYP|ILE|LEU|LYN|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL"
)

print("\n\n\nPRINTING LIB FILES\n\n")

residues_dict = {}
full_string = []
residues = amber_dict["ForceField"]["Residues"]["Residue"]
for residue in residues:
    resname = residue["@name"]
    full_string.append(f"\n{{{resname}}}\n")
    full_string.append("    [atoms]\n")
    for idx, atom in enumerate(residue["Atom"], start=1):
        atom_name = atom["@name"]
        atom_type = atom["@type"].replace("protein-", "")
        atom_charge = float(atom["@charge"])
        try:
            atom_mass = atom_type_to_mass[atom_type]
        except KeyError:
            element.Element.getBySymbol(atom_type).mass.value_in_unit(simtk.unit.amu)
        full_string.append(f"{idx:>9}  {atom_name:<7}{atom_type:<14}{atom_charge:>9.6f}\n")
    full_string.append("    [bonds]\n")
    for bond in residue["Bond"]:
        atom1 = bond["@atomName1"]
        atom2 = bond["@atomName2"]
        full_string.append(f"        {atom1:<4}  {atom2:<4}\n")
    full_string.append("    [impropers]\n")
    full_string.append(
        "        -C    N     CA    H\n"  # present for all amino acids
        "        CA    C     +N    O\n"
    )
    full_string.append("    [connections]\n")
    external = residue["ExternalBond"]

    # In this case, n-terminus residue
    if aa_names.sub("", resname) == "N" or resname in ["NME", "NHE"]:
        full_string.append("        tail C\n")
    # In this case, c-terminus residue
    elif aa_names.sub("", resname) == "C" or resname == "ACE":
        full_string.append("        head N\n")
    elif resname == "CYX":
        full_string.append("        head N\n" "        tail C\n")
    else:
        assert (  # we make the disulfide directly on the Q topology
            len(external) == 2
        ), f"Something might be wrong! External bonds should be 2:\n{resname} - {external}"
        full_string.append("        head N\n" "        tail C\n")

residues = phoosa_dict["ForceField"]["Residues"]["Residue"]
for residue in residues:
    resname = residue["@name"]
    full_string.append(f"{{{resname}}}\n")
    full_string.append("    [atoms]\n")
    for idx, atom in enumerate(residue["Atom"], start=1):
        atom_name = atom["@name"]
        atom_type = atom["@type"]
        atom_charge = float(atom["@charge"])
        try:
            atom_mass = atom_type_to_mass[atom_type]
        except KeyError:
            element.Element.getBySymbol(atom_type).mass.value_in_unit(simtk.unit.amu)
        full_string.append(f"{idx:>9}  {atom_name:<7}{atom_type:<14}{atom_charge:>9.6f}\n")
    full_string.append("    [bonds]\n")
    for bond in residue["Bond"]:
        atom1 = bond["@atomName1"]
        atom2 = bond["@atomName2"]
        full_string.append(f"        {atom1:<4}  {atom2:<4}\n")
    full_string.append("    [impropers]\n")
    full_string.append(
        "        -C    N     CA    H\n"  # present for all amino acids
        "        CA    C     +N    O\n"
    )
    full_string.append("    [connections]\n")
    external = residue["ExternalBond"]
    assert len(external) == 2, "Something might be wrong! External bonds should be 2"
    full_string.append("        head N\n" "        tail C\n")
print("".join(full_string))


all_lines = []
header = (
    "*------------------------------------------------------------------------------------\n"
    "*\n"
    "* Q-FF parameters - ported from ff14SB & phosaa14SB in OpenMM by David Araripe (2024)\n"
    "*\n"
    "*------------------------------------------------------------------------------------\n"
    "[options]\n"
    "name                           Q-Amber14SB\n"
    "type                           AMBER\n"
    "vdw_rule                       arithmetic !vdW combination rule (geometric or arithmetic)\n"
    "scale_14                       0.8333 ! electrostatic 1-4 scaling factor\n"
    "switch_atoms                   off\n"
    "improper_potential             periodic\n"
    "improper_definition            explicit\n"
    "\n"
    "[atom_types]\n"
    "* tac               R*1        R*2   epsilon1        R*3 epsilon2&3       mass\n"
)
all_lines.append(header)


print("\n\n\nPRINTING PRM FILES\n\n")
atom_types = phoosa_dict["ForceField"]["NonbondedForce"]["Atom"]
all_lines.extend(extract_atom_types(atom_types))
atom_types = amber_dict["ForceField"]["NonbondedForce"]["Atom"]
all_lines.extend(extract_atom_types(atom_types))
all_lines.append("\n" "! Ligand vdW parameters\n" "! End ligand vdW parameters\n" "\n" "[bonds]\n")

bonds = phoosa_dict["ForceField"]["HarmonicBondForce"]["Bond"]
all_lines.extend(extract_bonds(bonds))
bonds = amber_dict["ForceField"]["HarmonicBondForce"]["Bond"]
all_lines.extend(extract_bonds(bonds))
all_lines.append("\n" "! Ligand vdW parameters\n" "! End ligand vdW parameters\n" "\n" "[angles]\n")

angles = phoosa_dict["ForceField"]["HarmonicAngleForce"]["Angle"]
all_lines.extend(extract_angles(angles))
angles = amber_dict["ForceField"]["HarmonicAngleForce"]["Angle"]
all_lines.extend(extract_angles(angles))
all_lines.append(
    "\n"
    "! Zero fc dual topology angles\n"
    "3C           CX           2C                  0.0      120.0\n"
    "\n"
    "! Ligand angle parameters\n"
    "! End ligand angle parameters\n"
    "\n"
    "[torsions]\n"
)

torsions = phoosa_dict["ForceField"]["PeriodicTorsionForce"]["Proper"]
all_lines.extend(extract_torsions(torsions))
torsions = amber_dict["ForceField"]["PeriodicTorsionForce"]["Proper"]
all_lines.extend(extract_torsions(torsions))

print("".join(all_lines))
