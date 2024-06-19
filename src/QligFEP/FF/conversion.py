"""Module containing helper functions to convert between different forcefield"""

import math


def array_calc(atom, sig, eps, mass=None):
    # calculate LJ types
    sig = float(sig)
    eps = float(eps)
    if mass is None:
        mass = 0.0
    return f"{atom:<13}{sig:>10}{0.0:>11}" f"{eps:>11}{sig:>11}{eps/2:>11.6f}{mass:>11.2f}"


def atom_type_forces_from_openmm(sig, eps):
    # convert eps from kJ/mol to kcal/mol
    sig = round(sig * 5.612, 4)
    # convert sig to A
    eps = round(eps * 23.9005736 / 100, 4)
    return sig, eps


def convert_kj_mol_nm2_to_kcal_mol_a2(value_kj_mol_nm2):
    conversion_factor = 23.9005736
    value_kcal_mol_a2 = round(value_kj_mol_nm2 * conversion_factor / 10000, 1)
    return value_kcal_mol_a2


def radians_to_degrees(radians):
    return radians * 180 / math.pi


# non-bonded forces from Amber14sb. Source:
# https://github.com/openmm/openmm/blob/8.1.1/wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml
diction = {
    "proteinC": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCA": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCB": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCC": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCN": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCR": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCT": {"eps": 0.4577296, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCV": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCW": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinStar": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinCX": {"eps": 0.4577296, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinH": {"eps": 0.06568879999999999, "sig": 0.10690784617684071, "mass": 1.008},
    "proteinHC": {"eps": 0.06568879999999999, "sig": 0.2649532787749369, "mass": 1.008},
    "proteinH1": {"eps": 0.06568879999999999, "sig": 0.2471353044121301, "mass": 1.008},
    "proteinHA": {"eps": 0.06276, "sig": 0.25996424595335105, "mass": 1.008},
    "proteinH4": {"eps": 0.06276, "sig": 0.2510552587719476, "mass": 1.008},
    "proteinH5": {"eps": 0.06276, "sig": 0.24214627159054422, "mass": 1.008},
    "proteinHO": {"eps": 0.0, "sig": 1.0, "mass": 1.008},
    "proteinHS": {"eps": 0.06568879999999999, "sig": 0.10690784617684071, "mass": 1.008},
    "proteinHP": {"eps": 0.06568879999999999, "sig": 0.19599771799087468, "mass": 1.008},
    "proteinN": {"eps": 0.7112800000000001, "sig": 0.3249998523775958, "mass": 14.01},
    "proteinNA": {"eps": 0.7112800000000001, "sig": 0.3249998523775958, "mass": 14.01},
    "proteinNB": {"eps": 0.7112800000000001, "sig": 0.3249998523775958, "mass": 14.01},
    "proteinN2": {"eps": 0.7112800000000001, "sig": 0.3249998523775958, "mass": 14.01},
    "proteinN3": {"eps": 0.7112800000000001, "sig": 0.3249998523775958, "mass": 14.01},
    "proteinO": {"eps": 0.87864, "sig": 0.2959921901149463, "mass": 16.0},
    "proteinO2": {"eps": 0.87864, "sig": 0.2959921901149463, "mass": 16.0},
    "proteinOH": {"eps": 0.8803136, "sig": 0.3066473387839048, "mass": 16.0},
    "proteinS": {"eps": 1.046, "sig": 0.35635948725613575, "mass": 32.06},
    "proteinSH": {"eps": 1.046, "sig": 0.35635948725613575, "mass": 32.06},
    "proteinCO": {"eps": 0.359824, "sig": 0.3399669508423535, "mass": 12.01},
    "protein2C": {"eps": 0.4577296, "sig": 0.3399669508423535, "mass": 12.01},
    "protein3C": {"eps": 0.4577296, "sig": 0.3399669508423535, "mass": 12.01},
    "proteinC8": {"eps": 0.4577296, "sig": 0.3399669508423535, "mass": 12.01},
}

if __name__ == "__main__":
    print("Reproducing the forcefield conversion from OpenMM to Q implementation")
    print("[atom_types]")
    print("* tac               R*1        R*2   epsilon1        R*3 epsilon2&3       mass")
    for key in diction:
        sig, eps = atom_type_forces_from_openmm(diction[key]["sig"], diction[key]["eps"])
        string = array_calc(key.replace("protein", ""), sig, eps, diction[key]["mass"])
        print(string)
    # only one representative to show we reproduce the conversion...
    print("[bonds]")
    k = 259407.99999999994
    length = 0.1525
    type1 = "C"
    type2 = "C"  # noqa: E702
    new_k = convert_kj_mol_nm2_to_kcal_mol_a2(k)
    new_len = length * 10
    format_line = f"{type1:<13}{type2:<13}{new_k:>10}{new_len:>11}"
    print(format_line)
    print("[angles]")
    angle = 2.0943951023931953
    new_angle = round(radians_to_degrees(angle), 1)
    k = 669.44
    new_k = convert_kj_mol_nm2_to_kcal_mol_a2(k) * 100
    type1 = "C"
    type2 = "C"
    type3 = "O"
    format_line = f"{type1:<13}{type2:<13}{type3:<13}{new_k:>10}{new_angle:>11}"
    print(format_line)
    print("[torsions]")
    data = {  # data from PeriodicTorsionForce -> Proper in openMM
        "k1": 0.0,
        "k2": 1.6736000000000002,
        "k3": 0.8368000000000001,
        "k4": 0.8368000000000001,
        "periodicity1": 4,
        "periodicity2": 3,
        "periodicity3": 2,
        "periodicity4": 1,
        "phase1": 0.0,
        "phase2": 0.0,
        "phase3": 0.0,
        "phase4": 0.0,
        "type1": "protein-CT",
        "type2": "protein-CT",
        "type3": "protein-C",
        "type4": "protein-N",
    }
    k_values = [data["k1"], data["k2"], data["k3"], data["k4"]]
    periodicities = [data["periodicity1"], data["periodicity2"], data["periodicity3"], data["periodicity4"]]
    phases = [data["phase1"], data["phase2"], data["phase3"], data["phase4"]]
    output = []
    for i in range(4):
        k = k_values[i] / 4.184  # Force constant in kcal/mol
        periodicity = periodicities[i]
        phase = radians_to_degrees(phases[i])  # Phase in degrees
        num_paths = 1  # Given as 1 in the example

        # Check if more components follow
        if i < 3:
            periodicity = -periodicity

        # Append to the output
        types = [data[f"type{j}"].replace("protein-", "") for j in range(1, 5)]
        output.append(
            f"{types[0]:<13}{types[1]:<13}{types[2]:<13}{types[3]:<13}{k:>10.2f}{periodicity:>6.1f}{phase:>11.1f}{num_paths:>6.1f}"
        )
    print("\n".join(output))
