"""Replica exchange FEP for one leg and one replicate: command line and input directory.

The input directory holds the topology, the FEP file and the md template (see REQUIRED).
"""
import argparse
from dataclasses import dataclass
from pathlib import Path

# name in Inputs -> file in the input directory
REQUIRED = {
    "topology": "dualtop.top",
    "fep": "FEP1.fep",
    "template": "md_0500_0500.inp",
}


@dataclass(frozen=True)
class Inputs:
    """The files of one input directory."""
    directory: Path
    topology: Path
    fep: Path
    template: Path


def read_inputs(directory):
    """Check that the input directory is complete and return its files."""
    directory = Path(directory).resolve()
    if not directory.is_dir():
        raise SystemExit(f"{directory} is not a directory")
    missing = [name for name in REQUIRED.values() if not (directory / name).is_file()]
    if missing:
        raise SystemExit(f"{directory} lacks: {', '.join(missing)}")
    return Inputs(directory, **{key: directory / name for key, name in REQUIRED.items()})


def _scales(text):
    scales = tuple(float(x) for x in text.split(",") if x.strip())
    if any(not 0.0 < x < 1.0 for x in scales):
        raise argparse.ArgumentTypeError("REST2 scales must be between 0 and 1")
    return scales


def _atoms(text):
    """'1435,1496' or '1430:1440,1496', 1-based and inclusive -> 0-based indices."""
    atoms = []
    for item in text.split(","):
        first, _, last = item.partition(":")
        atoms.extend(range(int(first) - 1, int(last or first)))
    return tuple(atoms)


def q_atoms(fep_file):
    """0-based topology indices of the Q atoms: the [atoms] section of a Q FEP file."""
    atoms, inside = [], False
    for line in open(fep_file):
        line = line.split("!")[0].strip()
        if line.startswith("["):
            inside = line.lower() == "[atoms]"
        elif inside and line:
            atoms.append(int(line.split()[1]) - 1)
    return tuple(atoms)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input_dir", type=read_inputs,
                        help="folder with " + ", ".join(REQUIRED.values()))
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--precision", default="spfp", choices=("spfp", "double"))
    parser.add_argument("--windows", type=int, required=True,
                        help="number of lambda windows, 0 and 1 included")
    parser.add_argument("--replicate", type=int, default=1,
                        help="which replicate this run is, from 1; it selects the seed")
    parser.add_argument("--out", type=Path, default=Path("."),
                        help="folder for the results; the run writes into its subfolder rep_<replicate>")
    parser.add_argument("--settle-cycles", type=int, default=100,
                        help="swap cycles before production, thrown away")
    parser.add_argument("--cycles", type=int, default=1000, help="swap cycles of production")
    parser.add_argument("--block", type=int, default=500, help="MD steps between swap attempts")
    parser.add_argument("--rest", type=_scales, nargs="+", default=(),
                        help="REST2 scales of extra windows at both end states, e.g. 0.8,0.64; each end "
                             "gets a branch of these windows next to it (default: none)")
    parser.add_argument("--release-rotors", action="store_true",
                        help="keep a restraint between the two ligand copies only where both atoms are in "
                             "a ring or bonded to one (drops e.g. the methyl of a methoxy group)")
    parser.add_argument("--release-joints", action="store_true",
                        help="release restraints between the ligand copies on atoms past an sp3 joint (ether "
                             "O/S, 4-connected C) seen from the perturbed ring: methyls, CH2-linked rings")
    parser.add_argument("--restrain-core", action="store_true",
                        help="keep restraints between the ligand copies only in the ring system that carries "
                             "the changed element and the atoms bonded to it (replaces --release-rotors)")
    parser.add_argument("--refine", action="store_true",
                        help="after settling, place the windows again from the measured work and settle again")
    parser.add_argument("--save-atoms", type=_atoms, default=(),
                        help="extra atoms whose coordinates are saved, 1-based: 1435,1496 or 1430:1440; "
                             "the Q atoms are always saved")
    args = parser.parse_args(argv)
    args.rest = tuple(x for group in args.rest for x in group)    # "0.8,0.64" or "0.8 0.64"
    if args.windows < 2:
        parser.error("--windows must be at least 2")
    if args.replicate < 1:
        parser.error("--replicate counts from 1")
    return args

