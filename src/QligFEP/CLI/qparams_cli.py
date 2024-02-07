"""Module containing the QligFEP parameter writing functionalities."""

from ..openff2Q import OpenFF2Q
import argparse
# TODO: implement the other parameter writing functionalities

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Write parameters for QligFEP")
    parser.add_argument('-i', '--input', required = True, help = "Input (sdf) file.")
    parser.add_argument(
        '-FF',
        '--forcefield',
        dest="forcefield",
        default='OpenFF',
        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST', 'OpenFF'],
        help = "Forcefield to be used. Defaults to OpenFF."
    )
    return parser.parse_args()

def main(args: argparse.Namespace) -> None:
    if args.forcefield == 'OpenFF':
        openff2q = OpenFF2Q(args.input)
        openff2q.process_ligands()
    else:
        raise NotImplementedError("Forcefield not supported through this CLI yet")
    
def main_exe():
    args = parse_arguments()
    main(args)
    
if __name__ == '__main__':
    main_exe()