from pathlib import Path

# openFF modules
from openff.toolkit import Molecule
from openff.toolkit.utils import UndefinedStereochemistryError

# QligFEP modules
from .logger import logger

class MoleculeIO(object):
    """A class to handle the input/output of molecules. Ligands are usually initialized
    from .sdf files and then processed into individual molecules.
    
    Attributes:
        self.lig: the input ligand file.
        self.molecules: a list of Molecule objects.
        self.lig_names: a list of ligand names.
        self.sdf_contents: a dictionary of the sdf contents for each ligand.
        self.single_ligand: a boolean to indicate if the input file contains a single ligand.
    """    
    def __init__(self, lig, *args, **kwargs):
        """Initialize a Molecule Input/Output object. This helper class has a base functionality
        used for handling `.sdf`, like reading it, outputting separate `.sdf` files (required by lomap),
        and storing the molecules & their names into a single object. 
        """        
        self.lig = lig
        self.molecules = []
        self.lig_names = []
        self.setup_mols_and_names()
        self.sdf_contents = {} # store the sdf content for each entry
        self.parse_sdf_contents() # add the sdf content to the dictionary
        
    def setup_mols_and_names(self):
        try:
            mols = Molecule.from_file(self.lig)
        except UndefinedStereochemistryError:
            logger.warning('Undefined stereochemistry in the input file!! Will try to process the ligands anyway.')
            mols = Molecule.from_file(self.lig, allow_undefined_stereo=True)
        if not isinstance(mols, list):
            self.single_ligand = True
            mols = [mols]
        else:
            self.single_ligand = False
        
        if self.single_ligand:
            self.lig_names = [self.lig.split('.')[0]]
        else:
            self.lig_names = [ # in some cases the name is empty so we give it a name
                mol.name if mol.name != '' else f'lig_{idx}' for idx, mol in enumerate(mols)
            ]
        self.molecules = mols

    def parse_sdf_contents(self):
        """Parse the SDF content into individual entries and have them saved in a dictionary."""
        with open(self.lig) as infile:
            content = infile.readlines()
            ligands = []
            current_ligand = []

            for line in content:
                if line.strip() == '$$$$':  # End of a ligand entry
                    ligands.append(current_ligand)
                    current_ligand = []  # Start a new ligand entry
                else:
                    current_ligand.append(line.strip())

        for lname, sdf_content in zip(self.lig_names, ligands):
            self.sdf_contents.update({lname : sdf_content})
            
    def write_sdf_separate(self, output_dir):
        """Function to write the separate multiple molecules within a sdf file into their own
        .sdf, placed under `output_dir`.

        Args:
            output_dir: the directory to save the single sdf molecules.

        Raises:
            ValueError: if `output_dir` is neither a `str` or a `pathlib.Path` object.
        """        
        if isinstance(output_dir, str):
            output_dir = Path(output_dir)
        elif not isinstance(output_dir, Path):
            raise ValueError('output_dir must be either a str or pathlib.Path object')
        if not output_dir.exists():
            output_dir.mkdir(exist_ok=True)
        for idx, lname in enumerate(self.lig_names):
            fpath = str(output_dir / f'{lname}.sdf')
            self.molecules[idx].to_file(file_path=fpath, file_format='sdf')
        