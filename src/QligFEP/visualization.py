"""Module to hold utility visualization functions for QligFEP."""

from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import py3Dmol
from kartograf import SmallMoleculeComponent
from matplotlib.colors import rgb2hex
from openff.toolkit import Molecule
from rdkit import Chem

from .CLI.utils import get_avail_restraint_methods
from .logger import logger
from .restraints.restraint_setter import RestraintSetter


def mol_to_molblock(mol: Union[Molecule, Chem.Mol]) -> str:
    if isinstance(mol, Molecule):
        mol = mol.to_rdkit()
    return Chem.MolToMolBlock(mol)


def render_system(
    molecules: list[Union[Molecule, Chem.Mol]],
    protein_path: Optional[Union[Path, str]] = None,
    protein_style: str = "stick",
    size: tuple[int, int] = (600, 500),
) -> None:
    """Render a system containing molecules and a protein structure using py3Dmol.

    Args:
        molecules: list of Molecule or Chem.Mol objects to render.
        protein_path: path to a protein structure file to render alongside the molecules.
            Defaults to None.
        protein_style: style to render the protein structure in. Can be 'stick' or 'cartoon'.
            Defaults to 'stick'.

    Raises:
        ValueError: If protein_path is not a string or Path object.
        FileNotFoundError: If the protein file does not exist at the given path.
    """

    if protein_path is not None:
        if isinstance(protein_path, str):
            protein_path = Path(protein_path)
        if not isinstance(protein_path, Path):
            raise ValueError("protein_path should be a string or a Path object")
        elif not protein_path.exists():
            raise FileNotFoundError(f"Could not find the protein file at {protein_path}")

        view = py3Dmol.view(width=size[0], height=size[1])
        view.addModel(protein_path.read_text(), "pdb")
        if protein_style == "stick":
            view.setStyle({"model": -1}, {"stick": {"radius": 0.04}})
        elif protein_style == "cartoon":
            view.setStyle({"model": -1}, {"cartoon": {"color": "spectrum"}})
        view.addSurface(py3Dmol.VDW, {"opacity": 0.5}, {"model": -1})

    for i, mol in enumerate(molecules):
        view.addModel(mol_to_molblock(mol), "sdf")
        view.setStyle({"model": i + 1}, {"stick": {"colorscheme": "cyanCarbon"}})

    view.zoomTo({"model": -1})
    view.setHoverable(
        {},
        True,
        """function(atom,viewer,event,container) {
        if(!atom.label) {
        atom.label = viewer.addLabel(atom.resn+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
        }}""",
        """function(atom,viewer) {
        if(atom.label) {
        viewer.removeLabel(atom.label);
        delete atom.label;
        }
    }""",
    )
    view.render()
    view.show()


# Credit to: https://github.com/OpenFreeEnergy/openfe/blob/main/openfe/utils/visualization_3D.py
# And to: https://github.com/OpenFreeEnergy/kartograf/blob/main/src/kartograf/utils/mapping_visualization_widget.py
def render_ligand_restraints(
    ligand1: Union[Chem.Mol, Molecule, SmallMoleculeComponent],
    ligand2: Union[Chem.Mol, Molecule, SmallMoleculeComponent],
    restraint_mapping: dict[int, int],
    show_atom_idxs: bool = True,
    size: tuple[int, int] = (900, 500),
    sphere_palette: str = "hsv",
    sphere_radius: float = 0.6,
    sphere_alpha: float = 0.8,
    render: bool = True,
) -> None:
    """
    Visualize two ligands with their atom mappings in three linked windows. The first
    and the third windows show the ligands same-colored spheres rendered on the atoms that
    are restrained. The second window shows the overlapped ligands, with a lower opacity
    to show the space overlap.

    Args:
        ligand1: First ligand molecule (RDKit Mol or custom Molecule object).
        ligand2: Second ligand molecule (RDKit Mol or custom Molecule object).
        restraint_mapping: Dictionary of atom indexes; [0] - atom idxs in lig1 that are mapped
            into atom idxs in lig2 [1].
        show_atom_idxs: Whether to show atom idxs in the visualization.
        size: Tuple of (width, height) for the entire viewer.
        sphere_palette: Matplotlib colormap to use for coloring the spheres.
        sphere_radius: Radius of the spheres rendered in the restrained atoms. Defaults to 0.6.
        sphere_alpha: Alpha value for the spheres. Defaults to 0.8.
        render: Whether to render the visualization. Defaults to True.

    Returns:
        None
    """

    def add_spheres(view, mol, mapping, is_ligand1=True, viewer=(0, 0)):
        cmap = plt.get_cmap(sphere_palette, len(mapping))
        for i, (idx1, idx2) in enumerate(mapping.items()):
            idx = idx1 if is_ligand1 else idx2
            p = mol.GetConformer().GetAtomPosition(idx)
            color = rgb2hex(cmap(i))
            view.addSphere(
                {
                    "center": {"x": p.x, "y": p.y, "z": p.z},
                    "radius": sphere_radius,
                    "color": color,
                    "alpha": sphere_alpha,
                },
                viewer=viewer,
            )

    ligand1 = ligand1.to_rdkit() if hasattr(ligand1, "to_rdkit") else ligand1
    ligand2 = ligand2.to_rdkit() if hasattr(ligand2, "to_rdkit") else ligand2

    view = py3Dmol.view(width=size[0], height=size[1], viewergrid=(1, 3), linked=True)

    # Ligand 1 view
    view.addModel(mol_to_molblock(ligand1), "ligand1", viewer=(0, 0))
    view.setStyle({"stick": {}}, viewer=(0, 0))
    add_spheres(view, ligand1, restraint_mapping, is_ligand1=True, viewer=(0, 0))

    # Overlapped view
    view.addModel(mol_to_molblock(ligand1), "ligand1_overlap", viewer=(0, 1))
    view.addModel(mol_to_molblock(ligand2), "ligand2_overlap", viewer=(0, 1))
    view.setStyle({"model": 0}, {"stick": {"opacity": 0.7}}, viewer=(0, 1))
    view.setStyle(
        {"model": 1}, {"stick": {"opacity": 0.7, "colorscheme": "lightskyblueCarbon"}}, viewer=(0, 1)
    )  # available colorschemes: https://github.com/3dmol/3Dmol.js/blob/master/src/colors.ts

    # Ligand 2 view
    view.addModel(mol_to_molblock(ligand2), "ligand2", viewer=(0, 2))
    view.setStyle({"stick": {"colorscheme": "lightskyblueCarbon"}}, viewer=(0, 2))
    add_spheres(view, ligand2, restraint_mapping, is_ligand1=False, viewer=(0, 2))

    if show_atom_idxs:
        for v in [(0, 0), (0, 2)]:
            view.addPropertyLabels(
                "index",
                {},
                {
                    "fontColor": "black",
                    "font": "sans-serif",
                    "fontSize": "12",
                    "showBackground": False,
                    "alignment": "center",
                },
                viewer=v,
            )

    view.zoomTo()
    view.render()
    if render:
        view.show()
    return view


def apply_and_render_restraint(
    ligand1: Union[Molecule, Chem.Mol, str, Path],
    ligand2: Union[Molecule, Chem.Mol, str, Path],
    restraint_method: str = "hybridization_p",
    show_atom_idxs: bool = True,
    size: tuple[int, int] = (900, 500),
    sphere_palette: str = "hsv",
    sphere_radius: float = 0.6,
    sphere_alpha: float = 0.8,
) -> dict[int, int]:
    """Apply and render restraints between two ligands according to restraint methods
    available in QligFEP.

    Args:
        ligand1: First ligand (RDKit Mol or custom Molecule object).
        ligand2: Second ligand (RDKit Mol or custom Molecule object).
        restraint_method: Chosen method for setting the ligand restraints. Defaults to "hybridization_p".
        show_atom_idxs: Whether to show atom idxs in the visualization.
        size: Tuple of (width, height) for the entire viewer.
        sphere_palette: Matplotlib colormap to use for coloring the spheres.
        sphere_radius: Radius of the spheres rendered in the restrained atoms. Defaults to 0.6.
        sphere_alpha: Alpha value for the spheres. Defaults to 0.8.

    Returns:
        dict[int, int]: Dictionary of atom indexes containing the atoms to be restrained.
            [0] - atom indexes in ligand1;
            [1] - atom indexes in ligand2.

    Note: using the `show_atom_idxs` parameter will add the indexes to the rendering.
    """
    ligand1 = RestraintSetter.input_to_small_molecule_component(ligand1)
    ligand2 = RestraintSetter.input_to_small_molecule_component(ligand2)

    if restraint_method not in get_avail_restraint_methods():
        raise ValueError(
            f"restraint_method should be one of {get_avail_restraint_methods()}, got {restraint_method}"
        )
    elif restraint_method == "overlap":
        raise ValueError("Overlap method is not supported by this method yet, use `kartograf` instead.")

    rsetter = RestraintSetter(ligand1, ligand2)

    if restraint_method == "kartograf":
        restraint_dict = rsetter.set_restraints(
            kartograf_native=True,  # other arguments ignored
        )
    else:
        atom_compare_method, permissiveness_lvl = restraint_method.split("_")
        if permissiveness_lvl == "p":
            # in this case, ignore_surround_atom_type is not needed
            params = {"strict_surround": False}
        elif permissiveness_lvl == "ls":
            params = {"strict_surround": True, "ignore_surround_atom_type": True}
        elif permissiveness_lvl == "strict":
            params = {"strict_surround": True, "ignore_surround_atom_type": False}
        restraint_dict = rsetter.set_restraints(atom_compare_method=atom_compare_method, **params)
        logger.debug(f"Restraints set using {restraint_method} method. Parameters: {params}")

    render_ligand_restraints(
        ligand1,
        ligand2,
        restraint_dict,
        show_atom_idxs=show_atom_idxs,
        size=size,
        sphere_palette=sphere_palette,
        sphere_radius=sphere_radius,
        sphere_alpha=sphere_alpha,
    )

    return restraint_dict
