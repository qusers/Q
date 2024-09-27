"""Module to hold utility visualization functions for QligFEP."""

from pathlib import Path
from typing import Optional, Union

import py3Dmol
from openff.toolkit import Molecule
from rdkit import Chem


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
    """

    def mols_to_molblock(mol):
        if isinstance(mol, Molecule):
            mol = mol.to_rdkit()
        return Chem.MolToMolBlock(mol)

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
        view.addModel(mols_to_molblock(mol), "sdf")
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
