"""Module to hold utility visualization functions for QligFEP."""

from pathlib import Path
from typing import Optional, Union

import py3Dmol
from openff.toolkit import Molecule
from rdkit import Chem


def render_system(
    molecules: list[Union[Molecule, Chem.Mol]],
    protein_path: Optional[Union[Path, str]] = None,
) -> None:

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

        view = py3Dmol.view(width=600, height=500)
        view.addModel(protein_path.read_text(), "pdb")
        view.setStyle({"model": -1}, {"stick": {"radius": 0.04}})
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
    return view.render()
