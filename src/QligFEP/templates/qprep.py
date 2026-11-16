"""qprep.inp templates for topology generation.

Provides templates for qprep input files used in FEP setup (dual topology)
and protein preparation workflows.
"""

from dataclasses import dataclass


@dataclass
class QprepFEPParameters:
    """Parameters for FEP dual-topology qprep.inp."""

    ff_lib: str  # Force field library (e.g., "qamber14.lib")
    lig1_lib: str  # Ligand 1 library file
    lig2_lib: str  # Ligand 2 library file (renumbered)
    ligand_prm: str  # Merged parameter file
    ligand_pdb: str  # Combined ligand PDB
    center: str  # Center of geometry "x y z"
    sphere_radius: str  # Sphere radius
    solute_density: str  # Solute density value
    solvent: str  # Solvent specification (e.g., "1 HOH" or "4 water.pdb")
    cysbonds: str = ""  # Cysbond lines (empty string if none)
    solvate: bool = True  # Whether to solvate (False for vacuum)


@dataclass
class QprepResFEPParameters:
    """Parameters for a residue-mutation (QresFEP) qprep.inp.

    Differs from the ligand case in taking an open-ended list of libraries: the
    hybrid residue's own library has to be read before the force field's, and any
    cofactors add more.
    """

    libraries: list[str]  # Library files, read in order (hybrid residue first)
    prm: str  # Merged parameter file
    pdb: str  # Complex PDB (protein + hybrid residue, or the reference tripeptide)
    center: str  # Sphere center "x y z"
    sphere_radius: str  # Sphere radius
    solute_density: str  # Solute density value
    solvent: str  # Solvent specification ("4 water.pdb" or "1 HOH")
    solvent_pack: str = "2.5"
    cysbonds: str = ""  # addbond lines (empty string if none)
    solvate: bool = True  # False for vacuum


@dataclass
class QprepProteinParameters:
    """Parameters for protein prep qprep.inp."""

    ff_lib_path: str  # Force field library path
    ff_prm_path: str  # Force field parameter path
    pdb_file_path: str  # Input PDB file path
    cog: str  # Center of geometry "x y z"
    sphere_radius: str  # Sphere radius
    solvent_pack: str  # Solvent packing value
    cysbonds: str = ""  # Cysbond lines (empty string if none)


def render_qprep_fep_input(params: QprepFEPParameters) -> str:
    """Render qprep.inp for FEP dual-topology setup.

    Args:
        params: FEP qprep parameters

    Returns:
        Complete qprep.inp file content as string
    """
    lines = [
        f"rl {params.ff_lib}",
        f"rl {params.lig1_lib}",
        f"rl {params.lig2_lib}",
        f"rprm {params.ligand_prm}",
        "! TO DO Change if protein system is used",
        f"rp {params.ligand_pdb}",
        "! solvation - density adapted to membrane proteins",
        "set solvent_pack 2.3",
        f"set solute_density {params.solute_density}",
        f"boundary 1 {params.center} {params.sphere_radius}",
    ]

    # Solvate line (commented out for vacuum)
    if params.solvate:
        lines.append(f"solvate {params.center} {params.sphere_radius} {params.solvent}")
    else:
        lines.append(f"!solvate {params.center} {params.sphere_radius} {params.solvent}")

    # Cysbonds (if any) appear before maketop
    if params.cysbonds:
        lines.append(params.cysbonds.rstrip())

    lines.extend(
        [
            "maketop MKC_p",
            "writetop dualtop.top",
            "wp top_p.pdb y",
            "rt dualtop.top",
            "mask none",
            "mask not excluded",
            "wp complexnotexcluded.pdb y",
            "q",
        ]
    )

    return "\n".join(lines)


def render_qprep_resfep_input(params: QprepResFEPParameters) -> str:
    """Render qprep.inp for a QresFEP hybrid-residue topology.

    Args:
        params: ResFEP qprep parameters

    Returns:
        Complete qprep.inp file content as string
    """
    lines = [f"rl {library}" for library in params.libraries]
    lines += [
        f"rprm {params.prm}",
        f"rp {params.pdb}",
        f"set solvent_pack {params.solvent_pack}",
        f"set solute_density {params.solute_density}",
        f"boundary 1 {params.center} {params.sphere_radius}",
    ]

    solvate_line = f"solvate {params.center} {params.sphere_radius} {params.solvent}"
    lines.append(solvate_line if params.solvate else f"!{solvate_line}")

    if params.cysbonds:
        lines.append(params.cysbonds.rstrip())

    lines.extend(
        [
            "maketop MKC_p",
            "writetop dualtop.top",
            "wp top_p.pdb y",
            "rt dualtop.top",
            "mask none",
            "mask not excluded",
            "wp complexnotexcluded.pdb y",
            "q",
        ]
    )

    return "\n".join(lines) + "\n"


def render_qprep_protein_input(params: QprepProteinParameters) -> str:
    """Render qprep.inp for protein preparation.

    Args:
        params: Protein prep qprep parameters

    Returns:
        Complete qprep.inp file content as string
    """
    lines = [
        f"rl {params.ff_lib_path}",
        f"rprm {params.ff_prm_path}",
        "! TO DO Change if protein system is used",
        f"rp {params.pdb_file_path}",
        "! set solute_density 0.05794",
        "! NOTE, this is now large for water system, change for protein system",
        f"set solvent_pack {params.solvent_pack}",
        f"boundary 1 {params.cog} {params.sphere_radius}",
        f"solvate {params.cog} {params.sphere_radius} 1 HOH",
    ]

    # Cysbonds (if any) appear before maketop
    if params.cysbonds:
        lines.append(params.cysbonds.rstrip())

    lines.extend(
        [
            "maketop MKC_p",
            "writetop dualtop.top",
            "wp top_p.pdb y",
            "rt dualtop.top",
            "mask none",
            "mask not excluded",
            "wp complexnotexcluded.pdb y",
            "q",
        ]
    )

    return "\n".join(lines)
