"""Set up a series of amino-acid mutations, one prepared sphere per mutation.

A ligand series shares a binding site, so the ligand pipeline prepares one sphere
and reuses it for every edge. Mutations do not share anything: they sit wherever
the residues sit, and a sphere built around one of them describes the others
badly or not at all. Two things depend on the centre, not just one -- whether the
mutated residue is inside the simulated region, and which *other* charges get
neutralised as out-of-sphere. Both have to be re-derived per mutation.

So the unit of work here is a directory per mutation, each holding its own
``qprep_prot`` run, and the series is validated as a whole before any of them is
created. The centre comes from one of two places:

``centre on the mutated residue`` (default)
    The sphere is built around the residue's CB (HA3 for glycine), which is what
    the published protocol does. The mutation is at the centre by construction,
    so it can never be out of sphere or in the restrained shell.

``centre fixed by the caller`` (``center=``)
    One centre for the whole series, for scanning residues around a bound ligand
    that must stay in the same sphere. Nothing then guarantees a residue is
    properly simulated, so every mutation is checked against the sphere up front
    and the series is refused if any of them reaches beyond the boundary.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import amino_acids
from .IO import get_force_field_paths, parse_lib
from .logger import logger
from .pdb_utils import read_pdb_to_dataframe
from .qresfep import parse_mutation

#: Atom the sphere is centred on when centring on the mutated residue. Glycine
#: has no CB, so the published protocol falls back to HA3.
CENTRE_ATOMS = ("CB", "HA3")

#: Residues that carry no force field parameters of their own and are therefore
#: not expected in the library: qprep adds the sphere's own solvent.
SOLVENT_RESIDUES = frozenset({"HOH", "SOL", "WAT", "TIP3", "T3P"})

#: Residues PyMOL's mutagenesis wizard can build. Protonation variants (ASH,
#: GLH, HIP, LYN, ...) are Q's names, not PyMOL's, so they have to be supplied
#: as ready-made residue PDBs instead.
PYMOL_RESIDUES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
)

_MUTAGENESIS_BLOCK = """
reinitialize
cmd.load('{structure}')
cmd.wizard('mutagenesis')
cmd.do('refresh_wizard')
cmd.get_wizard().set_mode('{mutant}')
cmd.get_wizard().do_select('chain {chain} and resi {position}')
cmd.get_wizard().apply()
save {mutant}{position}.pdb, chain {chain} and resi {position}
"""


class SetupError(Exception):
    """A mutation series that cannot be set up as asked."""


@dataclass(frozen=True)
class Mutation:
    """One requested point mutation, in the input PDB's numbering."""

    wild_type: str
    position: int
    mutant: str
    chain: str

    @classmethod
    def from_string(cls, mutation: str, chain: str) -> Mutation:
        """Build from ``LEU39ALA`` or ``L39A``.

        Raises:
            SetupError: If the string does not read as a mutation.
        """
        try:
            wild_type, position, mutant = parse_mutation(mutation)
        except Exception as error:  # parse_mutation raises MutationError
            raise SetupError(str(error)) from error
        return cls(wild_type=wild_type, position=position, mutant=mutant, chain=chain)

    @property
    def name(self) -> str:
        """Directory-and-file name of the mutation, e.g. ``LEU39ALA``."""
        return f"{self.wild_type}{self.position}{self.mutant}"

    @property
    def mutant_pdb_name(self) -> str:
        """Name of the PDB holding the mutant residue alone, e.g. ``ALA39.pdb``."""
        return f"{self.mutant}{self.position}.pdb"


def read_mutations(path: Path, chain: str) -> list[Mutation]:
    """Read a mutation list: one mutation per line, blanks and ``#`` ignored.

    Raises:
        SetupError: If any entry does not read as a mutation.
    """
    entries = []
    for line in Path(path).read_text().splitlines():
        entry = line.split("#")[0].strip()
        if entry:
            entries.append(entry)
    if not entries:
        raise SetupError(f"{path} lists no mutations")
    return [Mutation.from_string(entry, chain) for entry in entries]


# ----------------------------------------------------------------------
# Reading the structure
# ----------------------------------------------------------------------


def residue_atoms(frame: pd.DataFrame, position: int, chain: str) -> pd.DataFrame:
    """Return the atoms of one residue of the input structure."""
    return frame[
        (frame["residue_seq_number"] == position) & (frame["chain_id"].str.strip() == chain)
    ]


def residue_center(frame: pd.DataFrame, position: int, chain: str) -> list[float]:
    """Return the sphere centre for a mutation of this residue.

    The centre is the residue's CB, or HA3 when it has none -- the atom the
    published protocol centres on.

    Raises:
        SetupError: If the residue is absent, or carries neither atom.
    """
    atoms = residue_atoms(frame, position, chain)
    if atoms.empty:
        raise SetupError(f"residue {position} of chain {chain} is not in the structure")
    for name in CENTRE_ATOMS:
        hit = atoms[atoms["atom_name"].str.strip() == name]
        if not hit.empty:
            return [float(hit.iloc[0][axis]) for axis in ("x", "y", "z")]
    raise SetupError(
        f"residue {position} of chain {chain} has no {' or '.join(CENTRE_ATOMS)} atom to "
        "centre the sphere on"
    )


def residue_reach(frame: pd.DataFrame, position: int, chain: str, center: list[float]) -> float:
    """Return how far the *farthest* atom of a residue sits from a centre.

    The farthest atom rather than the nearest: the side chain is what is being
    perturbed, so a residue whose tip reaches past the boundary is not one this
    sphere simulates properly, however close its backbone sits.
    """
    atoms = residue_atoms(frame, position, chain)
    if atoms.empty:
        return float("nan")
    coordinates = atoms[["x", "y", "z"]].to_numpy(dtype=float)
    return float(np.linalg.norm(coordinates - np.asarray(center, dtype=float), axis=1).max())


def library_residues(*library_paths: Path) -> set[str]:
    """Return every residue name defined in the given Q library files."""
    names: set[str] = set()
    for path in library_paths:
        names.update(re.findall(r"^\{(\w+)\}", Path(path).read_text(), flags=re.MULTILINE))
    return names


# ----------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------


def unknown_residues(frame: pd.DataFrame, known: set[str]) -> dict[str, list[str]]:
    """Return the structure's residues that no library defines.

    ``qprep`` reports these one atom at a time, deep inside a run that is already
    building topologies, so they are worth catching before any directory exists.

    Returns:
        Residue name to the ``chain:number`` labels it appears under.
    """
    missing: dict[str, list[str]] = {}
    solute = frame[~frame["residue_name"].str.strip().isin(SOLVENT_RESIDUES)]
    keys = ["residue_name", "chain_id", "residue_seq_number"]
    for name, chain, number in solute[keys].drop_duplicates().itertuples(index=False):
        name = str(name).strip()
        if name not in known:
            missing.setdefault(name, []).append(f"{str(chain).strip()}:{number}")
    return missing


def validate_series(
    frame: pd.DataFrame,
    mutations: list[Mutation],
    force_field: str,
    radius: float,
    center: list[float] | None = None,
    neutralization_offset: float = 3.0,
    cofactor_libraries: list[Path] | None = None,
    mutant_pdb_dir: Path | None = None,
) -> None:
    """Check a whole series before any of it is set up.

    Every problem is reported at once. A series is usually run unattended over
    twenty mutations, and finding out about the fourteenth after the thirteenth
    has finished is not worth the setup time.

    Args:
        frame: The input structure, already parsed.
        mutations: The requested mutations.
        force_field: Force field name or path.
        radius: Sphere radius in Angstrom.
        center: Fixed sphere centre, or None to centre on each mutated residue.
        neutralization_offset: ``qprep_prot``'s neutralisation offset; charged
            residues beyond ``radius - offset`` are prepared in a neutral form.
        cofactor_libraries: ``.lib`` files whose residues are also expected in
            the structure.
        mutant_pdb_dir: Directory of ready-made mutant residue PDBs, when the
            caller supplies them instead of having PyMOL build them.

    Raises:
        SetupError: If anything about the series would not set up correctly.
    """
    lib_path, _ = get_force_field_paths(force_field)
    known = set(parse_lib(force_field))
    if cofactor_libraries:
        known |= library_residues(*cofactor_libraries)

    problems: list[str] = []

    missing = unknown_residues(frame, known)
    if missing:
        listed = "; ".join(
            f"{name} ({', '.join(where[:4])}{', ...' if len(where) > 4 else ''})"
            for name, where in sorted(missing.items())
        )
        problems.append(
            f"{len(missing)} residue name(s) in the structure are not defined in "
            f"{Path(lib_path).name} or any cofactor library: {listed}. qprep cannot build a "
            "topology for these -- rename them to the force field's own names (HIS as "
            "HID/HIE/HIP, for instance), pass their libraries with -cof, or run pdb2amber."
        )

    boundary = radius - neutralization_offset
    seen: set[str] = set()
    for mutation in mutations:
        if mutation.name in seen:
            problems.append(f"{mutation.name} is listed more than once")
            continue
        seen.add(mutation.name)

        atoms = residue_atoms(frame, mutation.position, mutation.chain)
        if atoms.empty:
            problems.append(
                f"{mutation.name}: residue {mutation.position} of chain {mutation.chain} is "
                "not in the structure"
            )
            continue

        prepared_name = str(atoms.iloc[0]["residue_name"]).strip()
        if prepared_name != mutation.wild_type:
            problems.append(
                f"{mutation.name}: residue {mutation.position} of chain {mutation.chain} is "
                f"{prepared_name}, not {mutation.wild_type}"
            )

        if mutation.mutant not in amino_acids.SIDE_CHAINS:
            problems.append(
                f"{mutation.name}: {mutation.mutant} is not a residue QresFEP can build a "
                "hybrid for"
            )
        elif mutation.mutant not in known:
            problems.append(
                f"{mutation.name}: {mutation.mutant} is not defined in {Path(lib_path).name}"
            )

        if mutant_pdb_dir is not None:
            if not (Path(mutant_pdb_dir) / mutation.mutant_pdb_name).exists():
                problems.append(
                    f"{mutation.name}: {mutation.mutant_pdb_name} not found in {mutant_pdb_dir}"
                )
        elif mutation.mutant not in PYMOL_RESIDUES:
            problems.append(
                f"{mutation.name}: PyMOL's mutagenesis wizard cannot build {mutation.mutant}. "
                f"Build {mutation.mutant_pdb_name} yourself and pass its directory with "
                "--mutant-pdbs."
            )

        if center is None:
            # Centring on the residue: the only thing that can go wrong is having
            # no atom to centre on.
            try:
                residue_center(frame, mutation.position, mutation.chain)
            except SetupError as error:
                problems.append(f"{mutation.name}: {error}")
            continue

        reach = residue_reach(frame, mutation.position, mutation.chain, center)
        if reach > radius:
            problems.append(
                f"{mutation.name}: reaches {reach:.1f} A from the sphere centre, beyond the "
                f"{radius:.0f} A radius, so part of it is not in the simulated system"
            )
        elif reach > boundary:
            problems.append(
                f"{mutation.name}: reaches {reach:.1f} A from the sphere centre, into the "
                f"{neutralization_offset:.0f} A restrained shell below the {radius:.0f} A "
                "radius, where the boundary dominates the environment and charged residues "
                "are neutralised"
            )

    if problems:
        summary = "\n  ".join(problems)
        hint = ""
        if center is not None and any("A from the sphere centre" in p for p in problems):
            hint = (
                "\n\nThose residues need a sphere of their own. Drop the fixed centre to "
                "centre each sphere on its own mutated residue, or enlarge the radius."
            )
        raise SetupError(f"{len(problems)} problem(s) with this mutation series:\n  {summary}{hint}")


# ----------------------------------------------------------------------
# Mutant residues
# ----------------------------------------------------------------------


def write_mutagenesis_script(structure: Path, mutations: list[Mutation], path: Path) -> Path:
    """Write the PyMOL script that builds each mutant side chain as its own PDB.

    QresFEP needs the mutant residue positioned on the wild-type backbone, which
    is what PyMOL's mutagenesis wizard produces.
    """
    blocks = [
        _MUTAGENESIS_BLOCK.format(
            structure=Path(structure).name,
            mutant=mutation.mutant,
            position=mutation.position,
            chain=mutation.chain,
        )
        for mutation in dict.fromkeys(mutations)
    ]
    blocks.append("cmd.quit()\n")
    Path(path).write_text("".join(blocks))
    return Path(path)


def build_mutant_residues(structure: Path, mutations: list[Mutation], directory: Path) -> Path:
    """Build one PDB per mutant residue in `directory`, using PyMOL.

    Returns:
        The directory the residue PDBs were written to.

    Raises:
        SetupError: If PyMOL is missing, fails, or writes nothing for a mutation.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if shutil.which("pymol") is None:
        raise SetupError(
            "pymol is not on PATH. It builds the mutant side chains; install it, or build "
            "the residue PDBs yourself and pass their directory with --mutant-pdbs."
        )

    shutil.copy(structure, directory / Path(structure).name)
    script = write_mutagenesis_script(structure, mutations, directory / "mutagenesis.pml")
    logger.info(f"Building {len({m.mutant_pdb_name for m in mutations})} mutant residue(s) with PyMOL")
    result = subprocess.run(
        ["pymol", "-cq", script.name],
        cwd=directory,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise SetupError(f"pymol failed while building the mutant residues:\n{result.stderr.strip()}")

    missing = [m.name for m in mutations if not (directory / m.mutant_pdb_name).exists()]
    if missing:
        raise SetupError(
            f"PyMOL wrote no residue PDB for: {', '.join(missing)}.\n{result.stdout.strip()}"
        )
    return directory


# ----------------------------------------------------------------------
# Running the series
# ----------------------------------------------------------------------


@dataclass
class MutationOutcome:
    """What happened to one mutation of the series."""

    mutation: Mutation
    center: list[float]
    directory: Path
    legs: list[str]
    error: str | None = None

    @property
    def ok(self) -> bool:
        return self.error is None


class MutationSeries:
    """Set up every mutation of a series, each in its own prepared sphere."""

    def __init__(
        self,
        structure: Path,
        mutations: list[Mutation],
        force_field: str,
        cluster: str,
        radius: float = 25.0,
        center: list[float] | None = None,
        neutralization_offset: float = 3.0,
        legs: tuple[str, ...] = ("protein", "tripeptide"),
        cofactors: list[str] | None = None,
        mutant_pdb_dir: Path | None = None,
        workdir: Path | None = None,
        qresfep_options: list[str] | None = None,
    ):
        """
        Args:
            structure: The prepared protein PDB the mutations are numbered against.
            mutations: The requested mutations.
            force_field: Force field name or path, used for both preparation and setup.
            cluster: Cluster profile the run scripts are written for.
            radius: Sphere radius in Angstrom.
            center: Fixed sphere centre, or None to centre on each mutated residue.
            neutralization_offset: Passed to ``qprep_prot`` and used to decide
                which residues sit in the restrained shell.
            legs: Which legs of the thermodynamic cycle to set up.
            cofactors: Cofactor basenames; each needs ``.pdb``, ``.lib`` and
                ``.prm`` beside the structure.
            mutant_pdb_dir: Directory of ready-made mutant residue PDBs. When
                None, PyMOL builds them.
            workdir: Where the series is written. Defaults to the cwd.
            qresfep_options: Extra arguments passed through to ``qresfep``.
        """
        self.structure = Path(structure).resolve()
        self.mutations = list(mutations)
        self.force_field = force_field
        self.cluster = cluster
        self.radius = float(radius)
        self.center = [float(c) for c in center] if center else None
        self.neutralization_offset = float(neutralization_offset)
        self.legs = tuple(legs)
        self.cofactors = list(cofactors or [])
        self.mutant_pdb_dir = Path(mutant_pdb_dir).resolve() if mutant_pdb_dir else None
        self.workdir = Path(workdir or Path.cwd()).resolve()
        self.qresfep_options = list(qresfep_options or [])

        self.frame = read_pdb_to_dataframe(str(self.structure))

    @property
    def cofactor_files(self) -> list[Path]:
        """Every file the cofactors contribute, beside the structure."""
        return [
            self.structure.parent / f"{cofactor}{extension}"
            for cofactor in self.cofactors
            for extension in (".pdb", ".lib", ".prm")
        ]

    def validate(self) -> None:
        """Check the whole series, and the files it needs, before setting any of it up.

        Raises:
            SetupError: If anything would not set up correctly.
        """
        missing = [str(path) for path in self.cofactor_files if not path.exists()]
        if missing:
            raise SetupError(f"Missing cofactor file(s): {', '.join(missing)}")

        validate_series(
            self.frame,
            self.mutations,
            force_field=self.force_field,
            radius=self.radius,
            center=self.center,
            neutralization_offset=self.neutralization_offset,
            cofactor_libraries=[p for p in self.cofactor_files if p.suffix == ".lib"],
            mutant_pdb_dir=self.mutant_pdb_dir,
        )

    def center_for(self, mutation: Mutation) -> list[float]:
        """Return the sphere centre this mutation is prepared around."""
        if self.center is not None:
            return self.center
        return residue_center(self.frame, mutation.position, mutation.chain)

    def run(self) -> list[MutationOutcome]:
        """Validate, build the mutant residues, and set up every mutation.

        A mutation that fails does not stop the series -- the rest are still worth
        having -- but the failure is recorded and reported by the caller.
        """
        self.validate()

        residues = self.mutant_pdb_dir
        if residues is None:
            residues = build_mutant_residues(
                self.structure, self.mutations, self.workdir / "mutants"
            )

        for leg in self.legs:
            (self.workdir / leg).mkdir(parents=True, exist_ok=True)

        outcomes = []
        for mutation in self.mutations:
            outcomes.append(self._setup_one(mutation, residues))
        return outcomes

    def _setup_one(self, mutation: Mutation, residues: Path) -> MutationOutcome:
        """Prepare a sphere for one mutation and set up its legs in it."""
        directory = self.workdir / "work" / mutation.name
        center = self.center_for(mutation)
        outcome = MutationOutcome(mutation=mutation, center=center, directory=directory, legs=[])

        logger.info(
            f"=== {mutation.name} (chain {mutation.chain}) -- sphere centred on "
            f"[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]"
        )
        try:
            if directory.exists():
                logger.info(f"Replacing the existing {directory}")
                shutil.rmtree(directory)
            directory.mkdir(parents=True)
            shutil.copy(self.structure, directory / "protein.pdb")
            shutil.copy(residues / mutation.mutant_pdb_name, directory)
            for path in self.cofactor_files:
                shutil.copy(path, directory)

            self._run(self._qprep_command(center), directory)
            for leg in self.legs:
                self._run(self._qresfep_command(mutation, leg), directory)
                # Both legs are written as FEP_<mutation>, so each has to be moved
                # out before the next one is set up.
                destination = self.workdir / leg / f"FEP_{mutation.name}"
                if destination.exists():
                    logger.info(f"Replacing the existing {destination}")
                    shutil.rmtree(destination)
                shutil.move(str(directory / f"FEP_{mutation.name}"), str(destination))
                outcome.legs.append(leg)
        except (SetupError, OSError, subprocess.SubprocessError) as error:
            outcome.error = str(error)
            logger.error(f"{mutation.name} failed: {error}")
        return outcome

    def _qprep_command(self, center: list[float]) -> list[str]:
        command = [
            "qprep_prot",
            "-i", "protein.pdb",
            "-FF", self.force_field,
            "-r", f"{self.radius}",
            "-nbo", f"{self.neutralization_offset}",
            "-cog", *[f"{value:.3f}" for value in center],
        ]
        if self.cofactors:
            command += ["-cof", *[f"{cofactor}.pdb" for cofactor in self.cofactors]]
        return command

    def _qresfep_command(self, mutation: Mutation, leg: str) -> list[str]:
        command = [
            "qresfep",
            "-m", mutation.name,
            "-mc", mutation.chain,
            "-S", leg,
            "-FF", self.force_field,
            "-c", self.cluster,
        ]
        if self.cofactors:
            command += ["-cof", *self.cofactors]
        return command + self.qresfep_options

    @staticmethod
    def _run(command: list[str], directory: Path) -> None:
        """Run one setup step, failing loudly.

        Raises:
            SetupError: If the command is missing or exits non-zero.
        """
        logger.debug(f"[{directory.name}] {' '.join(command)}")
        try:
            result = subprocess.run(command, cwd=directory, capture_output=True, text=True)
        except FileNotFoundError as error:
            raise SetupError(
                f"{command[0]} is not on PATH -- reinstall the package (`pip install -e .`)"
            ) from error
        if result.returncode != 0:
            tail = (result.stderr or result.stdout).strip().splitlines()[-15:]
            raise SetupError(
                f"`{' '.join(command)}` failed in {directory}:\n" + "\n".join(tail)
            )
