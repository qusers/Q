"""Set up free energy perturbations for amino-acid mutations (QresFEP).

Follows the dual-topology protocol of Koenekoop et al. (2025): the mutated
position carries both side chains at once inside a single hybrid residue, and
the perturbation runs as two consecutive FEP stages.

  stage 1  the wild-type side chain is discharged while the mutant one is a
           chargeless soft-core ghost
  stage 2  the mutant side chain is charged and grown in, and the wild-type one
           is turned into a ghost

Stage 2 starts from stage 1's final coordinates, so the pair is one calculation
split in two. The shift in folding free energy comes from running the same
mutation twice -- once in the protein, once in a small capped reference peptide
-- and subtracting.

Inputs come from a finished ``qprep_prot`` run in the working directory
(the prepared protein PDB, ``water.pdb`` and ``prep.json``) plus a PDB of the
mutant residue alone, positioned on the wild-type backbone.
"""

from __future__ import annotations

import glob
import itertools
import json
import math
import os
import re
import shutil
import stat
import subprocess
from pathlib import Path

from . import amino_acids
from .functions import euclidian_overlap, resfep_lambda_ladder, sphere_solute_density
from .IO import AA, get_force_field_paths, read_prm, replace, run_qprep
from .logger import logger
from .pdb_utils import pdb_parse_in, pdb_parse_out
from .peptide_caps import PeptideBuildError, build_reference_peptide
from .resfep_protocols import DEFAULT_PRODUCTION_STEPS
from .settings.settings import CLUSTER_DICT, CONFIGS
from .sphere_prep import SpherePrep, residue_distance_from_center
from .templates import (
    QprepResFEPParameters,
    format_sequence_restraint,
    format_wall_restraints,
    get_production_config,
    get_resfep_equilibration_configs,
    render_md_input,
    render_qfep_input,
    render_qprep_resfep_input,
)

#: The two stages of a dual-topology mutation, in the order they must run.
FEP_STAGES = ("FEP1.fep", "FEP2.fep")

#: Residues carrying a formal charge. A mutation touching one changes the net
#: charge of the sphere, so the reference leg has to be charge-matched.
CHARGED_RESIDUES = frozenset({"ASP", "GLU", "LYS", "ARG", "HIP"})

#: Soft-core maximum potential for side-chain atoms in their ghost state, in
#: kcal/mol -- Q reads it as a potential because the FEP file sets
#: `softcore_use_max_potential on`.
SOFTCORE_MAX_POTENTIAL = "20"

#: Distance-restraint window and force holding the two side chains superimposed.
#: Pairs closer than `low` are free; the force acts beyond it.
RESTRAINT_LOW = 0.0
RESTRAINT_HIGH = 0.5
RESTRAINT_FORCE = 10.0

#: Alpha-carbon charges the OPLS libraries give glycine and every other residue.
#: A mutation to or from glycine moves CA between the two, in whichever stage
#: also handles that residue's HA atoms.
_CA_CHARGE_GLY = "0.080"
_CA_CHARGE_NON_GLY = "0.140"
_HA_CHARGE = "0.060"
_HA_TYPE = "H140"

#: Residues qprep adds as solvent, after the solute it was given.
_SOLVENT_RESIDUES = frozenset({"HOH", "SOL", "WAT"})

#: Parameter-file sections copied into both the temporary and final merged files.
_PARAMETER_SECTIONS = ("[options]", "[atom_types]", "[bonds]", "[angles]", "[torsions]", "[impropers]")

#: Local profiles cannot run the two chained FEP stages.
_LOCAL_CLUSTERS = frozenset({"LOCAL", "LOCALP"})


class MutationError(Exception):
    """Raised when a requested mutation cannot be set up."""


def parse_mutation(mutation: str) -> tuple[str, int, str]:
    """Split a mutation string into wild-type, position and mutant.

    Accepts one- or three-letter residue codes on either side, so ``L39A`` and
    ``LEU39ALA`` name the same mutation.

    Args:
        mutation: The mutation, e.g. ``LEU39ALA``.

    Returns:
        ``(wild_type, position, mutant)``, residues as three-letter codes.

    Raises:
        MutationError: If the string does not name a position, or either residue
            is not a known amino acid.
    """
    parts = re.split(r"(\d+)", mutation.strip())
    if len(parts) != 3 or not parts[1].isdigit():
        raise MutationError(
            f"Could not read {mutation!r}; expected <wild-type><position><mutant>, e.g. LEU39ALA"
        )

    residues = []
    for code in (parts[0], parts[2]):
        code = code.upper()
        if len(code) == 1:
            try:
                residues.append(AA(code))
            except KeyError as exc:
                raise MutationError(f"{code!r} in {mutation!r} is not an amino acid") from exc
        elif code in amino_acids.SIDE_CHAINS:
            residues.append(code)
        else:
            raise MutationError(f"{code!r} in {mutation!r} is not an amino acid")

    return residues[0], int(parts[1]), residues[1]


class QresFEP:
    """Generate the input files for one amino-acid mutation FEP."""

    def __init__(
        self,
        mutation: str,
        chain: str,
        system: str,
        force_field: str,
        cluster: str,
        windows: int = 25,
        sampling: str = "exponential",
        start: str = "1",
        timestep: str = "2fs",
        temperature: str = "298",
        replicates: int = 10,
        tripeptide_flanks: str = "A",
        eq5_steps: int | None = None,
        production_steps: int = DEFAULT_PRODUCTION_STEPS,
        cofactors: list[str] | None = None,
        seeds: list[int] | None = None,
        to_clean: list[str] | None = None,
        write_trajectories: bool = True,
        separate_scaling: bool = True,
        workdir: Path | None = None,
    ):
        """
        Args:
            mutation: Mutation in the input PDB's numbering, e.g. ``LEU39ALA``.
            chain: Chain of the mutated residue in the input PDB.
            system: ``protein`` for the bound leg, ``tripeptide`` for the reference.
            force_field: Force field name or path; must match the preparation.
            cluster: Cluster profile from ``CLUSTER_DICT``.
            windows: Lambda windows *per stage*; the total is twice this.
            sampling: Lambda spacing. ``exponential`` gives stage 1 exponential and
                stage 2 reverse-exponential spacing, concentrating windows where
                each stage is steep; any other choice uses one spacing for both.
            start: Lambda endpoint the simulations start from.
            timestep: ``1fs`` or ``2fs``.
            temperature: Temperature in K, or a comma-separated list.
            replicates: Independent repeats.
            tripeptide_flanks: Residues flanking the mutated one in the reference
                peptide: ``A`` (Ala), ``G`` (Gly), ``X`` (none) or ``Z`` (the
                native neighbours).
            eq5_steps: Length of the final equilibration in steps.
            production_steps: Production length per lambda window in steps.
            cofactors: Cofactor basenames; each needs a ``.pdb``, ``.lib`` and ``.prm``.
            seeds: Random seeds, one per replicate.
            to_clean: File suffixes to delete after each run.
            write_trajectories: Write DCD trajectories. Disable for quota-safe
                production when trajectories are not needed for analysis.
            separate_scaling: Scale solute and solvent independently.
            workdir: Directory holding the preparation output. Defaults to cwd.
        """
        self.workdir = Path(workdir or Path.cwd())
        self.wild_type, self.position, self.mutant = parse_mutation(mutation)
        self.chain = chain
        self.system = system
        self.force_field = force_field
        self.cluster = cluster
        self.windows = int(windows)
        self.sampling = sampling
        self.start = start
        self.timestep = timestep
        self.temperature = str(temperature)
        self.replicates = int(replicates)
        self.tripeptide_flanks = tripeptide_flanks
        self.eq5_steps = int(eq5_steps) if eq5_steps else None
        self.production_steps = int(production_steps)
        if self.production_steps <= 0:
            raise MutationError("Production steps must be a positive integer")
        self.cofactors = list(cofactors or [])
        self.to_clean = to_clean
        self.write_trajectories = bool(write_trajectories)
        self.separate_scaling = bool(separate_scaling)
        self.seeds = list(seeds) if seeds is not None else list(range(1, self.replicates + 1))
        if len(self.seeds) != self.replicates:
            raise MutationError(
                f"Expected {self.replicates} random seeds, received {len(self.seeds)}"
            )
        if any(not isinstance(seed, int) or not 1 <= seed <= 9999 for seed in self.seeds):
            raise MutationError("Random seeds must be integers from 1 through 9999")

        self.from_gly = self.wild_type == "GLY"
        self.to_gly = self.mutant == "GLY"
        self.backbone = amino_acids.backbone_atoms(self.from_gly, self.to_gly)
        self.hybrid_name = f"{AA(self.wild_type)}2{AA(self.mutant)}"

        self.prep: SpherePrep | None = None
        self.pdb: dict[int, list] = {}
        self.atom_ids: dict[str, str] = {}
        self.merged_lib: dict[str, list[str]] = {}
        self.zero_force: dict[str, list[str]] = {}
        self.anchor_atom: int | None = None
        self.counter_ions = 0
        self.system_size = 0

        self.name = f"{self.wild_type}{self.position}{self.mutant}"
        self.directory = self.workdir / f"FEP_{self.name}"
        self.inputfiles = self.directory / "inputfiles"

    # ------------------------------------------------------------------
    # Inputs
    # ------------------------------------------------------------------

    @property
    def mutant_pdb_path(self) -> Path:
        """PDB holding the mutant residue alone, on the wild-type backbone."""
        return self.workdir / f"{self.mutant}{self.position}.pdb"

    def read_prep(self) -> None:
        """Load the sphere preparation and resolve the mutation into Q numbering.

        Raises:
            FileNotFoundError: If a required input is missing.
            MutationError: If the residue at the requested position is not the
                stated wild type.
        """
        self.prep = SpherePrep.read(self.workdir)
        self.protein_pdb = self.workdir / self.prep.prepared_pdb
        self.radius = self.prep.radius
        self.center = self.prep.center
        self.shell_radius = self.radius

        # Check the mutation describes the prepared structure before checking for
        # files: naming the wrong wild type is the more fundamental mistake, and
        # reporting a missing PDB first would send the user off to build one that
        # was never going to work.
        self.q_position = self.prep.q_number(self.position, self.chain)
        self._check_residue_is_in_the_sphere()

        prepared_name = self.prep.residue_name(self.position, self.chain)
        if prepared_name != self.wild_type:
            raise MutationError(self._wrong_wild_type_message(prepared_name))

        required = [self.protein_pdb, self.workdir / "water.pdb", self.mutant_pdb_path]
        for cofactor in self.cofactors:
            required += [self.workdir / f"{cofactor}{ext}" for ext in (".pdb", ".lib", ".prm")]
        missing = [path.name for path in required if not path.exists()]
        if missing:
            raise FileNotFoundError(f"Missing required input file(s): {', '.join(missing)}")

        logger.info(
            f"{self.name} (chain {self.chain}) is Q residue {self.q_position}; "
            f"hybrid residue {self.hybrid_name}"
        )

    def _check_residue_is_in_the_sphere(self) -> None:
        """Refuse to mutate a residue the prepared sphere does not properly simulate.

        The sphere is built around one centre, and everything beyond its radius is
        outside the simulated system. Mutating a residue out there produces files
        that look ordinary and describe nothing: the perturbed side chain is not in
        the region being sampled.

        This has to be checked by geometry rather than by inference from the
        residue name. A charged residue outside the boundary is at least visible,
        because preparation renames it to a neutral form -- but a neutral residue
        is renamed by nothing, and would otherwise pass unnoticed.

        Mutating a set of residues therefore means preparing a sphere per residue,
        which also re-derives which *other* charges get neutralised: a charge left
        standing outside the sphere would poison the region that is simulated.

        Raises:
            MutationError: If the residue lies beyond the sphere radius.
        """
        distance = residue_distance_from_center(self.protein_pdb, self.q_position, self.center)
        if math.isnan(distance):
            logger.debug(f"Could not locate Q residue {self.q_position} to check its distance")
            return

        if distance > self.radius:
            raise MutationError(
                f"Residue {self.position} of chain {self.chain} sits {distance:.1f} A from the "
                f"sphere centre {self.center}, beyond the {self.radius:.0f} A radius, so it is "
                "not part of the simulated system. Prepare a sphere centred on this residue "
                "before mutating it.\n"
                "Mutating several residues means one `qprep_prot` run per residue: the sphere "
                "has to enclose the residue being mutated, and each centre also decides which "
                "other charges are neutralised."
            )

        if distance > self.prep.boundary_radius:
            logger.warning(
                f"Residue {self.position} of chain {self.chain} sits {distance:.1f} A from the "
                f"sphere centre, inside the {self.radius:.0f} A radius but within the "
                f"{self.prep.neutralization_offset:.0f} A restrained boundary shell. Its "
                "environment is dominated by the boundary rather than by the protein, and "
                "nearby charges there were neutralised during preparation. Re-centre the "
                "sphere on this residue unless you specifically want this."
            )

    def _wrong_wild_type_message(self, prepared_name: str) -> str:
        """Explain that the prepared residue is not the one the mutation names.

        Reached only for residues inside the sphere -- anything beyond it is
        rejected earlier, by :meth:`_check_residue_is_in_the_sphere`. So the cause
        here is always neutralisation: the residue lies in the restrained boundary
        shell, where preparation strips formal charges.
        """
        distance = residue_distance_from_center(self.protein_pdb, self.q_position, self.center)
        distance_note = (
            f" It sits {distance:.1f} A from the sphere centre, within the "
            f"{self.prep.neutralization_offset:.0f} A restrained boundary shell, where "
            "`qprep_prot` neutralises charged residues."
            if not math.isnan(distance)
            else ""
        )
        return (
            f"Residue {self.position} of chain {self.chain} is {prepared_name} in "
            f"{self.prep.prepared_pdb}, not {self.wild_type}.{distance_note} Re-centre the "
            f"sphere on this residue so it keeps its charge. Asking for {prepared_name}"
            f"{self.position}{self.mutant} instead would run, but it is a different "
            f"perturbation from {self.wild_type}{self.position}{self.mutant} -- only do that "
            "if the neutral form is what you mean to mutate."
        )

    def create_environment(self) -> None:
        """Create the FEP directory and copy in the libraries."""
        # Both legs of a mutation take the same directory name, so setting one up
        # on top of the other silently destroys the first. The legs are moved apart
        # (into protein/ and tripeptide/) between runs; warn if that has not happened.
        previous = self.inputfiles / "resfep_config.json"
        if previous.exists():
            try:
                earlier_system = json.loads(previous.read_text()).get("system")
            except json.JSONDecodeError:
                earlier_system = None
            if earlier_system and earlier_system != self.system:
                logger.warning(
                    f"{self.directory} already holds the {earlier_system} leg and is about to be "
                    f"overwritten with the {self.system} leg. Move the finished leg out of the way "
                    "first -- the analysis expects protein/ and tripeptide/ to hold one leg each."
                )

        self.inputfiles.mkdir(parents=True, exist_ok=True)
        lib_path, _ = get_force_field_paths(self.force_field)
        self.ff_lib_name = Path(lib_path).name
        shutil.copy(lib_path, self.inputfiles / self.ff_lib_name)
        for cofactor in self.cofactors:
            shutil.copy(self.workdir / f"{cofactor}.lib", self.inputfiles / f"{cofactor}.lib")

    def read_pdb(self) -> None:
        """Build the complex with the mutant side chain grafted onto the wild type.

        The mutant's side-chain atoms are lower-cased and inserted straight after
        the wild-type residue's backbone oxygen, so both topologies live in one
        residue that ``qprep`` reads through a single library entry.
        """
        mutant_side_chain = []
        for line in self.mutant_pdb_path.read_text().splitlines():
            if not line.startswith("ATOM"):
                continue
            # Mutant residues usually come from PyMOL, which numbers branched
            # hydrogens the other way round from the force field libraries.
            atom = pdb_parse_in(amino_acids.normalize_pdb_atom_name(line))
            if atom[2] in self.backbone:
                continue
            mutant_side_chain.append(atom)

        if not mutant_side_chain:
            raise MutationError(
                f"{self.mutant_pdb_path.name} holds no atoms outside the backbone "
                f"{self.backbone}. It should be the mutant residue only, aligned onto the "
                "wild-type backbone."
            )

        pdb_files = [self.protein_pdb] + [self.workdir / f"{c}.pdb" for c in self.cofactors]
        serial = itertools.count(1)
        residue_counter = itertools.count(1)
        residue_numbers: dict[tuple, int] = {}

        for pdb_file in pdb_files:
            for line in pdb_file.read_text().splitlines():
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                atom = pdb_parse_in(line)

                # Crystallographic waters belong to the prepared solvent sphere,
                # not to the solute complex. They are supplied again through
                # water.pdb below, so qprep can discard any that clash with the
                # hybrid side chain just like it does for grid solvent.
                if atom[4] in _SOLVENT_RESIDUES:
                    continue

                # Q numbers residues positionally, so the complex is renumbered as
                # it is read rather than trusting the numbers in the file.
                key = (pdb_file.name, atom[6], atom[5])
                if key not in residue_numbers:
                    residue_numbers[key] = next(residue_counter)
                residue = residue_numbers[key]
                atom[6] = residue
                atom[1] = next(serial)

                if residue == self.q_position:
                    atom[4] = self.hybrid_name

                self.pdb.setdefault(residue, []).append(atom)
                self.system_size += 1

                if residue == self.q_position and atom[2] == "O":
                    for mutant_atom in mutant_side_chain:
                        mutant_atom[1] = next(serial)
                        mutant_atom[2] = mutant_atom[2].lower()
                        mutant_atom[4] = self.hybrid_name
                        mutant_atom[6] = residue
                        self.pdb[residue].append(mutant_atom)
                        self.system_size += 1

        if self.q_position not in self.pdb:
            raise MutationError(
                f"Q residue {self.q_position} was not found in {self.protein_pdb.name}"
            )

    # ------------------------------------------------------------------
    # Hybrid residue parameters
    # ------------------------------------------------------------------

    def write_hybrid_lib(self) -> None:
        """Write the library entry for the hybrid residue.

        The wild-type entry is taken whole and the mutant's side chain appended
        with lower-cased names, minus the backbone definitions the two share. A
        bond from CA to the mutant's first side-chain atom attaches the second
        topology to the shared backbone.
        """
        headers = ["[atoms]", "[bonds]", "[connections]", "[impropers]", "[info]"]
        if Path(self.force_field).name.startswith("OPLS"):
            headers.append("[charge_groups]")

        lib_path, _ = get_force_field_paths(self.force_field)
        entries = self._read_library_entries(Path(lib_path), {self.wild_type, self.mutant})
        for residue in (self.wild_type, self.mutant):
            if residue not in entries:
                raise MutationError(f"{residue} has no entry in {Path(lib_path).name}")

        header = None
        for line in entries[self.wild_type]:
            # Glycine's alpha carbon carries its own type; the hybrid always has a
            # side chain, so it takes the general one.
            line = line.replace("C223", "C224")
            if any(line.lstrip().startswith(candidate) for candidate in headers):
                header = line.strip().split()[0]
                self.merged_lib[header] = []
            elif header is not None:
                self.merged_lib[header].append(line)

        self.merged_lib.setdefault("[bonds]", [])
        if self.from_gly:
            self.merged_lib["[bonds]"].append("\tCA    ha\n")
        if self.to_gly:
            self.merged_lib["[bonds]"].extend(["\tCA    ha2\n", "\tCA    ha3\n"])
        else:
            self.merged_lib["[bonds]"].append("\tCA    cb\n")

        lowercase: dict[str, str] = {}
        header = None
        for line in entries[self.mutant]:
            if any(line.lstrip().startswith(candidate) for candidate in headers):
                header = line.strip().split()[0]
                continue
            tokens = line.strip().split()
            if not tokens:
                continue
            if header == "[atoms]" and len(tokens) > 1 and tokens[1] not in self.backbone:
                lowercase[tokens[1]] = tokens[1].lower()
            if header == "[connections]":  # the hybrid joins the chain only once
                continue
            if len(tokens) > 1 and (tokens[0] in self.backbone or tokens[1] in self.backbone):
                continue
            self.merged_lib.setdefault(header, []).append(replace(line, lowercase))

        # Q atoms are handled entirely by the FEP file, so the hybrid residue
        # declares no charge groups of its own.
        self.merged_lib["[charge_groups]"] = []

        self.hybrid_lib_name = f"{self.hybrid_name}.lib"
        with open(self.inputfiles / self.hybrid_lib_name, "w") as out:
            out.write(f"{{{self.hybrid_name}}}\n")
            for header in headers:
                out.write(f"\n    {header}\n")
                for index, line in enumerate(self.merged_lib.get(header, []), 1):
                    if header == "[atoms]":
                        tokens = line.split()
                        out.write(f"    {index:2d} {tokens[1]:<10}{tokens[2]:<10}{tokens[3]:<10}\n")
                    else:
                        out.write(line)

    @staticmethod
    def _read_library_entries(lib_path: Path, wanted: set[str]) -> dict[str, list[str]]:
        """Return the raw lines of the named residues in a Q library file."""
        entries: dict[str, list[str]] = {}
        current = None
        for line in lib_path.read_text().splitlines(keepends=True):
            if line.startswith("*"):
                continue
            if line.startswith("{"):
                current = line[1 : line.index("}")]
                continue
            if current in wanted:
                entries.setdefault(current, []).append(line)
        return entries

    def find_missing_bonded(self) -> None:
        """Give a zero force constant to the bonded terms that span the topologies.

        Bonds, angles and torsions connecting the wild-type side chain to the
        mutant one have no physical parameters and must not exert force. Rather
        than enumerate them by hand, ``qprep`` is run once on a throwaway topology
        and asked what it is missing; each term it reports is then defined with a
        zero force constant.
        """
        parameters = read_prm([str(path) for path in self._parameter_files()])

        with open(self.inputfiles / "tmp.pdb", "w") as out:
            for atoms in self.pdb.values():
                for atom in atoms:
                    out.write(f"{pdb_parse_out(atom)}\n")

        with open(self.inputfiles / "tmp.prm", "w") as out:
            for section in _PARAMETER_SECTIONS:
                out.write(f"{section}\n")
                out.writelines(parameters[section])

        (self.inputfiles / "tmp.inp").write_text(
            "\n".join(
                [
                    f"rl {self.hybrid_lib_name}",
                    f"rl {self.ff_lib_name}",
                    "rprm tmp.prm",
                    "rp tmp.pdb",
                    f"boundary 1 {' '.join(str(c) for c in self.center)} {self.radius}",
                    "maketop tmp.top",
                    "q",
                ]
            )
            + "\n"
        )

        self._run_qprep("tmp.inp", "tmp.out", check=False)

        self.zero_force = {"[bonds]": [], "[angles]": [], "[torsions]": [], "[impropers]": []}
        for line in (self.inputfiles / "tmp.out").read_text().splitlines():
            tokens = line.split()
            if len(tokens) < 3 or tokens[1] != "Missing":
                continue
            if tokens[2] == "bond":
                target = "[bonds]"
                entry = f"{tokens[-2]:11}{tokens[-1]:11}{0.0: 8.2f}{0.0:>12.7}\n"
            elif tokens[2] == "angle":
                target = "[angles]"
                entry = f"{tokens[-3]:11}{tokens[-2]:11}{tokens[-1]:11}{0.0: 8.2f}{110.0:>12.7}\n"
            elif tokens[2] == "torsion":
                target = "[torsions]"
                entry = (
                    f"{tokens[-4]:11}{tokens[-3]:11}{tokens[-2]:11}{tokens[-1]:11}"
                    f"{0.0:<10.3f}{1:2d}.000{'0.000':>10}{'1':>10}\n"
                )
            else:
                continue
            if entry not in self.zero_force[target]:
                self.zero_force[target].append(entry)

        counts = {key: len(value) for key, value in self.zero_force.items() if value}
        logger.info(f"Zero-force terms added for the hybrid residue: {counts or 'none'}")

        for leftover in glob.glob(str(self.inputfiles / "tmp*")):
            os.remove(leftover)

    def _parameter_files(self) -> list[Path]:
        """Return the parameter files the hybrid topology is built from."""
        _, prm_path = get_force_field_paths(self.force_field)
        return [Path(prm_path)] + [self.workdir / f"{c}.prm" for c in self.cofactors]

    def write_merged_prm(self) -> None:
        """Write the parameter file, including the hybrid residue's zero-force terms."""
        parameters = read_prm([str(path) for path in self._parameter_files()])

        self.prm_name = f"{Path(self.force_field).name}_merged.prm"
        with open(self.inputfiles / self.prm_name, "w") as out:
            for section in _PARAMETER_SECTIONS:
                out.write(f"{section}\n")
                out.writelines(parameters[section])
                if self.zero_force.get(section):
                    out.write(f"! Zero order {section.strip('[]')} dual topology\n")
                    out.writelines(self.zero_force[section])
                    out.write("\n")

    # ------------------------------------------------------------------
    # Coordinates
    # ------------------------------------------------------------------

    def write_complex(self) -> None:
        """Write the coordinates ``qprep`` builds the topology from."""
        if self.system == "protein":
            with open(self.inputfiles / "complex.pdb", "w") as out:
                for atoms in self.pdb.values():
                    for atom in atoms:
                        out.write(f"{pdb_parse_out(atom)}\n")
        else:
            self._write_reference_peptide()

    def _write_reference_peptide(self) -> None:
        """Cut the mutated residue out and cap it into a small reference peptide.

        The reference leg reproduces the same mutation with no protein around it,
        so the residue is extracted -- with its native neighbours or with
        synthetic Ala/Gly flanks -- and capped with ACE and NMA.
        """
        if self.tripeptide_flanks == "Z":  # keep the native neighbours
            source_residues = [self.q_position - 1, self.q_position, self.q_position + 1]
        else:
            source_residues = [self.q_position]

        missing = [residue for residue in source_residues if residue not in self.pdb]
        if missing:
            raise MutationError(
                f"Cannot build a native-sequence reference peptide: Q residue(s) {missing} "
                "are not in the prepared sphere. Use -t A, -t G or -t X instead."
            )

        try:
            peptide = build_reference_peptide(
                [self.pdb[residue] for residue in source_residues], self.tripeptide_flanks
            )
        except PeptideBuildError as error:
            raise MutationError(str(error)) from error

        serial = itertools.count(1)
        for atom in peptide:
            atom[1] = next(serial)

        self.system_size = len(peptide)
        complex_pdb = self.inputfiles / "complex.pdb"
        with open(complex_pdb, "w") as out:
            for atom in peptide:
                out.write(f"{pdb_parse_out(atom)}\n")
            self.counter_ions = self._append_counter_ions(out, len(peptide))

    def _append_counter_ions(self, out, first_free_serial: int) -> int:
        """Charge-match the reference sphere to the protein sphere.

        A mutation that changes the net charge leaves the two legs at different
        total charge, and that difference does not cancel in the thermodynamic
        cycle. Ions are placed on a shell inside the boundary to bring the
        reference sphere to the charge the protein sphere carries.

        Returns:
            The number of ions written.
        """
        if not (CHARGED_RESIDUES & {self.wild_type, self.mutant}):
            return 0

        from .counter_ions import minimize_coulomb_on_sphere

        charge = self.prep.total_charge
        if self.wild_type in ("ASP", "GLU"):
            charge += 1
        elif self.wild_type in ("ARG", "LYS", "HIP"):
            charge -= 1

        if charge == 0:
            return 0

        # Kept clear of the boundary shell, where the radial restraints act.
        placement_radius = int(self.radius) - 11
        positions = minimize_coulomb_on_sphere(abs(charge), placement_radius, self.center)
        ion = "CLA" if charge < 0 else "SOD"
        serial = first_free_serial
        for residue, (x, y, z) in enumerate(positions, start=6):
            serial += 1
            out.write(
                f"ATOM   {serial:4}  {ion} {ion}   {residue:3}     "
                f"{x:7.3f} {y:7.3f} {z:7.3f}  0.00  0.00           \n"
            )
        logger.info(f"Added {len(positions)} {ion} ion(s) to match the protein sphere charge")
        return len(positions)

    def write_water(self) -> None:
        """Copy the prepared water sphere, dropping molecules that clash with cofactors."""
        cofactor_coordinates = []
        for cofactor in self.cofactors:
            for line in (self.workdir / f"{cofactor}.pdb").read_text().splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    atom = pdb_parse_in(line)
                    cofactor_coordinates.append([atom[8], atom[9], atom[10]])

        clashing: set[int] = set()
        waters: dict[int, list[str]] = {}
        for line in (self.workdir / "water.pdb").read_text().splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue  # qprep rejects TITLE lines
            # Read the columns directly: qprep writes water.pdb without the
            # occupancy and B-factor fields, so a full PDB parse would fail.
            residue = int(line[22:26])
            waters.setdefault(residue, []).append(line)
            if not cofactor_coordinates:
                continue
            position = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if any(euclidian_overlap(position, other, 1.6) for other in cofactor_coordinates):
                clashing.add(residue)

        with open(self.inputfiles / "water.pdb", "w") as out:
            out.write(f"{self.radius:.1f} SPHERE\n")
            for residue, lines in waters.items():
                if residue not in clashing:
                    out.write("\n".join(lines) + "\n")

        if clashing:
            logger.info(f"Removed {len(clashing)} water molecule(s) clashing with cofactors")

    # ------------------------------------------------------------------
    # Topology
    # ------------------------------------------------------------------

    def write_qprep(self) -> None:
        """Write the qprep input that builds the hybrid topology."""
        libraries = [self.hybrid_lib_name, self.ff_lib_name]
        libraries += [f"{cofactor}.lib" for cofactor in self.cofactors]

        # The protein leg reuses the water sphere built during preparation; the
        # reference peptide is small enough that qprep solvates it from scratch.
        solvent = "4 water.pdb" if self.system == "protein" else "1 HOH"

        cysbonds = ""
        if self.system == "protein":
            cysbonds = "".join(
                f"addbond {first}:SG {second}:SG y\n" for first, second in self.prep.disulfides
            )

        density = sphere_solute_density(self.protein_pdb, self.center, self.radius)
        params = QprepResFEPParameters(
            libraries=libraries,
            prm=self.prm_name,
            pdb="complex.pdb",
            center=" ".join(str(c) for c in self.center),
            sphere_radius=f"{self.radius:.1f}",
            solute_density=f"{density:.5f}",
            solvent=solvent,
            cysbonds=cysbonds,
        )
        (self.inputfiles / "qprep.inp").write_text(render_qprep_resfep_input(params))

    def run_qprep(self) -> None:
        """Build the topology and read back the atom numbering Q assigned."""
        self._run_qprep("qprep.inp", "qprep.out")
        self._read_topology_atom_ids()
        self._read_topology_size()

    def _read_topology_size(self) -> None:
        """Take the solute size from the topology rather than from the PDB.

        ``qprep`` builds atoms the input PDB does not carry -- the proton of a
        residue prepared in a neutral form (``GLU`` as ``GLH``), a missing
        N-terminal hydrogen -- so the complex PDB undercounts the solute. The
        whole-solute restraint of eq1--eq4 spans ``1..system_size``, and anything
        past that count moves unrestrained while the rest of the protein is held.

        Counter-ions are not part of it: they sit at the end of the solute, are
        held by their own wall restraint, and :meth:`_wall_restraints` addresses
        them as the atoms just past ``system_size``.
        """
        solute = 0
        for line in (self.inputfiles / "top_p.pdb").read_text().splitlines():
            if line.startswith("ATOM") and line[17:20].strip() not in _SOLVENT_RESIDUES:
                solute += 1

        size = solute - self.counter_ions
        if size != self.system_size:
            logger.debug(
                f"qprep built {solute - self.system_size - self.counter_ions} atom(s) the "
                f"input PDB did not carry; restraining {size} solute atoms, not "
                f"{self.system_size}"
            )
        self.system_size = size

    def _run_qprep(self, input_file: str, output_file: str, check: bool = True) -> None:
        qprep_path = CLUSTER_DICT[self.cluster]["QPREP"]
        cwd = Path.cwd()
        try:
            os.chdir(self.inputfiles)
            if check:
                run_qprep(qprep_path, input_file, output_file, self.force_field)
            else:
                # The probe run is *expected* to report missing parameters, so
                # run_qprep's error check would reject it. qprep reads its
                # directives from stdin, hence the redirection.
                with open(input_file) as stdin, open(output_file, "w") as stdout:
                    subprocess.run([qprep_path], stdin=stdin, stdout=stdout, stderr=subprocess.STDOUT)
        finally:
            os.chdir(cwd)

    def _read_topology_atom_ids(self) -> None:
        """Record the topology atom index of every hybrid-residue atom.

        Everything downstream -- the FEP files, the distance restraints, the
        reference peptide's anchor -- addresses atoms by their index in the
        topology, which only ``qprep`` can assign. The hybrid residue is found by
        its ``<WT>2<MUT>`` name.
        """
        hybrid = self.pdb[self.q_position]
        position = itertools.count()
        self.anchor_atom = None

        for line in (self.inputfiles / "top_p.pdb").read_text().splitlines():
            if not line.startswith("ATOM"):
                continue
            tokens = line.split()
            serial, name, residue_name = tokens[1], tokens[2], tokens[3]
            if residue_name != self.hybrid_name:
                continue

            index = next(position)
            if index < len(hybrid):
                hybrid[index][1] = int(serial)
                hybrid[index][2] = name
                hybrid[index][8] = tokens[5]
                hybrid[index][9] = tokens[6]
                hybrid[index][10] = tokens[7]
            if name.upper() in amino_acids.SIDE_CHAIN_ATOM_NAMES:
                self.atom_ids[name] = serial
            if name == "CA":
                self.anchor_atom = int(serial)

        if not self.atom_ids:
            raise MutationError(
                f"No {self.hybrid_name} residue found in top_p.pdb. Check inputfiles/qprep.out."
            )

    # ------------------------------------------------------------------
    # Simulation inputs
    # ------------------------------------------------------------------

    def lambda_schedule(self) -> list[list[str]]:
        """Return the lambda values of each stage, in the order they are run."""
        if self.sampling == "exponential":
            # Each stage is steep at the end where its own topology switches, so
            # the two get mirrored spacing.
            return [
                resfep_lambda_ladder(self.windows, "exponential"),
                resfep_lambda_ladder(self.windows, "reverse_exponential"),
            ]
        schedule = resfep_lambda_ladder(self.windows, self.sampling)
        return [list(schedule), list(schedule)]

    def _distance_restraints(self) -> str:
        """Restrain equivalent heavy atoms of the two side chains onto each other."""
        wild_type_atoms, mutant_atoms = amino_acids.side_chain_pair(self.wild_type, self.mutant)
        pairs = amino_acids.match_heavy_atoms(
            wild_type_atoms, mutant_atoms, self.pdb[self.q_position]
        )
        lines = [
            f"{self.atom_ids[wt]} {self.atom_ids[mut]} "
            f"{RESTRAINT_LOW} {RESTRAINT_HIGH} {RESTRAINT_FORCE} 0"
            for wt, mut in pairs
            if wt in self.atom_ids and mut in self.atom_ids
        ]
        logger.info(f"{len(lines)} distance restraint(s) between the two side chains")
        return "\n".join(lines)

    def _anchor_restraint(self) -> str:
        """Hold the reference peptide at the centre of its sphere.

        With no protein around it the capped peptide would drift to the boundary,
        so its alpha carbon is pinned. The protein leg needs no such anchor.
        """
        if self.system == "protein" or self.anchor_atom is None:
            return ""
        return f"{self.anchor_atom} {self.anchor_atom} 1.0 0 0"

    def _wall_restraints(self) -> str:
        """Keep the charge-matching ions away from the boundary shell."""
        if not self.counter_ions:
            return ""
        first = self.system_size + 1
        return format_wall_restraints(
            first, first + self.counter_ions - 1, radius=int(self.radius) - 5, force=1.0
        )

    def _equilibration_lambdas(self) -> tuple[str, str]:
        """Return the lambda pair the equilibration runs at.

        Equilibration sits at the endpoint the first stage starts from. A mutation
        out of glycine runs in the opposite direction, since the transformation
        has to begin from the larger side chain.
        """
        if self.from_gly or self.start == "0.5":
            return "0.000", "1.000"
        return "1.000", "0.000"

    def write_equilibration(self) -> None:
        """Write eq1 through eq5."""
        distance_restraints = self._distance_restraints()
        anchor = self._anchor_restraint()
        wall = self._wall_restraints()
        lambda1, lambda2 = self._equilibration_lambdas()

        configs = get_resfep_equilibration_configs(
            self.timestep,
            self.shell_radius,
            self.eq5_steps,
            separate_scaling=self.separate_scaling,
        )
        for index, config in enumerate(configs):
            if config.use_water_restraint:
                sequence_restraints = anchor
            else:
                sequence_restraints = format_sequence_restraint(
                    1, self.system_size, config.sequence_restraint_force
                )
            content = render_md_input(
                params=config.params,
                lambda1=lambda1,
                lambda2=lambda2,
                trajectory_file=f"{config.name}.dcd" if self.write_trajectories else None,
                final_file=f"{config.name}.re",
                restart_file=None if index == 0 else f"{configs[index - 1].name}.re",
                distance_restraints=distance_restraints,
                sequence_restraints=sequence_restraints,
                wall_restraints=wall,
                is_eq1=(index == 0),
            )
            (self.inputfiles / f"{config.name}.inp").write_text(content)

    def write_production(self) -> list[list[str]]:
        """Write the production MD inputs of both stages.

        Returns:
            The MD basenames of each stage, in run order.
        """
        distance_restraints = self._distance_restraints()
        anchor = self._anchor_restraint()
        wall = self._wall_restraints()
        config = get_production_config(
            self.timestep,
            self.shell_radius,
            steps=self.production_steps,
            interval_output=5,
            separate_scaling=self.separate_scaling,
        )

        stage_files = []
        for stage, lambdas in enumerate(self.lambda_schedule(), start=1):
            names = []
            previous = "eq5.re"  # stage 2 inherits stage 1's endpoint under this name
            for window in range(self.windows):
                if self.start == "0.5" and stage == 1:
                    value = 1.000 - float(lambdas[-window - 1])
                else:
                    value = float(lambdas[window])

                first, second = f"{value:.3f}", f"{1.000 - value:.3f}"
                basename = f"md_{stage}_{first}_{second}".replace(".", "")
                # Glycine mutations and mid-lambda starts run the transformation the
                # other way round; the file name keeps the canonical order so the
                # analysis and the run script can still find the windows.
                if self.from_gly or (self.start == "0.5" and stage == 1):
                    first, second = second, first

                content = render_md_input(
                    params=config.params,
                    lambda1=first,
                    lambda2=second,
                    trajectory_file=f"{basename}.dcd" if self.write_trajectories else None,
                    final_file=f"{basename}.re",
                    restart_file=previous,
                    energy_file=f"{basename}.en",
                    distance_restraints=distance_restraints,
                    sequence_restraints=anchor,
                    wall_restraints=wall,
                )
                (self.inputfiles / f"{basename}.inp").write_text(content)
                names.append(basename)
                previous = f"{basename}.re"
            stage_files.append(names)

        return stage_files

    def write_qfep(self) -> None:
        """Write the qfep header; the run script appends each stage's energy files."""
        content = render_qfep_input(
            total_lambdas=self.windows,
            temperature=float(self.temperature.split(",")[0]),
            windows=self.windows,
            energy_files=[],
        )
        (self.inputfiles / "qfep.inp").write_text(content + "\n")

    # ------------------------------------------------------------------
    # FEP files
    # ------------------------------------------------------------------

    def write_fep_files(self) -> None:
        """Write FEP1.fep and FEP2.fep.

        A mutation out of glycine has its two stages swapped, so that the
        transformation always begins from the larger side chain.
        """
        _, prm_path = get_force_field_paths(self.force_field)
        vdw_parameters = read_prm([str(prm_path)])["[atom_types]"]

        filenames = list(FEP_STAGES)
        if self.from_gly:
            filenames.reverse()

        for stage_name, filename in zip(FEP_STAGES, filenames):
            content = self._render_fep_file(stage_name == "FEP1.fep", vdw_parameters)
            (self.inputfiles / filename).write_text(content)

    def _render_fep_file(self, is_first_stage: bool, vdw_parameters: list[str]) -> str:
        hybrid = self.pdb[self.q_position]
        lines = [
            f"! Dual topology QresFEP: {self.hybrid_name}",
            "",
            "[FEP]",
            "   states 2",
            "   softcore_use_max_potential on",
            "",
            "[atoms]",
        ]
        for index, atom in enumerate(hybrid, start=1):
            lines.append(f"{index:4d} {atom[1]:4d} ! {atom[2]:<4s}")

        lines += ["", "[change_charges]"]
        for index, entry in enumerate(self.merged_lib["[atoms]"], start=1):
            tokens = entry.split()
            first, second = self._charge_pair(tokens[1], tokens[3], is_first_stage)
            lines.append(f"{index:4d}{first:>10s}{second:>10s}")

        lines += ["", "[atom_types]"]
        for entry in vdw_parameters:
            tokens = entry.split()
            if len(tokens) >= 7:
                lines.append(
                    f"{tokens[0]:<8s}{tokens[1]:>10s}{tokens[3]:>10s}  0.0  0.0"
                    f"{tokens[4]:>10s}{tokens[5]:>10s}{tokens[6]:>10s}"
                )

        lines += ["", "[change_atoms]"]
        for index, entry in enumerate(self.merged_lib["[atoms]"], start=1):
            tokens = entry.split()
            first, second = self._type_pair(tokens[1], tokens[2], is_first_stage)
            lines.append(f"{index:4d}{first:>10s}{second:>10s}")

        lines += ["", "[softcore]"]
        for index, entry in enumerate(self.merged_lib["[atoms]"], start=1):
            name = entry.split()[1]
            if name in self.backbone:
                pair = ("0", "0")  # the shared backbone is never softcored
            elif is_first_stage:
                pair = ("0", SOFTCORE_MAX_POTENTIAL)
            else:
                pair = (SOFTCORE_MAX_POTENTIAL, "0")
            lines.append(f"{index:4d}{pair[0]:>10s}{pair[1]:>10s}")

        # The two side chains occupy the same space and must not see each other.
        lines += ["", "[excluded_pairs]"]
        side_chain = [(atom[1], atom[2]) for atom in hybrid if atom[2] not in self.backbone]
        wild_type_serials = [serial for serial, name in side_chain if not name.islower()]
        mutant_serials = [serial for serial, name in side_chain if name.islower()]
        for first_serial in wild_type_serials:
            for second_serial in mutant_serials:
                lines.append(f"{first_serial:6d}{second_serial:6d}{'1':>5s}{'1':>5s}")

        # No bond spans the two topologies, so nothing needs zeroing here.
        lines += ["", "[bond_types]", "", "[change_bonds]", ""]
        lines += ["[angle_types]", "", "[change_angles]"] + self._zero_angles()
        lines += ["", "[torsion_types]", "", "[change_torsions]"] + self._zero_torsions()

        return "\n".join(lines) + "\n"

    def _charge_pair(self, name: str, charge: str, is_first_stage: bool) -> tuple[str, str]:
        """Return one hybrid atom's charge in both states of a stage."""
        if name.islower():  # mutant side chain: charged in stage 2
            if self.from_gly and name == "ha":
                return ("0.000", _HA_CHARGE) if is_first_stage else (_HA_CHARGE, _HA_CHARGE)
            return ("0.000", "0.000") if is_first_stage else ("0.000", charge)

        # wild-type side chain: discharged in stage 1
        if self.to_gly and name == "HA":
            return (_HA_CHARGE, _HA_CHARGE) if is_first_stage else (_HA_CHARGE, "0.000")
        if name in self.backbone:
            if name == "CA" and self.from_gly:
                return (
                    (_CA_CHARGE_GLY, _CA_CHARGE_NON_GLY)
                    if is_first_stage
                    else (_CA_CHARGE_NON_GLY, _CA_CHARGE_NON_GLY)
                )
            if name == "CA" and self.to_gly:
                return (
                    (_CA_CHARGE_NON_GLY, _CA_CHARGE_NON_GLY)
                    if is_first_stage
                    else (_CA_CHARGE_NON_GLY, _CA_CHARGE_GLY)
                )
            return (charge, charge)  # the backbone is shared and never perturbed
        return (charge, "0.000") if is_first_stage else ("0.000", "0.000")

    def _type_pair(self, name: str, atom_type: str, is_first_stage: bool) -> tuple[str, str]:
        """Return one hybrid atom's type in both states of a stage."""
        if name.islower():  # mutant side chain: grown in during stage 2
            if self.from_gly and name == "ha":
                return ("DUM", _HA_TYPE) if is_first_stage else (_HA_TYPE, _HA_TYPE)
            return ("DUM", atom_type) if is_first_stage else (atom_type, atom_type)

        if self.to_gly and name == "HA":
            return (_HA_TYPE, _HA_TYPE) if is_first_stage else (_HA_TYPE, "DUM")
        if name in self.backbone:
            return (atom_type, atom_type)
        return (atom_type, atom_type) if is_first_stage else (atom_type, "DUM")

    def _zero_angles(self) -> list[str]:
        """Zero the angles that would couple the two side chains through CA."""
        if self.from_gly:
            angles = [
                ("HA3", "CA", "cb"), ("HA2", "CA", "cb"),
                ("HA3", "CA", "ha"), ("HA2", "CA", "ha"),
            ]
        elif self.to_gly:
            angles = [
                ("ha3", "CA", "CB"), ("ha2", "CA", "CB"),
                ("ha3", "CA", "HA"), ("ha2", "CA", "HA"),
            ]
        else:
            angles = [("CB", "CA", "cb")]

        lines = []
        for names in angles:
            if not all(name in self.atom_ids for name in names):
                logger.debug(f"Skipping zero-angle {names}: not every atom is in the topology")
                continue
            ids = [self.atom_ids[name] for name in names]
            lines.append(f"{ids[0]:>6s}{ids[1]:>6s}{ids[2]:>6s}{0:>5d}{0:>5d}")
        return lines

    def _zero_torsions(self) -> list[str]:
        """Zero the chi1 torsions that would couple the two side chains."""
        lines = []
        for residue in (self.wild_type, self.mutant):
            atoms = amino_acids.chi1_atoms(residue)
            if not atoms:  # glycine has no chi1
                continue

            is_mutant = residue == self.mutant
            atoms = [atom.lower() if is_mutant else atom.upper() for atom in atoms]
            core = ["cb", "CA"] if is_mutant else ["CB", "CA"]
            if self.from_gly and is_mutant:
                dummies = ["HA2", "HA3"]
            elif self.to_gly and not is_mutant:
                dummies = ["ha2", "ha3"]
            elif is_mutant:
                dummies = ["CB"]
            else:
                dummies = ["cb"]

            for atom in atoms:
                for dummy in dummies:
                    names = [atom, core[0], core[1], dummy]
                    if not all(name in self.atom_ids for name in names):
                        logger.debug(
                            f"Skipping zero-torsion {names}: not every atom is in the topology"
                        )
                        continue
                    ids = [self.atom_ids[name] for name in names]
                    lines.append(f"{ids[0]:>6s}{ids[1]:>6s}{ids[2]:>6s}{ids[3]:>6s}{0:>4d}{0:>4d}")
        return lines

    # ------------------------------------------------------------------
    # Job scripts
    # ------------------------------------------------------------------

    def write_runfile(self, stage_files: list[list[str]]) -> None:
        """Write the cluster run script driving both FEP stages.

        Raises:
            MutationError: If asked for a local profile. The two stages have to be
                chained through a job that carries stage 1's endpoint into stage 2,
                which only the SLURM script does; emitting it for LOCAL would
                produce a script that cannot run.
        """
        if self.cluster in _LOCAL_CLUSTERS:
            raise MutationError(
                "QresFEP has no local run script: the two FEP stages must be chained, "
                "which only the cluster script does. Pick a cluster profile from "
                "settings.py, or run the generated eq*.inp / md_*.inp by hand."
            )

        template = Path(CONFIGS["INPUT_DIR"]) / "run_resfep.sh"
        target = self.inputfiles / f"run{self.cluster}.sh"

        substitutions = dict(CLUSTER_DICT[self.cluster])
        substitutions["FEPS"] = " ".join(FEP_STAGES)
        substitutions["JOBNAME"] = f"{self.system[0]}_{self.name}"
        substitutions["TOTAL_JOBS"] = str(self.replicates * len(self.temperature.split(",")))
        substitutions["TEMP_VAR"] = self.temperature.replace(",", " ")
        substitutions["RANDOM_SEEDS"] = " ".join(str(seed) for seed in self.seeds)
        # Stage 2 restarts from stage 1's last window.
        substitutions["RESTARTFILE"] = f"{stage_files[0][-1]}.re"

        mpirun = "time mpirun -n $SLURM_NTASKS --map-by core --bind-to core $qdyn"
        equilibration = [
            Path(path).stem for path in sorted(glob.glob(str(self.inputfiles / "eq*.inp")))
        ]

        with open(template) as source, open(target, "w") as out:
            for line in source:
                if line.strip() == "#SBATCH -A ACCOUNT" and "ACCOUNT" not in substitutions:
                    continue
                out.write(replace(line, substitutions))

                if line.strip() == "#EQ_FILES":
                    for name in equilibration:
                        out.write(f"{mpirun} {name}.inp > {name}.log\n")
                for marker, stage in (("#RUN_1_FILES", 0), ("#RUN_2_FILES", 1)):
                    if line.strip() == marker:
                        for name in stage_files[stage]:
                            out.write(f"echo {name}\n{mpirun} {name}.inp > {name}.log\n")
                if line.strip() == "#CLEANUP" and self.to_clean:
                    out.write(
                        "rm -f -- "
                        + " ".join(f"*{suffix}" for suffix in self.to_clean)
                        + "\n"
                    )

        target.chmod(target.stat().st_mode | stat.S_IEXEC)

    def write_submitfile(self) -> None:
        """Write the FEP_submit.sh wrapper that queues the run script."""
        if self.cluster in _LOCAL_CLUSTERS:
            return
        template = Path(CONFIGS["ROOT_DIR"]) / "INPUTS" / "FEP_submit.sh"
        target = self.directory / "FEP_submit.sh"
        target.write_text(replace(template.read_text(), {"RUNFILE": f"run{self.cluster}.sh"}))
        target.chmod(target.stat().st_mode | stat.S_IEXEC)

    # ------------------------------------------------------------------

    def run(self) -> Path:
        """Generate every input file for this mutation.

        Returns:
            The FEP directory that was created.
        """
        self.read_prep()
        self.create_environment()
        self.read_pdb()
        self.write_hybrid_lib()
        self.find_missing_bonded()
        self.write_merged_prm()
        self.write_complex()
        self.write_water()
        self.write_qprep()
        self.run_qprep()
        self.write_equilibration()
        stage_files = self.write_production()
        self.write_fep_files()
        self.write_qfep()
        self.write_runfile(stage_files)
        self.write_submitfile()
        logger.info(f"Wrote {self.directory}")
        return self.directory
