"""Command-line handling of execution profiles and simulation protocols.

Precedence, lowest to highest: packaged defaults, the user file, the project file, and
finally flags given on the command line. Protocol-governed flags therefore default to
``None`` in the parsers and are filled in here from the selected protocol.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys

from ..config import ConfigError, Profile, Protocol, load_config
from ..logger import logger

#: NEQ stage lengths (MD steps) used when neither a flag nor the protocol sets them.
NEQ_STEP_DEFAULTS = {"neq_steps": 50000, "neq_eq_steps": 1000, "neq_relax_steps": 5000}
#: Flag -> protocol stage whose length it sets. An explicit flag beats the protocol.
NEQ_STEP_STAGES = {"neq_steps": "switch", "neq_eq_steps": "spacing", "neq_relax_steps": "relax"}


def add_config_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the profile/protocol selection flags shared by ``qligfep`` and ``setupFEP``."""
    parser.add_argument(
        "-c",
        "--cluster",
        "--profile",
        dest="cluster",
        default=None,
        help=(
            "Execution profile: where and how the jobs run (Slurm array or local, engine, modules, "
            "cores). Defaults to `default_profile` from your configuration. List the available "
            "profiles with `qligfep_config show`; create your own with `qligfep_config init`."
        ),
    )
    parser.add_argument(
        "--protocol",
        dest="protocol",
        default=None,
        help=(
            "Simulation protocol from the project's qligfep.yaml: method, lambda windows, stage "
            "lengths and qdyn input overrides. Flags given on the command line take precedence over "
            "the protocol. Defaults to `default_protocol` (the packaged `default`)."
        ),
    )
    parser.add_argument(
        "--config",
        dest="config",
        default=None,
        help="Project configuration file. Defaults to qligfep.yaml in the working directory.",
    )
    parser.add_argument(
        "--account",
        dest="account",
        default=None,
        help=(
            "Slurm account (project allocation) to charge, overriding the profile's `slurm.account`. "
            "Without either, jobs are submitted without -A and Slurm charges your default account."
        ),
    )
    parser.add_argument(
        "--partition",
        dest="partition",
        default=None,
        help="Slurm partition to submit to, overriding the profile's `slurm.partition`.",
    )
    parser.add_argument(
        "--q-bin",
        dest="q_bin",
        default=None,
        help=(
            "Directory holding the Q binaries (qprep, qdyn, qdynp, qfep) to use instead of the "
            "profile's `q_bin` or the ones compiled inside this installation."
        ),
    )


def _fill(args: argparse.Namespace, name: str, value) -> None:
    if getattr(args, name, None) is None:
        setattr(args, name, value)


def resolve_run_settings(args: argparse.Namespace, system: str | None = None) -> tuple[Profile, Protocol]:
    """Resolve the profile and protocol, and fill unset protocol-governed flags in ``args``.

    Also sets ``args.cluster`` to the resolved profile name, ``args.stage_overrides`` to
    the protocol's stage overrides, and ``args.explicit_neq_steps`` to the NEQ length
    flags that were given on the command line.

    Exits with a readable message on a configuration error.
    """
    try:
        config = load_config(project_file=getattr(args, "config", None))
        profile = config.get_profile(
            getattr(args, "cluster", None), system=system, q_bin=getattr(args, "q_bin", None)
        )
        protocol = config.get_protocol(getattr(args, "protocol", None))
    except ConfigError as error:
        sys.exit(f"Configuration error: {error}")

    account, partition = getattr(args, "account", None), getattr(args, "partition", None)
    if account or partition:
        if profile.is_local:
            logger.warning(f"--account/--partition have no effect on the local profile '{profile.name}'.")
        slurm = dataclasses.replace(
            profile.slurm,
            account=account or profile.slurm.account,
            partition=partition or profile.slurm.partition,
        )
        profile = dataclasses.replace(profile, slurm=slurm)

    args.cluster = profile.name
    args.protocol = protocol.name
    _fill(args, "sphereradius", str(protocol.sphere_radius))
    _fill(args, "start", protocol.start)
    _fill(args, "temperature", protocol.temperature)
    _fill(args, "replicates", str(protocol.replicates))
    _fill(args, "sampling", protocol.sampling)
    _fill(args, "windows", str(protocol.windows))
    _fill(args, "timestep", protocol.timestep)
    _fill(args, "neq", protocol.method == "nonequilibrium")
    _fill(args, "neq_reps", protocol.switch_pairs)
    _fill(args, "neq_L", protocol.steepness)
    _fill(args, "neq_schedule", protocol.schedule)

    stages = dict(protocol.stages)
    args.explicit_neq_steps = [name for name in NEQ_STEP_DEFAULTS if getattr(args, name, None) is not None]
    for name in args.explicit_neq_steps:
        stage = NEQ_STEP_STAGES[name]
        override = stages.get(stage)
        if override is not None and (override.time_fs is not None or override.steps is not None):
            logger.warning(
                f"--{name.replace('_', '-')} overrides the length protocol '{protocol.name}' sets for "
                f"stage '{stage}'."
            )
            stages[stage] = dataclasses.replace(override, time_fs=None, steps=None)
    for name, default in NEQ_STEP_DEFAULTS.items():
        # Report the protocol's stage length in the run record instead of the unused default.
        stage = NEQ_STEP_STAGES[name]
        override = stages.get(stage)
        if getattr(args, name, None) is None and override is not None:
            custom_stepsize = any("stepsize" in keys for sec, keys in override.inp.items() if sec.lower() == "md")
            if not custom_stepsize:
                try:
                    steps = override.n_steps(float(args.timestep[0]), f"protocols.{protocol.name}.stages.{stage}")
                except ConfigError as error:
                    sys.exit(f"Configuration error: {error}")
                if steps is not None:
                    setattr(args, name, steps)
        _fill(args, name, default)
    args.stage_overrides = stages

    method = "nonequilibrium" if args.neq else "equilibrium"
    if method != protocol.method:
        logger.warning(
            f"Running the {method} method although protocol '{protocol.name}' is written for "
            f"{protocol.method}; its {protocol.method}-only stage overrides do not apply."
        )
    return profile, protocol


def resolved_snapshot(profile: Profile, protocol: Protocol, args: argparse.Namespace) -> dict:
    """What was actually used to set up a run, for the provenance files."""
    effective = protocol.as_dict()
    effective.update(
        method="nonequilibrium" if args.neq else "equilibrium",
        timestep=args.timestep,
        temperature=str(args.temperature),
        replicates=int(args.replicates),
        sphere_radius=int(args.sphereradius),
        equilibrium={"windows": int(args.windows), "sampling": args.sampling, "start": str(args.start)},
        nonequilibrium={
            "switch_pairs": args.neq_reps,
            "schedule": args.neq_schedule,
            "steepness": args.neq_L,
            "switch_steps": args.neq_steps,
            "spacing_steps": args.neq_eq_steps,
            "relax_steps": args.neq_relax_steps,
        },
        stages={name: stage.as_dict() for name, stage in args.stage_overrides.items()},
    )
    return {"profile": profile.as_dict(), "protocol": effective}
