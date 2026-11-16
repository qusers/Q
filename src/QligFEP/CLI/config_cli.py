"""``qligfep_config``: create, inspect and check QligFEP configuration.

Subcommands:
    init        write a commented user configuration for this machine
    show        list profiles and protocols, or print one fully resolved
    check       validate the configuration and test that a profile can run here
    build-info  record the toolchain of the Q binaries (for ``modules: from-build``)
"""

from __future__ import annotations

import argparse
import sys

import yaml

from ..config import ConfigError, load_config, user_config_path
from ..config.detect import check_profile, detect_environment, render_user_config, write_build_info
from ..config.schema import resolve_q_bin


def _add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        default=None,
        help="Project configuration file. Defaults to qligfep.yaml in the working directory.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="qligfep_config",
        description=(
            "Create, inspect and check QligFEP configuration. Execution profiles (where and how jobs "
            "run) are read from the packaged defaults, then the user file, then the project's "
            "qligfep.yaml; protocols (what is simulated) from the packaged defaults and the project file."
        ),
    )
    commands = parser.add_subparsers(dest="command", required=True)

    init = commands.add_parser(
        "init",
        help="write a commented user configuration for this machine",
        description=(
            "Detect the scheduler, cores, loaded modules and Slurm accounts, and write a commented "
            "profile to the user configuration file. An existing file is never modified: the proposal "
            "is printed instead."
        ),
    )
    init.add_argument("--name", default=None, help="Profile name. Defaults to `cluster` or `local`.")
    init.add_argument(
        "--execution",
        choices=["auto", "slurm-array", "local"],
        default="auto",
        help="Execution type. `auto` picks slurm-array when sbatch is on PATH.",
    )
    init.add_argument("--ntasks", default=None, help="Cores per replicate (or `auto` for local execution).")
    init.add_argument("--time", default=None, help='Wall time per job, e.g. "0-08:00:00".')
    init.add_argument("--account", default=None, help="Slurm account to charge.")
    init.add_argument("--partition", default=None, help="Slurm partition.")
    init.add_argument(
        "--modules",
        nargs="*",
        default=None,
        help=(
            "Modules to load in jobs. Defaults to the modules recorded at build time when available "
            "(`from-build`), otherwise none."
        ),
    )
    init.add_argument("--q-bin", default=None, help="Directory with the Q binaries to inspect.")
    init.add_argument("--output", default=None, help="Write to this file instead of the user configuration file.")
    init.add_argument("--print", action="store_true", dest="print_only", help="Print the file instead of writing it.")

    show = commands.add_parser(
        "show",
        help="list profiles and protocols, or print one fully resolved",
        description="Without options, list what is available and where it comes from.",
    )
    _add_common(show)
    show.add_argument("--profile", "-c", default=None, help="Print this profile fully resolved.")
    show.add_argument("--protocol", default=None, help="Print this protocol fully resolved.")
    show.add_argument(
        "--system",
        choices=["water", "protein", "vacuum"],
        default=None,
        help="Apply the profile's overrides for this system.",
    )

    check = commands.add_parser(
        "check",
        help="validate the configuration and test that a profile can run here",
        description=(
            "Validate every profile and protocol, then check the selected profile on this machine: "
            "binaries, sbatch, and that the engine's libraries and mpirun are found after the "
            "profile's setup lines and modules have run."
        ),
    )
    _add_common(check)
    check.add_argument("--profile", "-c", default=None, help="Profile to test. Defaults to default_profile.")
    check.add_argument("--q-bin", default=None, help="Directory with the Q binaries to test.")

    build = commands.add_parser(
        "build-info",
        help="record the toolchain of the Q binaries (for `modules: from-build`)",
        description=(
            "Write build-info.json next to the Q binaries with the modules loaded in this shell and "
            "the compiler and MPI versions. Run it in the environment the binaries were compiled in."
        ),
    )
    build.add_argument("--q-bin", default=None, help="Directory with the Q binaries. Defaults to this installation's.")
    build.add_argument("--modules", nargs="*", default=None, help="Record these modules instead of the loaded ones.")
    return parser


def cmd_init(args: argparse.Namespace) -> int:
    q_bin, _ = resolve_q_bin(None, "init", args.q_bin)
    env = detect_environment(q_bin)
    execution = args.execution
    if execution == "auto":
        execution = "slurm-array" if env["has_slurm"] else "local"
    if execution == "slurm-array" and not env["has_slurm"]:
        print("Note: sbatch is not on PATH here; writing a Slurm profile anyway.", file=sys.stderr)
    local = execution == "local"

    if args.modules is not None:
        modules: list[str] | str = args.modules
    elif env["build_info"] is not None:
        modules = "from-build"
    else:
        # The modules of this shell say nothing about what the binaries were built with,
        # so they are reported below rather than written into the profile.
        modules = []
    account = args.account
    if account is None and not local and len(env["accounts"]) == 1:
        account = env["accounts"][0]

    text = render_user_config(
        env,
        name=args.name or ("local" if local else "cluster"),
        execution=execution,
        ntasks=args.ntasks or ("auto" if local else 8),
        time=args.time or ("0-24:00:00" if local else "0-08:00:00"),
        modules=modules,
        account=account,
        partition=args.partition,
    )
    try:  # never hand out a file that does not load
        yaml.safe_load(text)
    except yaml.YAMLError as error:
        print(f"Internal error: the generated configuration is not valid YAML ({error}).", file=sys.stderr)
        return 1

    notes = []
    missing = [name for name, present in env["binaries"].items() if not present]
    if missing:
        notes.append(f"Q binaries not found in {q_bin}: {', '.join(missing)}.")
    if modules == "from-build":
        notes.append(f"Modules come from {q_bin}/build-info.json: {', '.join(env['build_info'].get('modules', [])) or 'none'}.")
    elif args.modules is None and env["has_modules"]:
        notes.append(
            "No modules were written: it is not known which ones the Q binaries were compiled with "
            f"(loaded in this shell: {', '.join(env['loaded_modules']) or 'none'}). Load the build modules "
            "and run `qligfep_config build-info`, then set `modules: from-build`; or pass --modules."
        )
    if not local and account is None and env["accounts"]:
        notes.append(f"Several Slurm accounts are available ({', '.join(env['accounts'])}); set `slurm.account`.")

    target = None if args.print_only else (args.output or str(user_config_path()))
    if target is None:
        print(text)
    else:
        from pathlib import Path

        path = Path(target).expanduser()
        if path.exists():
            print(text)
            print(
                f"\n{path} already exists and was left untouched. The proposal for this machine is printed "
                "above; copy what you need into it.",
                file=sys.stderr,
            )
            return 1
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        print(f"Wrote {path}")
        notes.append("Next: `qligfep_config check` tests the profile on this machine.")
    for note in notes:
        print(f"Note: {note}", file=sys.stderr)
    return 0


def cmd_show(args: argparse.Namespace) -> int:
    config = load_config(project_file=args.config)
    if args.profile or args.protocol:
        data = {}
        if args.profile:
            data["profile"] = config.get_profile(args.profile, system=args.system).as_dict()
        if args.protocol:
            data["protocol"] = config.get_protocol(args.protocol).as_dict()
        print(yaml.safe_dump(data, sort_keys=False), end="")
        return 0

    print("Configuration files, lowest to highest priority:")
    for source in config.sources:
        print(f"  {source}")
    user_file = user_config_path()
    if str(user_file) not in config.sources:
        print(f"  (no user file at {user_file}; create one with `qligfep_config init`)")
    if len(config.sources) == 1 + int(str(user_file) in config.sources):
        print("  (no qligfep.yaml in this directory)")

    print(f"\nProfiles (default: {config.default_profile or 'none set'}):")
    for name in config.profile_names():
        try:
            profile = config.get_profile(name)
        except ConfigError as error:
            print(f"  {name:<14} INVALID: {error}")
            continue
        cores = f"{profile.resources.ntasks} core" + ("s" if profile.resources.ntasks != 1 else "")
        origin = profile.sources[-1] if profile.sources else ""
        print(f"  {name:<14} {profile.execution:<12} {profile.engine:<6} {cores:<9} {origin}")

    print(f"\nProtocols (default: {config.default_protocol}):")
    for name in config.protocol_names():
        try:
            protocol = config.get_protocol(name)
        except ConfigError as error:
            print(f"  {name:<24} INVALID: {error}")
            continue
        print(f"  {name:<24} {protocol.method:<15} {protocol.description}")
    print("\nPrint one in full with `qligfep_config show --profile NAME` or `--protocol NAME`.")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    config = load_config(project_file=args.config)
    problems = 0
    for name in config.profile_names():
        try:
            config.get_profile(name, q_bin=args.q_bin)
        except ConfigError as error:
            problems += 1
            print(f"[fail] profile {name}: {error}")
    for name in config.protocol_names():
        try:
            config.get_protocol(name)
        except ConfigError as error:
            problems += 1
            print(f"[fail] protocol {name}: {error}")
    if not problems:
        print(f"[ok]   {len(config.profile_names())} profiles and {len(config.protocol_names())} protocols are valid")

    if not (args.profile or config.default_profile):
        print("No profile selected and no default_profile set: skipping the machine checks (use --profile NAME).")
        return 1 if problems else 0
    profile = config.get_profile(args.profile, q_bin=args.q_bin)
    print(f"\nProfile '{profile.name}' ({profile.execution}, {profile.engine}, binaries from {profile.q_bin_source}):")
    for status, check, detail in check_profile(profile):
        problems += status == "fail"
        print(f"[{status}]".ljust(7) + f"{check}: {detail}")
    if profile.is_local is False:
        print("Note: these checks ran on this node. Compute nodes can differ; a first short job is the real test.")
    return 1 if problems else 0


def cmd_build_info(args: argparse.Namespace) -> int:
    q_bin, _ = resolve_q_bin(None, "build-info", args.q_bin)
    if not q_bin.is_dir():
        print(f"{q_bin} does not exist; compile the Q binaries first.", file=sys.stderr)
        return 1
    path = write_build_info(q_bin, args.modules)
    print(f"Wrote {path}")
    print(path.read_text(), end="")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handler = {"init": cmd_init, "show": cmd_show, "check": cmd_check, "build-info": cmd_build_info}[args.command]
    try:
        return handler(args)
    except ConfigError as error:
        print(f"Configuration error: {error}", file=sys.stderr)
        return 2


def main_exe():
    sys.exit(main())


if __name__ == "__main__":
    main_exe()
