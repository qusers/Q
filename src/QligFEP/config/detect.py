"""Inspect the machine: scheduler, modules, toolchain, and whether a profile can run.

Module names differ too much between sites to guess them reliably. Detection therefore
prefers facts over guesses: the modules that were loaded when the binaries were built
(recorded in ``build-info.json`` next to them) and the modules loaded right now.
"""

from __future__ import annotations

import base64
import datetime
import json
import os
import platform
import re
import shutil
import subprocess
from multiprocessing import cpu_count
from pathlib import Path

from .schema import BUILD_INFO_NAME, Profile, read_build_info


def _run(command: list[str], timeout: int = 20) -> str | None:
    """Return a command's stdout, or None if it is missing, fails or hangs."""
    if shutil.which(command[0]) is None:
        return None
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def loaded_modules() -> list[str]:
    """Environment modules loaded in this shell (``$LOADEDMODULES``), in load order."""
    return [module for module in os.environ.get("LOADEDMODULES", "").split(":") if module]


def requested_modules() -> list[str]:
    """Modules the user loaded, without the dependencies they pulled in.

    Loading one toolchain module typically loads a dozen more. Lmod records which ones
    were asked for (``stackDepth`` 0) in its module table; replaying only those in a job
    reproduces the environment without pinning every dependency. Falls back to every
    loaded module for other module systems.
    """
    try:
        size = int(os.environ.get("_ModuleTable_Sz_", "0"))  # noqa: SIM112 - Lmod's spelling
        encoded = "".join(os.environ.get(f"_ModuleTable{index:03d}_", "") for index in range(1, size + 1))
        table = base64.b64decode(encoded).decode() if encoded else ""
    except ValueError:
        table = ""
    entries = re.findall(
        r'fullName\s*=\s*"([^"]+)".*?loadOrder\s*=\s*(\d+).*?stackDepth\s*=\s*(\d+)', table, flags=re.DOTALL
    )
    if not entries:
        return loaded_modules()
    top_level = sorted((int(order), name) for name, order, depth in entries if int(depth) == 0)
    return [name for _, name in top_level]


def _first_line(command: list[str]) -> str | None:
    output = _run(command)
    return output.strip().splitlines()[0] if output and output.strip() else None


def slurm_accounts() -> list[str]:
    """Slurm accounts the current user can charge, if the accounting database answers."""
    user = os.environ.get("USER", "")
    output = _run(["sacctmgr", "-nP", "show", "associations", f"user={user}", "format=account"])
    return sorted({line.strip() for line in (output or "").splitlines() if line.strip()})


def slurm_partitions() -> list[str]:
    """Partition names from ``sinfo``; the cluster default is listed first."""
    output = _run(["sinfo", "-h", "-o", "%P"])
    names = [line.strip() for line in (output or "").splitlines() if line.strip()]
    default = [name.rstrip("*") for name in names if name.endswith("*")]
    return default + [name for name in names if not name.endswith("*")]


def detect_environment(q_bin: Path) -> dict:
    """Collect what ``qligfep_config init`` needs to propose a profile."""
    has_slurm = shutil.which("sbatch") is not None
    build_info = read_build_info(q_bin)
    return {
        "hostname": platform.node(),
        "system": platform.system(),
        "cpu_count": cpu_count(),
        "has_slurm": has_slurm,
        "has_modules": bool(os.environ.get("MODULESHOME") or os.environ.get("LMOD_CMD")),
        "loaded_modules": requested_modules(),
        "build_info": build_info,
        "mpirun": shutil.which("mpirun"),
        "accounts": slurm_accounts() if has_slurm else [],
        "partitions": slurm_partitions() if has_slurm else [],
        "binaries": {name: (Path(q_bin) / name).is_file() for name in ("qprep", "qdyn", "qdynp", "qfep")},
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
    }


def write_build_info(q_bin: Path, modules: list[str] | None = None) -> Path:
    """Record the toolchain active in this shell next to the Q binaries.

    Run it in the environment the binaries were compiled in. Profiles with
    ``modules: from-build`` then load exactly these modules in their jobs, which keeps
    several clones compiled with different toolchains apart without separate profiles.
    """
    info = {
        "written": datetime.datetime.now().isoformat(timespec="seconds"),
        "hostname": platform.node(),
        "modules": list(modules) if modules is not None else requested_modules(),
        "gfortran": _first_line(["gfortran", "--version"]),
        "mpif90": _first_line(["mpif90", "--version"]),
        "mpirun": _first_line(["mpirun", "--version"]),
        "conda_prefix": os.environ.get("CONDA_PREFIX"),
    }
    path = Path(q_bin) / BUILD_INFO_NAME
    path.write_text(json.dumps(info, indent=2) + "\n")
    return path


def _yaml_list(items: list[str]) -> str:
    return "[" + ", ".join(json.dumps(str(item)) for item in items) + "]"


def render_user_config(
    env: dict,
    name: str,
    execution: str,
    ntasks: int | str,
    time: str,
    modules: list[str] | str,
    account: str | None,
    partition: str | None,
) -> str:
    """Render a commented user configuration file for this machine.

    The file is written once and never rewritten by QligFEP, so the comments (and any
    the user adds) are preserved.
    """
    local = execution == "local"
    engine = "qdynp" if (env["binaries"].get("qdynp") and (env["mpirun"] or not local)) else "qdyn"
    if local and engine == "qdyn":
        ntasks = 1
    modules_value = modules if isinstance(modules, str) else _yaml_list(modules)

    facts = [
        f"host {env['hostname']}",
        "Slurm found" if env["has_slurm"] else "no Slurm (sbatch not on PATH)",
        f"{env['cpu_count']} cores",
        f"mpirun: {env['mpirun'] or 'not found'}",
    ]
    lines = [
        "# QligFEP user configuration: execution profiles for this machine.",
        f"# Generated by `qligfep_config init` ({'; '.join(facts)}).",
        "# QligFEP never rewrites this file, so edit it freely. Check it with",
        "# `qligfep_config check`, and see the result with `qligfep_config show`.",
        "#",
        "# Protocols (what is simulated) do not belong here: put them in a qligfep.yaml",
        "# inside the project directory, where they travel with the results.",
        "version: 1",
        "",
        "# Profile used when -c/--profile is not given on the command line.",
        f"default_profile: {name}",
        "",
        "profiles:",
        f"  {name}:",
        "    # slurm-array: each replicate is a task of a Slurm job array",
        "    # local:       replicates run on this machine, started by FEP_submit.sh",
        f"    execution: {execution}",
        "    # qdynp (MPI, needed for non-equilibrium FEP) | qdyn (serial)",
        f"    engine: {engine}",
        "    # Directory with qprep/qdyn/qdynp/qfep. Leave unset to use the binaries compiled",
        "    # inside the QligFEP installation that runs the command, so every clone uses its",
        "    # own engine. Set it (or $QLIGFEP_Q_BIN, or --q-bin) to use another clone's.",
        "    # q_bin: /path/to/Q/src/q6/bin/q6",
        "    # Modules loaded in the job, in order. `from-build` loads the ones recorded next",
        "    # to the binaries by `qligfep_config build-info`.",
        f"    modules: {modules_value}",
        "    # Shell lines run in the job before the modules are loaded, e.g. `module purge`.",
        "    # Jobs inherit the submit environment (slurm.export: ALL), so activating your",
        "    # environment before submitting is normally enough and this can stay empty.",
        '    setup: ""',
        "    # Override how the MPI engine is launched (default: mpirun with core binding",
        "    # under Slurm). `$ncores` is available, e.g. \"srun -n $ncores\".",
        "    # mpi_command: null",
        "    resources:",
        f"      ntasks: {ntasks}" + ("            # cores per replicate; auto = all cores" if local else "              # cores per replicate"),
        f'      time: "{time}"' + ("" if not local else "   # ignored for local execution"),
        "      # MB per core. auto = 512 (equilibrium) / 1700 (non-equilibrium); null omits it.",
        "      mem_per_cpu: auto",
        "      # Local execution only: replicates run at the same time.",
        "      max_parallel_replicates: 1",
        "      # Non-equilibrium jobs fire 2 x switch_pairs switches concurrently.",
        "      # auto = enough cores for all of them at once, never fewer than ntasks.",
        "      # A number is used as is; fewer cores than switches run them in waves.",
        "      neq_cores: auto",
    ]
    if not local:
        accounts = ", ".join(env["accounts"]) or "none reported by sacctmgr"
        partitions = ", ".join(env["partitions"]) or "none reported by sinfo"
        lines += [
            "    slurm:",
            "      # null submits without -A, so Slurm charges your default account.",
            f"      # Accounts you can charge: {accounts}",
            f"      account: {account if account else 'null'}",
            f"      # Partitions: {partitions}",
            f"      partition: {partition if partition else 'null'}",
            "      # Pass the submit environment to the job. Set to null to leave it to the site.",
            "      export: ALL",
            "      # Any other #SBATCH options, e.g. [\"--qos=short\", \"--constraint=zen3\"].",
            "      extra_directives: []",
        ]
    lines += [
        "    # Different settings per leg, e.g. fewer cores for the small water system:",
        "    # systems:",
        "    #   water: {resources: {ntasks: 4}}",
        "",
        "  # A second profile for another engine build on the same machine:",
        f"  # {name}-dev:",
        f"  #   extends: {name}",
        "  #   q_bin: /path/to/other/clone/src/q6/bin/q6",
        "  #   modules: from-build",
        "",
    ]
    return "\n".join(lines)


def check_profile(profile: Profile) -> list[tuple[str, str, str]]:
    """Check that a profile can run here. Returns ``(status, check, detail)`` rows.

    ``status`` is ``ok``, ``warn`` or ``fail``. The environment checks run in a shell that
    executes the profile's setup lines and module loads first, the way a job would.
    """
    rows: list[tuple[str, str, str]] = []
    needed = ["qprep", "qfep", profile.engine] + (["qdyn"] if profile.use_mpi else [])
    for name in dict.fromkeys(needed):
        path = profile.binary(name)
        if path.is_file() and os.access(path, os.X_OK):
            rows.append(("ok", f"binary {name}", str(path)))
        else:
            rows.append(("fail", f"binary {name}", f"{path} is missing or not executable"))

    if not profile.is_local:
        found = shutil.which("sbatch")
        rows.append(("ok", "sbatch", found) if found else ("fail", "sbatch", "not on PATH; cannot submit from here"))

    engine = profile.binary(profile.engine)
    linker = "otool -L" if platform.system() == "Darwin" else "ldd"
    script = (
        profile.environment_block()
        + "echo \"MPIRUN=$(command -v mpirun)\"\n"
        + f"{linker} '{engine}' 2>&1 | grep -i 'not found' | sed 's/^/MISSINGLIB=/'\n"
        + "true\n"
    )
    try:
        result = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=120, check=False)
    except (OSError, subprocess.SubprocessError) as error:
        rows.append(("warn", "job environment", f"could not be probed: {error}"))
        return rows

    stderr = [line for line in result.stderr.splitlines() if line.strip()]
    if profile.modules or profile.setup.strip():
        if stderr:
            rows.append(("warn", "setup/modules", "; ".join(stderr[:3])))
        else:
            rows.append(("ok", "setup/modules", "ran without messages"))
    mpirun = next((line[7:] for line in result.stdout.splitlines() if line.startswith("MPIRUN=")), "")
    if profile.use_mpi and not profile.mpi_command:
        rows.append(("ok", "mpirun", mpirun) if mpirun else ("fail", "mpirun", "not found after setup/modules"))
    missing = [line[11:].strip() for line in result.stdout.splitlines() if line.startswith("MISSINGLIB=")]
    if missing and engine.is_file():
        rows.append(("fail", f"{profile.engine} libraries", "; ".join(missing[:3])))
    elif engine.is_file():
        rows.append(("ok", f"{profile.engine} libraries", "all shared libraries resolve"))
    return rows
