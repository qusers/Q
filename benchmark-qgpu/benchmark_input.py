from pathlib import Path
import shutil


def read_inp_settings(input_file):
    """Read benchmark-relevant settings from a native Q .inp file."""
    input_path = Path(input_file).expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Q input file not found: {input_path}")
    if input_path.suffix.lower() != ".inp":
        raise ValueError(f"Q input file must end with .inp: {input_path}")

    section = None
    settings = {}
    with open(input_path, encoding="utf-8") as inp_f:
        for raw_line in inp_f:
            line = raw_line.split("!", 1)[0].strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                section = line[1:-1].strip().lower().replace("_", "-")
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            key = fields[0].lower().replace("_", "-")
            if section == "md" and key in {"steps", "stepsize"}:
                settings[key] = fields[1]
            elif section == "files" and key == "topology":
                settings["topology"] = fields[1]

    missing = [key for key in ("steps", "stepsize", "topology") if key not in settings]
    if missing:
        raise RuntimeError(f"Missing {', '.join(missing)} in {input_path}")
    try:
        steps = int(settings["steps"])
        stepsize_fs = float(settings["stepsize"])
    except ValueError as exc:
        raise RuntimeError(f"Invalid steps or stepsize in {input_path}") from exc
    if steps < 1 or stepsize_fs <= 0:
        raise RuntimeError(f"steps and stepsize must be positive in {input_path}")

    topology_path = Path(settings["topology"]).expanduser()
    if not topology_path.is_absolute():
        topology_path = (input_path.parent / topology_path).resolve()
    return {
        "path": input_path,
        "steps": steps,
        "stepsize_fs": stepsize_fs,
        "topology": topology_path,
    }


def stage_inp_input(input_file, run_dir):
    """Copy an input file into an isolated output directory for one task."""
    source = read_inp_settings(input_file)["path"]
    staged_dir = Path(run_dir) / "input"
    staged_dir.mkdir(parents=True, exist_ok=True)
    staged_input = staged_dir / source.name
    shutil.copyfile(source, staged_input)
    return staged_input


def count_atoms_from_inp(input_file):
    topology_path = read_inp_settings(input_file)["topology"]
    if not topology_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_path}")
    with open(topology_path, encoding="utf-8") as top_f:
        for line in top_f:
            if "Total no. of atoms" in line:
                try:
                    return int(line.split()[0])
                except (IndexError, ValueError) as exc:
                    raise RuntimeError(
                        f"Invalid atom-count line in {topology_path}: {line.strip()}"
                    ) from exc
    raise RuntimeError(f"Could not find atom count in topology: {topology_path}")


def ns_per_day(steps, stepsize_fs, wall_seconds):
    if wall_seconds <= 0:
        return None
    return steps * stepsize_fs * 1e-6 * 86400 / wall_seconds
