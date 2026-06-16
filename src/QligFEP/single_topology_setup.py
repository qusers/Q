"""Set up single-topology single-Hamiltonian FEP edges (benchmark driver).

Mirrors setupFEP's directory layout (``2.protein/FEP_x``, ``1.water/FEP_x``) but
builds single-topology hybrids via :class:`SingleTopologyFEP` instead of the dual
two-ligand layout. The orchestration runs in the perturbation directory, which
must contain the per-ligand ``.lib/.prm/.pdb/.sdf`` plus ``protein.pdb`` and
``water.pdb`` (exactly what a setupFEP-ready directory already holds).

Each edge is validated before it is moved into place: the hybrid atom counts must
match the endpoint ligands, both states must be net-neutral, every bonded atom
must be defined, and qprep must resolve all junction boundary parameters and emit
a topology. Edges that fail (e.g. a mapping that breaks at a linker-element change)
are reported and skipped rather than silently producing a broken simulation.
"""
import json
import shutil
from pathlib import Path

from .logger import logger
from .single_topology_fep import SingleTopologyFEP


class BuildValidationError(RuntimeError):
    """Raised when a single-topology build fails its sanity checks."""


def _lib_atom_count(lib_path):
    """Number of atoms in a Q .lib [atoms] block."""
    n, block = 0, None
    for line in Path(lib_path).read_text().splitlines():
        s = line.split()
        if not s:
            continue
        if s[0].startswith("["):
            block = s[0]
            continue
        if block == "[atoms]" and s[0].isdigit():
            n += 1
    return n


def _validate_mapping(run):
    """Sanity-check the hybrid atom table built by read_files."""
    roles = [h["role"] for h in run.hyb]
    n_real_a = sum(r in ("shared", "vanish") for r in roles)
    n_real_b = sum(r in ("shared", "appear") for r in roles)
    n_a, n_b = _lib_atom_count(run.lig1 + ".lib"), _lib_atom_count(run.lig2 + ".lib")
    q_a = sum(h["sA"][1] for h in run.hyb)
    q_b = sum(h["sB"][1] for h in run.hyb)

    problems = []
    if n_real_a != n_a:
        problems.append(f"real-in-A {n_real_a} != lig1 atoms {n_a}")
    if n_real_b != n_b:
        problems.append(f"real-in-B {n_real_b} != lig2 atoms {n_b}")
    if abs(q_a) > 0.05:
        problems.append(f"state-A net charge {q_a:+.4f}")
    if abs(q_b) > 0.05:
        problems.append(f"state-B net charge {q_b:+.4f}")
    if not run.same_charge:
        problems.append(
            f"charge change {run.charge_lig1}->{run.charge_lig2}; single-H estimator unproven on charged edges"
        )
    defined = {h["name"] for h in run.hyb}
    undefined = {atom for term in run.bonds for atom in term} - defined
    if undefined:
        problems.append(f"bond atoms undefined: {sorted(undefined)}")

    return {
        "n_hybrid": len(roles),
        "shared": roles.count("shared"),
        "vanish": roles.count("vanish"),
        "appear": roles.count("appear"),
        "real_a": n_real_a,
        "real_b": n_real_b,
        "q_a": round(q_a, 4),
        "q_b": round(q_b, 4),
        "mapping_problems": problems,
    }


def _validate_qprep(run, inputdir):
    """Confirm qprep resolved every junction parameter and built the topology."""
    out = Path(inputdir) / "qprep.out"
    missing = run._parse_missing_params(out) if out.exists() else [("file", ("qprep.out missing",))]
    top_built = (Path(inputdir) / "singletop.top").exists()
    problems = []
    if missing:
        problems.append(f"{len(missing)} unresolved parameter(s), e.g. {missing[0]}")
    if not top_built:
        problems.append("singletop.top not built")
    return {"qprep_missing": len(missing), "top_built": top_built, "qprep_problems": problems}


def setup_single_topology_edge(
    lig1,
    lig2,
    system,
    *,
    FF="AMBER14sb",
    cluster="SNELLIUS",
    sphereradius="25",
    windows=21,
    replicates="5",
    temperature="298",
    timestep="2fs",
    cysbond="auto",
    to_clean=("dcd",),
    random_state=42,
    softcore_method="gapsys",
    strict=True,
):
    """Build one single-topology edge for one leg (protein/water) in the cwd.

    Returns a diagnostics dict (atom partition, net charges, qprep status, writedir).
    With ``strict`` (default), raises :class:`BuildValidationError` on any mapping or
    qprep problem so a broken edge never reaches MD.
    """
    run = SingleTopologyFEP(
        lig1=lig1,
        lig2=lig2,
        FF=FF,
        system=system,
        cluster=cluster,
        sphereradius=str(sphereradius),
        start="0.0",
        temperature=str(temperature),
        replicates=str(replicates),
        timestep=timestep,
        cysbond=(cysbond if system == "protein" else "none"),
        to_clean=list(to_clean),
        random_state=random_state,
    )
    writedir = run.makedir()
    inputdir = writedir + "/inputfiles"

    a = run.read_files()
    diag = _validate_mapping(run)
    if diag["mapping_problems"] and strict:
        raise BuildValidationError(f"{lig1}->{lig2} ({system}) mapping: {'; '.join(diag['mapping_problems'])}")

    run.change_lib(a[0][1], inputdir)
    fep_vdw = run.change_prm(a[0][1], inputdir)
    run.merge_pdbs(inputdir)
    if system == "protein":
        run.write_water_pdb(inputdir)
        run.avoid_water_protein_clashes(inputdir, header=f"{run.sphereradius}.0 SPHERE")

    run.write_qprep(inputdir)
    run.qprep(inputdir)
    diag.update(_validate_qprep(run, inputdir))
    if diag["qprep_problems"] and strict:
        raise BuildValidationError(f"{lig1}->{lig2} ({system}) qprep: {'; '.join(diag['qprep_problems'])}")

    run.write_FEP_file(
        a[1][0], a[1][1], fep_vdw, inputdir, a[2][0], a[2][1],
        softcore_method=softcore_method, single_hamiltonian=True,
    )
    overlapping = run.set_restraints(writedir, "overlap")
    lambdas = [str(i) for i in range(int(windows))]
    file_list = run.write_md_files(lambdas, inputdir, a[2][0], a[2][1], overlapping)
    run.write_runfile(inputdir, file_list)
    run.write_submitfile(writedir)

    diag["writedir"] = writedir
    return diag


def run_benchmark(json_map="mapping.json", systems=("protein", "water"), **opts):
    """Build every edge in ``json_map`` for the requested legs, mirroring setupFEP's
    ``2.protein``/``1.water`` layout. Writes a build summary and continues past edges
    that fail validation (recording them), so one bad mapping never aborts the set."""
    cwd = Path.cwd()
    sys_dir = {"protein": cwd / "2.protein", "water": cwd / "1.water"}
    edges = json.loads(Path(json_map).read_text())["edges"]

    summary = []
    for system in systems:
        sys_dir[system].mkdir(exist_ok=True)
        for edge in edges:
            lig1, lig2 = edge["from"], edge["to"]
            dst = sys_dir[system] / f"FEP_{lig1}_{lig2}"
            if dst.exists():
                logger.info(f"skip existing {dst}")
                continue
            try:
                diag = setup_single_topology_edge(lig1, lig2, system, **opts)
            except Exception as ex:  # noqa: BLE001 - record and continue
                logger.error(f"BUILD FAILED {lig1}->{lig2} ({system}): {ex}")
                summary.append({"edge": f"{lig1}_{lig2}", "system": system, "ok": False, "error": str(ex)})
                partial = cwd / f"FEP_{lig1}_{lig2}"
                if partial.exists():
                    shutil.rmtree(partial)
                continue
            shutil.move(diag["writedir"], dst)
            diag.update({"edge": f"{lig1}_{lig2}", "system": system, "ok": True})
            summary.append(diag)
            logger.info(
                f"built {lig1}->{lig2} ({system}): {diag['n_hybrid']} atoms "
                f"(shared {diag['shared']}, vanish {diag['vanish']}, appear {diag['appear']}), "
                f"qprep_missing={diag['qprep_missing']}"
            )

    (cwd / "single_topology_build_summary.json").write_text(json.dumps(summary, indent=2))
    n_ok = sum(s["ok"] for s in summary)
    logger.info(f"single-topology build: {n_ok}/{len(summary)} edge-legs succeeded")
    return summary
