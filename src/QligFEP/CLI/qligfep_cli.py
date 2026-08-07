"""Module containing the QligFEP command line interface."""

import argparse
import datetime
import json
from pathlib import Path
from typing import Optional

from QligFEP import __version__

from ..logger import logger, setup_logger
from ..pdb_utils import read_pdb_to_dataframe, residue_atom_serial_range
from ..qligfep import COUNTER_WATER_RESNAME, QligFEP
from ..templates.sections import format_wall_restraints
from .parser_base import parse_arguments


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> None:
    """Main function for qligfep_cli.py. Takes arguments from argparse and passes them
    to QligFEP class. If no arguments are given, the function will use the keyword arguments
    that are passed to it.

    Args:
        args: argparse.Namespace object containing all the arguments.
        kwargs: keyword arguments that will be passed to QligFEP class.
    """
    if args is not None:
        param_dict = {
            "lig1": args.lig1,
            "lig2": args.lig2,
            "FF": args.FF,
            "system": args.system,
            "cluster": args.cluster,
            "sphereradius": args.sphereradius,
            "cysbond": args.cysbond,
            "start": args.start,
            "temperature": args.temperature,
            "replicates": args.replicates,
            "sampling": args.sampling,
            "timestep": args.timestep,
            "minimize": args.minimize,
            "to_clean": args.to_clean,
            "water_thresh": args.water_thresh,
            "dr_force": args.dr_force,
            "random_state": args.random_state,
            "wath_ligand_only": args.wath_ligand_only,
            "neq": args.neq,
            "neq_reps": args.neq_reps,
            "neq_steps": args.neq_steps,
            "neq_eq_steps": args.neq_eq_steps,
            "neq_relax_steps": args.neq_relax_steps,
            "neq_L": args.neq_L,
            "neq_schedule": args.neq_schedule,
            "softcore_method": args.softcore_method,
            "charge_method": args.charge_method,
        }
        if args.protein_charge is not None:
            param_dict["protein_charge"] = args.protein_charge
    else:
        param_dict = {}
    param_dict.update(kwargs)
    setup_logger(level=args.log.upper())
    softcore_method = param_dict.pop("softcore_method", "standard")
    run = QligFEP(**param_dict)
    param_dict["softcore_method"] = softcore_method

    writedir = run.makedir()
    inputdir = writedir + "/inputfiles"

    command_str = "qligfep"
    for k, v in param_dict.items():
        if k == "FF":
            command_str += f" -{k} {v}"
        elif k == "to_clean":
            if v:
                command_str += f" --files-to-clean {' '.join(v)}"
        elif k == "water_thresh":
            command_str += f" --{k.replace('_', '-')} {v}"
        elif k == "wath_ligand_only":
            if v:
                command_str += f" --{k}".replace("_", "-")
        elif k == "dr_force":
            command_str += f" --{k} {v}".replace("dr_force", "distance_restraint_force")
        elif k in ("neq", "minimize"):
            if v:
                command_str += f" --{k}"
        elif k == "neq_L":
            command_str += f" --neq-steepness {v}"
        elif k in ("neq_reps", "neq_steps", "neq_eq_steps", "neq_relax_steps", "neq_schedule"):
            command_str += f" --{k.replace('_', '-')} {v}"
        elif k == "protein_charge":
            if v is not None:
                command_str += f" -pq {v}"
        elif k == "softcore_method":
            if v != "standard":
                command_str += f" --softcore-method {v}"
        elif k == "charge_method":
            if v != "ion_match":
                command_str += f" --charge-method {v}"
        else:
            command_str += f" --{k} {v}"
    command_str += f" --restraint_method {args.restraint_method}"
    time_now = datetime.datetime.now()
    date_str = time_now.strftime("%Y-%m-%d %H:%M:%S")
    (Path(inputdir) / "fep_config.json").write_text(
        json.dumps(
            {**param_dict, **{"time": date_str, "command_str": command_str, "QligFEP_version": __version__}},
            indent=4,
        )
    )

    a = run.read_files()
    changes_for_libfiles = a[0][1]
    changes_for_prmfiles = a[0][1]
    change_charges = a[1][0]
    change_vdw = a[1][1]
    lig_size1, lig_size2 = a[2][0], a[2][1]

    # Write the merged files
    logger.debug("Writing changes on files")
    run.change_lib(changes_for_libfiles, inputdir)
    logger.debug("Changing parameters")
    FEP_vdw = run.change_prm(changes_for_prmfiles, inputdir)
    logger.debug("Writing PDB files")
    run.merge_pdbs(inputdir)

    # Dispatch on --charge-method for the pre-qprep step.
    # 'ion_match' places real Cl-/Na+ ions in the water-leg PDB before qprep.
    # 'coalchemical_water' appends real waters as solute atoms before qprep
    #   so qprep includes them in the topology; their Q-atom indices are then
    #   refreshed from top_p.pdb after qprep renumbers atoms.
    # 'none' does nothing at setup time.
    charge_method = param_dict.get("charge_method", "ion_match")
    if charge_method == "ion_match":
        n_ions = run.place_counter_ions(inputdir)
        if n_ions > 0:
            logger.info(f"Placed {n_ions} {run.ion_type} counter-ion(s) for charge-changing perturbation")
    elif charge_method == "coalchemical_water":
        n_cw = run.place_counter_water(inputdir)
        if n_cw > 0:
            logger.info(f"Placed {n_cw} co-alchemical counter-water(s) ({run.ion_type} swap)")
    elif charge_method == "none" and run.same_charge is False:
        logger.warning("charge_method=none on a charge-changing edge: ddG will retain the Born artifact")

    run.write_water_pdb(inputdir)

    if not run.neq:
        logger.debug("Getting the lambdas")
        lambdas = run.get_lambdas(args.windows, args.sampling)
    logger.debug("Writing atom mapping for distance restraints")

    run.avoid_water_protein_clashes(
        inputdir, header=f"{run.sphereradius}.0 SPHERE", save_removed=(args.log in ["trace", "debug"])
    )

    logger.debug("Writing the QPREP files & running qprep")
    run.write_qprep(inputdir)
    run.qprep(inputdir)

    # Refresh the counter-water topology indices: qprep renumbers atoms.
    if charge_method == "coalchemical_water":
        run.update_counter_water_indices(inputdir)

    logger.debug("Writing FEP files")
    run.write_FEP_file(
        change_charges,
        change_vdw,
        FEP_vdw,
        inputdir,
        lig_size1,
        lig_size2,
        softcore_method=softcore_method,
    )
    overlapping_atoms = run.set_restraints(writedir, args.restraint_method, strict_check=True)

    if run.neq:
        logger.debug("Writing the non-equilibrium MD files")
        file_list = run.write_MD_neq(inputdir, lig_size1, lig_size2, overlapping_atoms)
        run.write_neq_runfile(inputdir, file_list)
        logger.debug(f"Generated files: {file_list}")
        logger.debug("Writing the submit files")
        run.write_submitfile(writedir)
        return
    # Build wall restraints for counter-ions or coalchemical waters (after qprep so top_p.pdb exists).
    wall_restraints_str = ""
    if run.n_counter_ions > 0:
        pdb_df = read_pdb_to_dataframe(Path(inputdir) / "top_p.pdb")
        atom_range = residue_atom_serial_range(pdb_df, [run.ion_type, COUNTER_WATER_RESNAME])
        if atom_range is not None:
            first_ion, last_ion = atom_range
            wall_radius = int(run.sphereradius) - 5
            wall_restraints_str = format_wall_restraints(first_ion, last_ion, wall_radius, force=1.0)
            logger.debug(f"Wall restraints for counter atoms: {first_ion}-{last_ion}, radius {wall_radius}")

    # Handling the correct offset here
    logger.debug("Writing the MD files")
    file_list = run.write_md_files(
        lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms, wall_restraints=wall_restraints_str
    )
    run.write_runfile(inputdir, file_list)
    logger.debug(f"Generated files: {file_list}")
    logger.debug("Writing the submit files")
    run.write_submitfile(writedir)
    logger.debug("Writing the QFEP files")
    run.write_qfep(args.windows, lambdas)


def main_exe():
    args = parse_arguments(program="QligFEP")
    main(args)


if __name__ == "__main__":
    main_exe()
