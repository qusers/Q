"""Module to generate (setup) all FEP files for the directory you're working on."""

import argparse
import json
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from ..logger import logger, setup_logger
from .parser_base import parse_arguments

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def ligpairs_from_json(json_file):
    with open(json_file) as infile:
        json_dict = json.load(infile)
    try:
        edges = json_dict["edges"]  # should be a list of dictionaries
    except KeyError as kerr:
        raise KeyError('Could not find "edges" in json file') from kerr
    ligpairs = [(e["from"], e["to"]) for e in edges]
    return ligpairs


def create_call(**kwargs):
    """Function to dynamically create a call to QligFEP.cli.main_exe() based on the kwargs."""
    template = (
        "qligfep -l1 '{lig1}' -l2 '{lig2}' -FF {FF} -s {system} -c {cluster} -R {replicates} "
        "-S {sampling} -r {sphereradius} -l {start} -w {windows} -T {temperature} -ts {timestep} "
        "-rest {rest} -drf {dr_force} -log {log}"
    )
    if "cysbond" in kwargs and kwargs["cysbond"] is not None:
        template += " -b {cysbond}"
    if "to_clean" in kwargs and kwargs["to_clean"] is not None:
        template += " -clean {to_clean}"
    if "random_state" in kwargs and kwargs["random_state"] is not None:
        template += " -rs {random_state}"
    if "water_thresh" in kwargs and kwargs["water_thresh"] != 1.4:
        template += " -wath {water_thresh}"
    return template.format(**kwargs)


def submit_command(command: str) -> None:
    """Function to submit a command using the subprocess module.

    Args:
        command: string with the command to be submitted.
    """
    try:
        result = subprocess.run(command, shell=True, text=True, capture_output=True, check=True)
        logging.info(f"Command executed successfully: {command}\nOutput: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {command}\nError: {e.stderr}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> None:
    # setup the logger with the desired log level
    setup_logger(level=args.log)

    cwd = Path.cwd()

    if args.protein_only:
        systems = ["protein"]
        sys_directories = [cwd / "2.protein"]
    elif args.water_only:
        systems = ["water"]
        sys_directories = [cwd / "1.water"]
    else:
        systems = ["protein", "water"]
        sys_directories = [cwd / "2.protein", cwd / "1.water"]

    # make sure that the default directories for running the FEP calculations are there
    for sys_dir in sys_directories:
        if not sys_dir.exists():
            sys_dir.mkdir()

    if args.json_map is None:  # Try to load a json file from cwd
        json_files = list(cwd.glob("*.json"))
        if len(json_files) == 1:
            args.json_map = json_files[0]
        else:
            raise FileNotFoundError("No QmapFEP json file found in the current directory")

    lig_pairs = ligpairs_from_json(args.json_map)
    for system, sys_dir in zip(systems, sys_directories):
        for pair in lig_pairs:
            lig1 = pair[0]
            lig2 = pair[1]

            temp_dir = cwd / f"FEP_{lig1}_{lig2}"
            to_clean = " ".join(args.to_clean) if args.to_clean is not None else None
            command = create_call(
                lig1=lig1,
                lig2=lig2,
                FF=args.FF,
                system=system,
                cluster=args.cluster,
                sphereradius=args.sphereradius,
                cysbond=(args.cysbond if system == "protein" else None),
                start=args.start,
                temperature=args.temperature,
                replicates=args.replicates,
                sampling=args.sampling,
                timestep=args.timestep,
                windows=args.windows,
                to_clean=to_clean,
                random_state=args.random_state,
                rest=args.restraint_method,
                dr_force=args.dr_force,
                water_thresh=args.water_thresh,
                log=args.log,
            )
            logger.info(f"Submitting the command:\n{command}")
            dst = sys_dir / f"FEP_{lig1}_{lig2}"
            os.system(command)
            shutil.move(temp_dir, dst)


def main_exe():
    args = parse_arguments(program="setupFEP")
    main(args)


if __name__ == "__main__":
    main_exe()
