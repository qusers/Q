"""The log of a run: one line per event, with the time, on standard error and in a file.

    from fep.log import log
    log.info("heating at lambda %.2f", lam)

setup() is called once by the program; without it the lines are not shown.
"""
import logging
import sys

log = logging.getLogger("fep")


def setup(path=None):
    """Show the log on standard error and, with `path`, append it to that file."""
    form = logging.Formatter("%(asctime)s  %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    handlers = [logging.StreamHandler(sys.stderr)]
    if path is not None:
        handlers.append(logging.FileHandler(path))
    log.handlers.clear()
    for handler in handlers:
        handler.setFormatter(form)
        log.addHandler(handler)
    log.setLevel(logging.INFO)
    log.propagate = False
