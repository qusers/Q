from pathlib import Path

def get_version():
    if (Path(__file__).parent / "_version.py").exists():
        from ._version import __version__  # noqa F401
    else:
        __version__ = "2.1.0"
    return __version__

SRC = Path(__file__).parents[1]

__version__ = get_version()