"""Execution profiles and simulation protocols loaded from YAML.

See ``defaults.yaml`` in this package for the annotated schema, and
``qligfep_config --help`` for the command that creates and inspects configuration.
"""

from .loader import (
    PROJECT_FILE_NAMES,
    Config,
    find_project_config,
    load_config,
    user_config_path,
)
from .schema import (
    ALL_STAGES,
    BUILD_INFO_NAME,
    STAGES,
    ConfigError,
    Profile,
    Protocol,
    StageOverride,
    read_build_info,
)

__all__ = [
    "ALL_STAGES",
    "BUILD_INFO_NAME",
    "PROJECT_FILE_NAMES",
    "STAGES",
    "Config",
    "ConfigError",
    "Profile",
    "Protocol",
    "StageOverride",
    "find_project_config",
    "load_config",
    "read_build_info",
    "user_config_path",
]
