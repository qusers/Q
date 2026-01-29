"""FEP input file template system.

Provides Python-defined templates and configuration dataclasses for generating
MD input files (.inp) with explicit parameter control.

Examples:
    Inspect parameter dictionaries directly::

        >>> from QligFEP.templates import EQ5_2FS_PARAMS
        >>> print(EQ5_2FS_PARAMS)
        {'steps': 50000, 'stepsize': 2.0, 'temperature': 'T_VAR', ...}

    Render a complete input file::

        >>> from QligFEP.templates import get_equilibration_configs, render_md_input
        >>> configs = get_equilibration_configs("2fs", shell_radius=25)
        >>> eq5 = configs[4]
        >>> content = render_md_input(
        ...     params=eq5.params,
        ...     lambda1="0.500",
        ...     lambda2="0.500",
        ...     trajectory_file="eq5.dcd",
        ...     final_file="eq5.re",
        ...     restart_file="eq4.re",
        ... )
        >>> print(content)
        [MD]
        steps                     50000
        stepsize                  2.0
        ...

    Render a production MD file::

        >>> from QligFEP.templates import get_production_config, render_md_input
        >>> config = get_production_config("2fs", shell_radius=25)
        >>> content = render_md_input(
        ...     params=config.params,
        ...     lambda1="0.500",
        ...     lambda2="0.500",
        ...     trajectory_file="md_0500_0500.dcd",
        ...     final_file="md_0500_0500.re",
        ...     restart_file="eq5.re",
        ...     energy_file="md_0500_0500.en",
        ... )
"""

from QligFEP.templates.equilibration import (
    EQ1_PARAMS,
    EQ2_1FS_PARAMS,
    EQ2_2FS_PARAMS,
    EQ3_1FS_PARAMS,
    EQ3_2FS_PARAMS,
    EQ4_1FS_PARAMS,
    EQ4_2FS_PARAMS,
    EQ5_1FS_PARAMS,
    EQ5_2FS_PARAMS,
    EquilibrationConfig,
    get_equilibration_configs,
)
from QligFEP.templates.md_template import MDParameters, render_md_input
from QligFEP.templates.production import (
    PRODUCTION_1FS_PARAMS,
    PRODUCTION_2FS_PARAMS,
    ProductionConfig,
    get_production_config,
)
from QligFEP.templates.sections import (
    format_distance_restraints,
    format_sequence_restraint,
)

__all__ = [
    # Core
    "MDParameters",
    "render_md_input",
    # Equilibration
    "EquilibrationConfig",
    "get_equilibration_configs",
    "EQ1_PARAMS",
    "EQ2_2FS_PARAMS",
    "EQ3_2FS_PARAMS",
    "EQ4_2FS_PARAMS",
    "EQ5_2FS_PARAMS",
    "EQ2_1FS_PARAMS",
    "EQ3_1FS_PARAMS",
    "EQ4_1FS_PARAMS",
    "EQ5_1FS_PARAMS",
    # Production
    "ProductionConfig",
    "get_production_config",
    "PRODUCTION_2FS_PARAMS",
    "PRODUCTION_1FS_PARAMS",
    # Sections
    "format_distance_restraints",
    "format_sequence_restraint",
]
