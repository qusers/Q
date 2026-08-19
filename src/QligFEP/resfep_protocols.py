"""Named protocol settings for residue FEP."""

DEFAULT_PRODUCTION_STEPS = 10_000

MANUSCRIPT_SEEDS = (
    971,
    2856,
    8253,
    1405,
    5342,
    7486,
    5374,
    1123,
    8397,
    8884,
)

MANUSCRIPT_SETTINGS = {
    "force_field": "OPLSAAM",
    "tripeptide_flanks": "A",
    "windows": 50,
    "sampling": "exponential",
    "start": "1",
    "timestep": "2fs",
    "temperature": "298",
    "replicates": len(MANUSCRIPT_SEEDS),
    "eq5_steps": 2_500_000,
    "production_steps": DEFAULT_PRODUCTION_STEPS,
    "separate_scaling": "off",
    "random_state": None,
    "seeds": MANUSCRIPT_SEEDS,
}

MANUSCRIPT_PREPARATION_SETTINGS = {
    "radius": 25.0,
    "center": None,
    "neutralization_offset": 3.0,
    "legs": ("protein", "tripeptide"),
}


def apply_manuscript_settings(namespace: object, *, include_preparation: bool = False) -> None:
    """Apply the manuscript protocol to a parsed CLI namespace."""
    settings = MANUSCRIPT_SETTINGS
    if include_preparation:
        settings = {**settings, **MANUSCRIPT_PREPARATION_SETTINGS}

    for name, value in settings.items():
        setattr(namespace, name, list(value) if isinstance(value, tuple) else value)
