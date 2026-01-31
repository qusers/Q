"""qfep.inp template for FEP analysis input files.

Generates qfep input files that specify parameters for free energy
perturbation analysis calculations.
"""


def calculate_kT(temperature: float) -> str:
    """Calculate kT value at given temperature.

    Args:
        temperature: Temperature in Kelvin

    Returns:
        kT value as formatted string with 3 decimal places
    """
    k = 0.0019872041  # kcal/(mol.K)
    return f"{k * temperature:.3f}"


def format_energy_files(lambdas: list[str]) -> list[str]:
    """Generate energy filenames from lambda values.

    Converts lambda values to corresponding energy file names using the
    pattern md_{lambda1}_{lambda2}.en where lambda values run forward
    and backward through the list.

    Args:
        lambdas: List of lambda values as strings (e.g., ["1.000", "0.500", "0.000"])

    Returns:
        List of energy filenames (e.g., ["md_1000_0000.en", "md_0500_0500.en", ...])
    """
    filenames = []
    total = len(lambdas)
    for i in range(total):
        j = -(i + 1)  # Reverse index
        lambda1 = lambdas[i].replace(".", "")
        lambda2 = lambdas[j].replace(".", "")
        filenames.append(f"md_{lambda1}_{lambda2}.en")
    return filenames


def render_qfep_input(
    total_lambdas: int,
    temperature: float,
    windows: int,
    energy_files: list[str],
) -> str:
    """Render qfep.inp content.

    Args:
        total_lambdas: Total number of lambda values
        temperature: Temperature in Kelvin for kT calculation
        windows: Number of windows for FEP analysis
        energy_files: List of energy file names to include

    Returns:
        Complete qfep.inp file content as string
    """
    kT_value = calculate_kT(temperature)

    lines = [
        str(total_lambdas),
        "2  0",
        f"{kT_value}  {windows}",
        str(windows),
        str(windows),
        str(windows),
        "0",
        "0",
        "1 0",
    ]

    # Append energy files
    lines.extend(energy_files)

    return "\n".join(lines)
