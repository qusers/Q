from ..logger import logger
from ..pdb_utils import write_dataframe_to_pdb


def atom_element_mapping():
    """
    Create a mapping of POPC atom names to their correct elements.

    Returns:
        dict: {atom_name: element_symbol}
    """
    # Based on AMBER14sb_popc-only.lib atom types and molecular weights
    element_mapping = {
        # Carbon atoms (molecular weight ~12.01)
        "C1": "C",
        "C2": "C",
        "C3": "C",
        "C11": "C",
        "C12": "C",
        "C13": "C",
        "C14": "C",
        "C15": "C",
        "C16": "C",
        "C17": "C",
        "C18": "C",
        "C19": "C",
        "C21": "C",
        "C22": "C",
        "C23": "C",
        "C24": "C",
        "C25": "C",
        "C26": "C",
        "C27": "C",
        "C28": "C",
        "C29": "C",
        "C31": "C",
        "C32": "C",
        "C33": "C",
        "C34": "C",
        "C35": "C",
        "C110": "C",
        "C111": "C",
        "C112": "C",
        "C113": "C",
        "C114": "C",
        "C115": "C",
        "C116": "C",
        "C210": "C",
        "C211": "C",
        "C212": "C",
        "C213": "C",
        "C214": "C",
        "C215": "C",
        "C216": "C",
        "C217": "C",
        "C218": "C",
        "CA1": "C",
        "CA2": "C",
        # Nitrogen atoms (molecular weight ~14.01)
        "N4": "N",
        "N31": "N",
        # Oxygen atoms (molecular weight ~16.00)
        "O7": "O",
        "O9": "O",
        "O10": "O",
        "O11": "O",
        "O12": "O",
        "O14": "O",
        "O16": "O",
        "O21": "O",
        "O22": "O",
        "O31": "O",
        "O32": "O",
        "O33": "O",
        "O34": "O",
        "O35": "O",
        # Phosphorus atoms (molecular weight ~30.97)
        "P8": "P",
        "P31": "P",
        # Hydrogen atoms (molecular weight ~1.01) - all H atoms start with H
        "HR": "H",
        "HS": "H",
        "HX": "H",
        "HA": "H",
        "HB": "H",
        "H1A": "H",
        "H1B": "H",
        "H2A": "H",
        "H2B": "H",
        "H2R": "H",
        "H2S": "H",
        "H3A": "H",
        "H3B": "H",
        "H3C": "H",
        "H3R": "H",
        "H3S": "H",
        "H4A": "H",
        "H4B": "H",
        "H4C": "H",
        "H4R": "H",
        "H4S": "H",
        "H5A": "H",
        "H5B": "H",
        "H5C": "H",
        "H5R": "H",
        "H5S": "H",
        "H6R": "H",
        "H6S": "H",
        "H7R": "H",
        "H7S": "H",
        "H8R": "H",
        "H8S": "H",
        "H9R": "H",
        "H9S": "H",
        "H10R": "H",
        "H10S": "H",
        "H11R": "H",
        "H11S": "H",
        "H12R": "H",
        "H12S": "H",
        "H13R": "H",
        "H13S": "H",
        "H14R": "H",
        "H14S": "H",
        "H15R": "H",
        "H15S": "H",
        "H16R": "H",
        "H16S": "H",
        "H16T": "H",
        "H22R": "H",
        "H22S": "H",
        "H23R": "H",
        "H23S": "H",
        "H24R": "H",
        "H24S": "H",
        "H25R": "H",
        "H25S": "H",
        "H26R": "H",
        "H26S": "H",
        "H27R": "H",
        "H27S": "H",
        "H28R": "H",
        "H28S": "H",
        "H29R": "H",
        "H20R": "H",
        "H21R": "H",
        "H21S": "H",
        "H22A": "H",
        "H22B": "H",
        "H23A": "H",
        "H23B": "H",
        "H24A": "H",
        "H24B": "H",
        "H25A": "H",
        "H25B": "H",
        "H26A": "H",
        "H26B": "H",
        "H26C": "H",
        "H27A": "H",
        "H27B": "H",
        "H28A": "H",
        "H28B": "H",
        "H28C": "H",
    }

    return element_mapping


def pymemdyn_to_amber_mapping():
    """
    Direct mapping from pymemdyn POPC atom names to unified format.
    This skips the buildH step and only handles heavy atoms.
    """
    return {
        # Choline head group - direct mapping from pymemdyn to unified
        "C1": "C34",  # Choline methyl carbon
        "C2": "C35",  # Choline methyl carbon
        "C3": "C33",  # Choline methyl carbon
        "N4": "N31",  # Choline nitrogen
        "C5": "C31",  # Choline CH2 carbon
        "C6": "C32",  # Choline CH2 carbon (connected to phosphate)
        # Phosphate group
        "O7": "O32",  # Phosphate-choline oxygen
        "P8": "P31",  # Phosphorus
        "O9": "O34",  # Phosphate oxygen
        "O10": "O33",  # Phosphate oxygen
        "O11": "O31",  # Phosphate-glycerol oxygen
        # Glycerol backbone
        "C12": "C3",  # Glycerol CH2 (sn-3 position)
        "C13": "C2",  # Glycerol CH (sn-2 position)
        "C32": "C1",  # Glycerol CH2 (sn-1 position)
        # sn-2 ester linkage and carbonyl
        "O14": "O21",  # sn-2 ester oxygen
        "C15": "C21",  # sn-2 ester carbonyl carbon
        "O16": "O22",  # sn-2 ester carbonyl oxygen
        # sn-1 ester linkage and carbonyl
        "O33": "O11",  # sn-1 ester oxygen
        "C34": "C11",  # sn-1 ester carbonyl carbon
        "O35": "O12",  # sn-1 ester carbonyl oxygen
        # sn-2 fatty acid chain (oleic acid chain)
        "C17": "C22",  # First carbon of sn-2 chain
        "C18": "C23",
        "C19": "C24",
        "C20": "C25",
        "C21": "C26",
        "C22": "C27",
        "C23": "C28",
        "C24": "C29",  # Before double bond
        "C25": "C210",  # Double bond carbon
        "C26": "C211",  # Double bond carbon
        "C27": "C212",  # After double bond
        "C28": "C213",
        "C29": "C214",
        "C30": "C215",
        "C31": "C216",
        "CA1": "C217",  # Additional carbons in oleic chain
        "CA2": "C218",  # Terminal methyl
        # sn-1 fatty acid chain (palmitic acid chain)
        "C36": "C12",  # First carbon of sn-1 chain
        "C37": "C13",
        "C38": "C14",
        "C39": "C15",
        "C40": "C16",
        "C41": "C17",
        "C42": "C18",
        "C43": "C19",
        "C44": "C110",
        "C45": "C111",
        "C46": "C112",
        "C47": "C113",
        "C48": "C114",
        "C49": "C115",
        "C50": "C116",  # Terminal methyl of palmitic chain
    }


def has_pop_residues(pdb_data):
    """
    Check if the PDB data contains any POP residues.

    Args:
        pdb_data: DataFrame containing PDB data

    Returns:
        tuple: (bool: True if POP residues are found, int: Number of POP residues found)
    """
    # Use value_counts for more efficient counting
    residue_counts = pdb_data["residue_name"].value_counts()
    has_pop = "POP" in residue_counts

    if has_pop:
        # Count unique residue instances (chain_id + residue_seq_number combinations)
        pop_residues = pdb_data[pdb_data["residue_name"] == "POP"]
        pop_count = len(pop_residues.groupby(["chain_id", "residue_seq_number"]))
        return True, pop_count
    else:
        return False, 0


def convert_pymemdyn_to_unified_dataframe(pdb_data):
    """Convert pymemdyn POPC format directly to unified format for DataFrame using vectorized operations."""
    logger.info("Converting POP residues from pymemdyn to unified format...")

    # Get mapping
    mapping = pymemdyn_to_amber_mapping()

    # Create a mask for POP residues
    pop_mask = pdb_data["residue_name"] == "POP"

    if not pop_mask.any():
        return pdb_data

    pop_residues_count = len(pdb_data[pop_mask].groupby(["chain_id", "residue_seq_number"]))
    pdb_data_corrected = pdb_data.copy()

    # Apply atom name mapping only to the POPC atoms
    pop_atoms = pdb_data_corrected.loc[pop_mask, "atom_name"]
    mapped_atoms = pop_atoms.map(mapping)

    # keep original for unmapped atoms
    unmapped_mask = mapped_atoms.isna()
    mapped_atoms.loc[unmapped_mask] = pop_atoms.loc[unmapped_mask]
    pdb_data_corrected.loc[pop_mask, "atom_name"] = mapped_atoms

    # Set default values for missing fields
    occupancy_mask = pop_mask & (
        pdb_data_corrected["occupancy"].isna() | (pdb_data_corrected["occupancy"] == "")
    )
    temp_factor_mask = pop_mask & (
        pdb_data_corrected["temp_factor"].isna() | (pdb_data_corrected["temp_factor"] == "")
    )
    charge_mask = pop_mask & (pdb_data_corrected["charge"].isna() | (pdb_data_corrected["charge"] == ""))

    pdb_data_corrected.loc[occupancy_mask, "occupancy"] = 1.00
    pdb_data_corrected.loc[temp_factor_mask, "temp_factor"] = 0.00
    pdb_data_corrected.loc[charge_mask, "charge"] = ""

    logger.debug(f"POP lipid conversion completed for {pop_residues_count} residue(s)")
    logger.debug("Element symbols will be corrected based on atom names during file writing")

    unmapped_atoms = set(pop_atoms.loc[unmapped_mask].unique()) if unmapped_mask.any() else set()
    if unmapped_atoms:
        logger.warning(f"Unmapped atoms (kept original names): {unmapped_atoms}")

    return pdb_data_corrected


def correct_popc_elements(df):
    """
    Correct element symbols for POP residues using vectorized operations.
    """
    element_mapping = atom_element_mapping()
    pop_mask = df["residue_name"] == "POP"

    if not pop_mask.any():
        return df

    df_corrected = df.copy()

    pop_atoms = df_corrected.loc[pop_mask, "atom_name"]
    corrected_elements = pop_atoms.map(element_mapping)

    # For atoms not in mapping, guess from first character (fallback)
    missing_mask = corrected_elements.isna()
    corrected_elements.loc[missing_mask] = pop_atoms.loc[missing_mask].str[0].fillna("C")

    df_corrected.loc[pop_mask, "element_symbol"] = corrected_elements
    return df_corrected


def df_to_pdb_corrected_element(df, output_file, header=None):
    """
    Write DataFrame to PDB with element correction for lipid atoms using vectorized operations.
    """
    has_pop, _ = has_pop_residues(df)
    if has_pop:
        df_corrected = correct_popc_elements(df)
        write_dataframe_to_pdb(df_corrected, output_file, header=header)
    else:
        write_dataframe_to_pdb(df, output_file, header=header)
