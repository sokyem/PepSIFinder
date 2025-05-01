def add_descriptor(candidates_df: pd.DataFrame, proteins_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add protein descriptions to the peptide candidates DataFrame using a separate protein information DataFrame.

    Parameters:
        candidates_df (pd.DataFrame): Candidate peptides indexed by ['protein', 'Peptide'].
        proteins_df (pd.DataFrame): DataFrame with protein 'Accession' and 'Description'.

    Returns:
        pd.DataFrame: Merged DataFrame with protein descriptions added.
    """
    df = candidates_df.reset_index().merge(
        proteins_df[["Accession", "Description"]],
        left_on="protein",
        right_on="Accession",
        how="left"
    )
    df = df.set_index(keys=["protein", "Peptide"])
    return df

def extract_modified_peptides(file_path, modifications=["+14.02"], peptide_column="Peptide", residue=None):
    """
    Extracts a sorted list of unique peptides that contain specified modifications.
    
    Parameters:
        file_path (str): Path to the CSV file.
        modifications (list): List of modification mass shift strings (e.g. ["+14", "-18"]).
        peptide_column (str): Name of the column containing peptide sequences.
        residue (str, optional): If provided, only peptides with modifications on this residue (e.g., "D") are returned.
        
    Returns:
        list: Sorted list of unique peptide sequences containing the specified modifications.
    """
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return []
    
    # Create a regex for the modification(s) provided.
    # The pattern accepts an optional decimal point and optional digits after it.
    mod_regex = "(" + "|".join(re.escape(mod) for mod in modifications) + ")"
    
    if residue:
        # Pattern matching a specific residue followed by the modification.
        # For example, if residue="D" and mod is "+14", it matches D(+14) or D(+14.02)
        pattern = re.compile(r"(" + re.escape(residue) + r")\(" + mod_regex + r"(?:\.\d*)?\)")
    else:
        # Pattern matching any residue with the specified modification.
        pattern = re.compile(r"([A-Z])\(" + mod_regex + r"(?:\.\d*)?\)")
    
    modified_peptides = []
    for peptide in df[peptide_column].dropna().astype(str):
        if pattern.search(peptide):
            modified_peptides.append(peptide)
    
    return sorted(set(modified_peptides))