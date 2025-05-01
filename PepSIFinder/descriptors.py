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