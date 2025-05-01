def slice_rt_range(df: pd.DataFrame, rt_min: float = 0, rt_max: float = 90) -> pd.DataFrame:
    """
    Filter peptides based on a defined retention time (RT) range.

    Parameters:
        df (pd.DataFrame): DataFrame containing peptide identifications with an 'RT' column.
        rt_min (float): Minimum allowed RT value (inclusive). Defaults to 0.
        rt_max (float): Maximum allowed RT value (exclusive). Defaults to 50.

    Returns:
        pd.DataFrame: Filtered DataFrame containing peptides within the specified RT range.
    """
    return df[(df["RT"] >= rt_min) & (df["RT"] < rt_max)]

def compare_peptide_df(df1: pd.DataFrame, df2: pd.DataFrame):
    """
    Compare two peptide DataFrames and identify common and unique peptides.

    Parameters:
        df1 (pd.DataFrame): First peptide DataFrame.
        df2 (pd.DataFrame): Second peptide DataFrame.

    Returns:
        tuple: Three DataFrames representing:
            - Peptides found in both dataframes.
            - Peptides only in df1.
            - Peptides only in df2.
    """
    merged_df = df1.merge(df2, on='Peptide', how='outer', indicator=True)
    common_peptides = merged_df[merged_df['_merge'] == 'both'].drop(columns=['_merge'])
    unique_to_df1 = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge'])
    unique_to_df2 = merged_df[merged_df['_merge'] == 'right_only'].drop(columns=['_merge'])
    return common_peptides, unique_to_df1, unique_to_df2