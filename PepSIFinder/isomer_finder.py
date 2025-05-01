import pandas as pd
import numpy as np

def get_peptide_rt_range(df):
    """
    Compute the first and last retention time and their difference for each unique peptide.
    
    Parameters:
        df (DataFrame): Input DataFrame containing at least the 'Peptide' and 'RT' columns.
        
    Returns:
        DataFrame: A DataFrame with columns 'Peptide', 'first_RT', 'last_RT', and 'max_RT_diff'.
    """
    # Ensure RT is numeric
    df["RT"] = pd.to_numeric(df["RT"], errors="coerce")
    
    # Group by peptide to compute min and max RT values
    peptide_rt = df.groupby("Peptide")["RT"].agg(["min", "max"]).reset_index()
    peptide_rt.rename(columns={"min": "first_RT", "max": "last_RT"}, inplace=True)
    peptide_rt["max_RT_diff"] = peptide_rt["last_RT"] - peptide_rt["first_RT"]
    return peptide_rt

def get_ims(df):
    """
    Compute mean 1/k0 values from the '1/k0 Range' column.
    
    Parameters:
        df (DataFrame): DataFrame with a '1/k0 Range' column (e.g., "0.8-1.0").
    
    Returns:
        Series: Mean values computed from the 1/k0 range.
    """
    temp_df = df["1/k0 Range"].str.split("-", expand=True).apply(pd.to_numeric, errors="coerce")
    return temp_df.mean(axis=1)

def get_mob_frac_variation(df):
    """
    Compute fractional variation of ion mobility (1/k0) per peptide and charge.
    
    Parameters:
        df (DataFrame): Input DataFrame with columns 'Peptide', 'z', and '1/k0'.
    
    Returns:
        DataFrame: A DataFrame with columns 'Peptide', 'z', and 'frac_variation'.
    """
    frac_df = df.groupby(["Peptide", "z"])["1/k0"].agg(
        lambda x: (x.max() - x.min()) / x.mean()
    ).reset_index(name="frac_variation")
    return frac_df

def get_candidates(df, RT_tol=0.5, ims_diff=0.05):
    """
    Identify isomer candidates based on retention time and ion mobility variation.
    
    Parameters:
        df (DataFrame): DataFrame from PEAKS output with relevant columns.
        RT_tol (float): Minimum retention time difference to qualify as a candidate.
        ims_diff (float): Minimum fractional mobility difference to qualify.
    
    Returns:
        DataFrame: A filtered DataFrame of isomer candidates (one row per peptide).
    """
    # Ensure RT is numeric
    df["RT"] = pd.to_numeric(df["RT"], errors="coerce")
    
    # Compute mean 1/k0 from the '1/k0 Range' column
    df["1/k0"] = get_ims(df)
    
    # Compute fractional mobility variation per (Peptide, z)
    mob_frac_df = get_mob_frac_variation(df)
    df = df.merge(mob_frac_df, on=["Peptide", "z"], how="left")
    
    # Compute the retention time range for each unique peptide from the full dataset
    peptide_rt_df = get_peptide_rt_range(df)
    
    # Merge RT range information back into the main DataFrame
    df = df.merge(peptide_rt_df, on="Peptide", how="left")
    
    # Filter rows based on fractional mobility variation and RT difference thresholds
    candidates = df[(df["frac_variation"] > ims_diff) & (df["max_RT_diff"] > RT_tol)]
    
    # Drop duplicates to get one row per peptide candidate
    candidates = candidates.drop_duplicates(subset=["Peptide"]).reset_index(drop=True)
    keep_columns = ["Peptide", "z", "first_RT", "last_RT", "max_RT_diff", "Precursor Id","1/k0","frac_variation","m/z", "Mass", "Scan", "Accession"]
    candidates = candidates[keep_columns]
    
    return candidates
