
import pandas as pd


def row_norm_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sum normalize a pandas DataFrame across rows.
    
    Parameters
    ----------
    df: pd.DataFrame
        Input (raw, un-normalized) df
    
    Returns
    -------
    df_norm: pd.DataFrame
        Row-sum normalized df
    """
    return df.div(df.sum(axis=1), axis=0)


def col_norm_df(df) -> pd.DataFrame:
    """
    Sum normalize a pandas DataFrame over columns.
    
    Parameters
    ----------
    df: pd.DataFrame
        Input (raw, un-normalized) df
    
    Returns
    -------
    df_norm: pd.DataFrame
        Column-sum normalized df
    """
    return df.div(df.sum(axis=0), axis=1)
