
# -- import packages: --------------------------------------------------------------------
import anndata
import numpy as np
import pandas as pd


# -- supporting function(s): -------------------------------------------------------------
def _has_t(df, time_key, query_t: list):
    """return if ALL specified values are contained in passed list"""
    return np.array([t in df[time_key].unique() for t in query_t])


# -- main function: ----------------------------------------------------------------------
def _time_occupance(
    adata: anndata.AnnData, lineage_key: str = "clone_idx", time_key: str = "Time point"
) -> pd.DataFrame:

    """
    Parameters:
    -----------
    adata
        type: anndata.AnnData

    lineage_key
        type: str

    time_key
        type: str

    Returns:
    --------
    time_occupance
        type: pandas.DataFrame
    """

    df = adata.obs.copy()
    t_query = sorted(df[time_key].unique())

    grouped_lineages = df.dropna(subset=[lineage_key]).groupby(lineage_key)

    return pd.DataFrame.from_dict(
        grouped_lineages.apply(_has_t, time_key, t_query).to_dict(),
        orient="index",
        columns=t_query,
    )

def time_occupance(adata: anndata.AnnData, lineage_key: str = "clone_idx", time_key: str = "Time point", return_df=False):
    
    time_occ_df = _time_occupance(adata, lineage_key, time_key)
    
    adata.uns['time_occupance'] = time_occ_df
    
    if return_df:
        return time_occ_df
