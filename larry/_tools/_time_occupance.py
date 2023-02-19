
# -- import packages: --------------------------------------------------------------------
import anndata
import numpy as np
import pandas as pd


# -- define types: -----------------------------------------------------------------------
from typing import Union
NoneType = type(None)


# -- supporting function(s): -------------------------------------------------------------
def _has_t(df, time_key, query_t: list):
    """return if ALL specified values are contained in passed list"""
    return np.array([t in df[time_key].unique() for t in query_t])
    
def filter_fate(
    adata,
    fate_time=[4, 6],
    fate=["undiff"],
    time_key="Time point",
    state_key="Annotation",
):
    df = adata.obs.copy()
    
    idx = ~(df[time_key].isin(fate_time) & (df[state_key].isin(fate)))
    return adata[idx].copy()
    
    
# -- main function: ----------------------------------------------------------------------
def time_occupance(
    adata: anndata.AnnData,
    lineage_key: str = "clone_idx",
    exclude_fate: tuple = ("Annotation", ["undiff"]),
    fate_time=[4, 6],
    time_key: str = "Time point",
    return_df=False,
) -> Union[pd.DataFrame, None]:

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
    
    if not isinstance(exclude_fate, NoneType):
        
        adata_ = filter_fate(
            adata,
            fate_time=fate_time,
            fate=exclude_fate[1],
            time_key=time_key,
            state_key=exclude_fate[0],
        )
    
    else:
        adata_ = adata
        
    df = adata_.obs.copy()
    t_query = sorted(df[time_key].unique())
    
    grouped_lineages = df.dropna(subset=[lineage_key]).groupby(lineage_key)

    adata.uns["time_occupance"] = time_occ_df = pd.DataFrame.from_dict(
        grouped_lineages.apply(_has_t, time_key, t_query).to_dict(),
        orient="index",
        columns=t_query,
    )
    msg = "- [NOTE] | Added lineage-time occupance to: adata.uns['time_occupance']"
    print(msg)

    if return_df:
        return time_occ_df