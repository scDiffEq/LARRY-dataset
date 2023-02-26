
__module_name__ = "_count_clones_by_timepoint.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import anndata
import pandas


# -----------------------------------------------------------------------------
def count_clones_by_timepoint(
    adata: anndata.AnnData,
    lineage_key: str = "clone_idx",
    time_key: str = "Time point",
    key_added: str = "lineage_time_counts",
    return_df: bool = False,
) -> pandas.DataFrame:

    """
    Annotates adata with pandas.DataFrame: adata.uns['clone_x_timepoint'].
    
    Parameters:
    -----------
    adata
        type: anndata.AnnData
        
    lineage_key
        type: str
        default: "clone_idx"
    time_key
        type: str
        default: "Time point"
    
    return_df
        type: bool
        default: False
    
    Returns:
    --------
    None or clone_x_timepoint
        type: NoneType or pandas.DataFrame
    """

    obs_df = adata.obs.copy()

    clone_x_timepoint = (
        obs_df.dropna()
        .groupby("Time point")["clone_idx"]
        .value_counts()
        .to_frame()
        .unstack()
        .fillna(0)
    ).T
    
    clone_x_timepoint = clone_x_timepoint.reset_index(drop=True)
    clone_x_timepoint.columns.name = None

    adata.uns["clone_x_timepoint"] = clone_x_timepoint

    if return_df:
        return clone_x_timepoint