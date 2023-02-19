
__module_name__ = "_get_lineage_obs.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import pandas
import anndata


# -----------------------------------------------------------------------------
def get_lineage_obs(
    adata: anndata.AnnData,
    lineage_idx: float,
    lineage_key: str = "clone_idx",
    obs_key: str = "Cell type annotation",
    time_key: str = "Time point",
) -> pandas.DataFrame:

    """
    Parameters:
    -----------
    adata
        type: anndata.AnnData

    lineage_idx
        type: float

    lineage_key
        type: str
        default: "clone_idx"

    obs_key
        type: str
        default: "Cell type annotation"

    time_key
        type: str
        default: "Time point"

    Returns:
    --------
    lineage_obs
        type: pandas.DataFrame
    """

    obs_df = adata.obs.copy()

    return (
        obs_df.loc[obs_df[lineage_key] == lineage_idx]
        .groupby(time_key)[obs_key]
        .value_counts()
        .unstack()
    )