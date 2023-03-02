
# -- import packages: ----------------------------------------------------------
import anndata


# -- API-facing function: ------------------------------------------------------
def subset_clonal(adata: anndata.AnnData, lineage_key: str = "clone_idx"):

    """
    Parameters:
    -----------
    adata
        type: anndata.AnnData

    lineage_key
        type: str
        default: "clone_idx"

    Returns:
    --------
    clonal_adata
    """

    df = adata.obs.copy()
    clonal_adata = adata[df.loc[df[lineage_key].notna()].index].copy()
    clonal_adata.obs = clonal_adata.obs.reset_index()
    clonal_adata.obs.index = clonal_adata.obs.index.astype(str)

    return clonal_adata