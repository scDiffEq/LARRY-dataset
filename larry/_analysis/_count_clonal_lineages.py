

# import packages: ------------------------------------------------------------
import pandas as pd
from anndata import AnnData


# -----------------------------------------------------------------------------
def count_clonal_lineages(
    adata: AnnData, lineage_key: str = "clone_idx", groupby_key: str = "Time point", return_df: bool = False
) -> pd.DataFrame:

    """
    Count clonal lineages over groups (e.g., time points).

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    lineage_key
        type: str
        default: "clone_idx"

    groupby_key
        type: str
        default: "Time point"

    Returns:
    --------
    count_df

    Notes:
    ------
    """

    meta_df = adata.obs.copy()
    meta_clonal = meta_df.dropna()

    lineage_grouped = meta_clonal.groupby(lineage_key)
    count_df = (
        lineage_grouped[groupby_key].value_counts().unstack().fillna(0).astype(int)
    )
    
    adata.uns['lineage_count_df'] = count_df
    
    if return_df:
        return count_df