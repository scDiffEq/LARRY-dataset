
__module_name__ = "_get_annotated_metadata_d2_lineage_cells.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


# import packages: ------------------------------------------------------------
import anndata
import pandas


# primary function: -----------------------------------------------------------
def get_annotated_metadata_d2_lineage_cells(
    adata: anndata.AnnData, return_df: bool = False
) -> pandas.DataFrame:

    """
    Get the annotated metadata table for d2 lineage cells

    Parameters:
    -----------
    adata

    Returns:
    --------
    d2_lin_cells_ann_meta
        type: pandas.DataFrame
    """

    meta_df = adata.obs.copy()
    is_d2 = meta_df["Time point"] == 2
    d2_lin_mask = adata.uns["d2_lin_mask"]
    major_fate_df = adata.uns["major_fate_df"]
    growth_rates = adata.uns['GrowthRateDict']

    d2_lin_cells_ann_meta = meta_df[is_d2][d2_lin_mask].merge(
        major_fate_df, on="clone_idx", how="left"
    )    
    
    for k, v in growth_rates.items():
        d2_lin_cells_ann_meta[k] = v
        
    adata.uns['d2_lineage_cells_metadata'] = d2_lin_cells_ann_meta
    
    if return_df:
        return d2_lin_cells_ann_meta