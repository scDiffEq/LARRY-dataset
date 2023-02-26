
# -- import packages: ----------------------------------------------------------
import anndata as a
import pandas as pd
import torch


# -- API-facing function: ------------------------------------------------------
def fetch_data(
    adata: a.AnnData, idx: pd.Index, n_sim: int = 1, use_key: str = "X_pca"
) -> torch.Tensor:

    """
    Fetch data as torch.Tensor using an index for adata.

    Parameters:
    -----------
    adata
        type: anndata.AnnData

    idx
        type: pd.Index

    n_sim:
        type: int
        default: 1

    use_key:
        type: str
        default: "X_pca"

    Returns:
    --------
    X_data
        type: torch.Tensor
    """

    n_dim = adata.obsm[use_key].shape[1]

    return torch.Tensor(adata[idx].obsm[use_key])[:, None, :].expand(
        len(idx), n_sim, n_dim
    )
