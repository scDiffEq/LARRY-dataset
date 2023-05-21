
# -- import packages: ----------------------------------------------------------
from autodevice import AutoDevice
import torch_adata
import anndata as a
import pandas as pd
import torch


# -- set typing: ---------------------------------------------------------------
NoneType = type(None)

def fetch(adata, use_key: str, device=AutoDevice()):
    return torch_adata.tl.fetch(adata, use_key).to(device)


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


def fetch_format(adata, use_key: str, idx=None, N=1, device=AutoDevice()):

    if isinstance(idx, NoneType):
        idx = range(len(adata))

    return fetch(adata[idx], use_key=use_key, device=device)[:, None].expand(
        -1, N, -1
    )



from ._data_format import DataFormat

class Fetch:
    def __init__(self, adata, use_key: str, device=autodevice.AutoDevice()):
        ...
        
def fetch(adata: anndata.AnnData, use_key: str)->torch.Tensor:
    
    f = Fetch(adata)
    return f.X(use_key)
            

def fetch(adata, use_key: str, device=AutoDevice()):
    return torch_adata.tl.fetch(adata, use_key).to(device)


# -- API-facing function: ------------------------------------------------------
def fetch_data(
    adata: a.AnnData, idx: pd.Index, n_sim: int = 1, use_key: str = "X_pca"
) -> torch.Tensor:

    n_dim = adata.obsm[use_key].shape[1]

    return torch.Tensor(adata[idx].obsm[use_key])[:, None, :].expand(
        len(idx), n_sim, n_dim
    )


def fetch_format(adata, use_key: str, idx=None, N=1, device=AutoDevice()):

    if isinstance(idx, NoneType):
        idx = range(len(adata))

    return fetch(adata[idx], use_key=use_key, device=device)[:, None].expand(
        -1, N, -1
    )
