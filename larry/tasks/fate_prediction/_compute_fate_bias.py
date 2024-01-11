
import autodevice
import torch
import sklearn.decomposition
import pandas as pd
import anndata
import lightning

from typing import Optional


def compute_fate_bias(
    adata: anndata.AnnData,
    DiffEq: lightning.LightningModule,
    t0_idx: Optional[pd.Index] = None,
    obs_key: str = "Cell type annotation",
    use_key: str = "X_pca",
    time_key: str = "Time point",
    N: int = 2000,
    kNN: Optional["kNN"] = None,
    PCA: sklearn.decomposition.PCA = None,
    device: torch.device = autodevice.AutoDevice(),
    *args,
    **kwargs,
):
    """Evaluate fate bias

    Args:
        adata (anndata.AnnData): adata obj

        DiffEq (lightning.LightningModule):  Lightning neural differential equation

        t0_idx (Optional[pd.Index]): index of cells from which simulations should be run.
        If not provided, all t0 cells will be used. **Default** = None

        obs_key (str): description. **Default**

        use_key (str): description. **Default**

        time_key (str): description. **Default**

        N (int): number of trajectories to sample for each simulation. **Default** = 2000

        kNN ('kNN'): nearest nieghbors graph.

        PCA (sklearn.decomposition.PCA): **Default** = None

        device (torch.device): **Default** = autodevice.AutoDevice()

    Returns:
        F_hat (pd.DataFrame): predicted fate bias

    Notes:
        1. Default argument values are catered towards the LARRY dataset.
    """

    fate_bias = larry.tasks.fate_prediction.FateBias(device=device)
    function_kwargs = ABCParse.function_kwargs(func=fate_bias.__call__, kwargs=locals())
    return fate_bias(**function_kwargs)
