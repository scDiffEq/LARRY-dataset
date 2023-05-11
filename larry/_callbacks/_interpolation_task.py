
__module_name__ = "_interpolation_task.py"
__doc__ = """Callback for passing the test set to measure interpolation of a withheld timepoint."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])


# -- import packages: ----------------------------------------------------------
from scdiffeq.core.lightning_models.base import SinkhornDivergence
from pytorch_lightning import Callback
from autodevice import AutoDevice
import numpy as np
import anndata
import torch


# -- Main callback class: ------------------------------------------------------
class InterpolationTask(Callback):
    def __init__(
        self,
        adata: anndata.AnnData,
        potential: bool = False,
        time_key: str = "Time point",
        lineage_key: str = "clone_idx",
        use_key: str = "X_pca",
        seed: int = 617,
        silent: bool = True,
        device = AutoDevice(),
    ):
        
        if isinstance(seed, int):
            torch.manual_seed(seed)
            np.random.seed(seed)

        self.adata = adata
        self.use_key = use_key
        self.time_key = time_key
        self.lineage_key = lineage_key
        self.Loss = SinkhornDivergence()
        self.potential = potential
        self.silent = silent
        self.device = device
        self.configure_eval_data()

    def configure_eval_data(self):

        self.df = self.adata.obs.copy()
        self.df_clonal = self.df.loc[self.df[self.lineage_key].notna()]
        self.t0_idx = (
            self.df_clonal.loc[self.df_clonal[self.time_key] == 2]
            .sample(10_000, replace=True)
            .index
        )
        self.X0 = torch.Tensor(self.adata[self.t0_idx].obsm[self.use_key]).to(self.device)
        self.t = torch.Tensor([2, 4, 6]).to(self.device)
        
        # test
        self.X_d4 = torch.Tensor(
            self.adata[self.df_clonal.loc[self.df_clonal[self.time_key] == 4].index].obsm[self.use_key]
        ).to(self.device)
        
        # train
        self.X_d6 = torch.Tensor(
            self.adata[self.df_clonal.loc[self.df_clonal[self.time_key] == 6].index].obsm[self.use_key]
        ).to(self.device)
        
    def forward_without_grad(self, DiffEq):
        """Forward integrate over the model without gradients."""
        with torch.no_grad():
            return DiffEq.forward(self.X0, self.t)
            
    def forward_with_grad(self, DiffEq):
        """Forward integrate over the model retaining gradients."""
        torch.set_grad_enabled(True)
        return DiffEq.forward(self.X0, self.t)
        

    def on_train_epoch_end(self, trainer, DiffEq):
        """Call"""
                
        if self.potential:
            X_hat = self.forward_with_grad(DiffEq)
        else:
            X_hat = self.forward_without_grad(DiffEq)
        
        d4_loss = self.Loss(X_hat[1], self.X_d4).item()
        d6_loss = self.Loss(X_hat[2], self.X_d6).item()
        
        if not self.silent:
            print("- Epoch: {:<5}| Day 4 loss: {:.2f} | Day 6 loss: {:.2f}".format(
                DiffEq.current_epoch,
                d4_loss,
                d6_loss,
            ),
                 )
            
        self.log("eval_d4_loss", d4_loss)
        self.log("eval_d6_loss", d6_loss)
