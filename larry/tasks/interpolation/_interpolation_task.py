
# -- import local dependencies: -----------------------------------------------
from ._interpolation_data import InterpolationData


# -- import packages: ---------------------------------------------------------
import ABCParse
import autodevice
import torch


# -- set typing: --------------------------------------------------------------
from typing import Tuple


# -- Operational class: -------------------------------------------------------
class InterpolationTask(ABCParse.ABCParse):
    def __init__(
        self,
        adata,
        time_key="Time point",
        use_key="X_pca",
        t0=2,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
        backend = "auto",
        silent = False,
        PCA = None,
        *args,
        **kwargs,
    ):
        self.__parse__(locals())

        self.data = InterpolationData(**self._DATA_KWARGS)
        
        from scdiffeq.core.lightning_models.base import SinkhornDivergence
        self._sinkhorn_fn = SinkhornDivergence
        
        self.SinkhornDivergence = self._sinkhorn_fn(**self._SINKHORN_KWARGS)

    @property
    def _DATA_KWARGS(self):
        return ABCParse.function_kwargs(
            func=InterpolationData, kwargs=self._PARAMS
        )
    @property
    def _SINKHORN_KWARGS(self):
        return ABCParse.function_kwargs(
            func=self._sinkhorn_fn, kwargs=self._PARAMS
        )

    def forward_without_grad(self, DiffEq):
        """Forward integrate over the model without gradients."""
        with torch.no_grad():
            X_hat = DiffEq.forward(self.data.X0, self.data.t)
            return self._parse_forward_out(X_hat)

    def forward_with_grad(self, DiffEq):
        """Forward integrate over the model retaining gradients."""
        torch.set_grad_enabled(True)
        X_hat = DiffEq.forward(self.data.X0, self.data.t)
        return self._parse_forward_out(X_hat)
    
    @property
    def potential(self):
        return "Potential" in str(self._DiffEq)
    
    def _parse_forward_out(self, X_hat):
        """to account for KLDiv"""
        if isinstance(X_hat, Tuple):
            return X_hat[0]
        return X_hat
    
    def _dimension_reduce_pca(self, X_hat):
        return torch.stack(
            [torch.Tensor(self._PCA.transform(x)) for x in X_hat.detach().cpu().numpy()]
        ).to(self.device)     


    def __call__(self, trainer, DiffEq, *args, **kwargs):
        
        self.__update__(locals())
        
        if self.potential:
            X_hat = self.forward_with_grad(DiffEq)
        else:
            X_hat = self.forward_without_grad(DiffEq)
            
        if not self._PCA is None:
            X_hat = self._dimension_reduce_pca(X_hat)
        
        d4_loss = self.SinkhornDivergence(X_hat[1], self.data.X_test_d4).item()
        d6_loss = self.SinkhornDivergence(X_hat[2], self.data.X_train_d6).item()

        if not self._silent:
            print(
                "- Epoch: {:<5}| Day 4 loss: {:.2f} | Day 6 loss: {:.2f}".format(
                    DiffEq.current_epoch, d4_loss, d6_loss,
                ),
            )

        return d4_loss, d6_loss
