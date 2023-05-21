

__module_name__ = "_interpolation_task.py"
__doc__ = """Callback for passing the test set to measure interpolation of a withheld timepoint."""
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu"])

# -- import local dependencies: -----
from ._interpolation_data import InterpolationData
from ... import utils


# -- import packages: ----
from scdiffeq.core.lightning_models.base import SinkhornDivergence
import autodevice
import torch


class InterpolationTask(utils.ABCParse):
    def __init__(
        self,
        adata,
        time_key="Time point",
        use_key="X_pca",
        t0=2,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
        *args,
        **kwargs,
    ):
        self.__parse__(locals())

        self.data = InterpolationData(**self._DATA_KWARGS)
        self.SinkhornDivergence = SinkhornDivergence(**self._SINKHORN_KWARGS)

    @property
    def _DATA_KWARGS(self):
        return utils.extract_func_kwargs(
            func=InterpolationData, kwargs=self._PARAMS
        )
    @property
    def _SINKHORN_KWARGS(self):
        return utils.extract_func_kwargs(
            func=SinkhornDivergence, kwargs=self._PARAMS
        )

    def forward_without_grad(self, DiffEq):
        """Forward integrate over the model without gradients."""
        with torch.no_grad():
            return DiffEq.forward(self.X0, self.t)

    def forward_with_grad(self, DiffEq):
        """Forward integrate over the model retaining gradients."""
        torch.set_grad_enabled(True)
        return DiffEq.forward(self.X0, self.t)

    def __call__(self, trainer, DiffEq, *args, **kwargs):
        if self.potential:
            X_hat = self.forward_with_grad(DiffEq)
        else:
            X_hat = self.forward_without_grad(DiffEq)

        d4_loss = self.Loss(X_hat[1], self.data.X_test_d4).item()
        d6_loss = self.Loss(X_hat[2], self.data.X_train_d6).item()

        if not self.silent:
            print(
                "- Epoch: {:<5}| Day 4 loss: {:.2f} | Day 6 loss: {:.2f}".format(
                    DiffEq.current_epoch, d4_loss, d6_loss,
                ),
            )

        return d4_loss, d6_loss
