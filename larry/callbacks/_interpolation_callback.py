# -- import packages: ---------------------------------------------------------
import ABCParse
import autodevice
import lightning
import sklearn.decomposition

# -- import local dependencies: -----------------------------------------------
from .. import tasks

# -- set type hints: ----------------------------------------------------------
from typing import Optional

# -- callback cls: ------------------------------------------------------------
class InterpolationCallback(lightning.Callback, ABCParse.ABCParse):

    def __init__(
        self,
        model,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
        backend="auto",
        PCA: Optional[sklearn.decomposition._pca.PCA] = None,  # for VAE
        silent=False,
        *args,
        **kwargs,
    ) -> None:

        adata = model.adata
        use_key = model._use_key
        time_key = model.t_config.time_key
        t0 = model.t[0].item()
        # PCA = model.reducer.PCA

        self.__parse__(locals())

        self.interpolation_task = tasks.interpolation.InterpolationTask(
            **self._TASK_KWARGS
        )

    def __repr__(self) -> str:
        return "InterpolationCallback()"

    @property
    def _TASK_KWARGS(self):
        return ABCParse.function_kwargs(
            func=tasks.interpolation.InterpolationTask,
            kwargs=self._PARAMS,
        )

    def on_train_epoch_end(self, trainer, DiffEq, *args, **kwargs) -> None:
        """"""
        d4_loss, d6_loss = self.interpolation_task(trainer, DiffEq, *args, **kwargs)

        self.log("interpolation_eval_d4_loss", d4_loss)
        self.log("interpolation_eval_d6_loss", d6_loss)
