
# -- import packages: -----
import autodevice
import lightning


from .. import utils
from .. import tasks


class InterpolationCallback(lightning.Callback, utils.ABCParse):
    def __init__(
        self,
        model,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
        backend = "auto",
        silent = False,
        *args,
        **kwargs,
    ):
        
        adata = model.adata
        use_key = model._use_key
        time_key = model._time_key
        t0 = model.t[0].item()
        PCA = model.reducer.PCA

        self.__parse__(locals())

        self.interpolation_task = tasks.interpolation.InterpolationTask(**self._TASK_KWARGS)

    @property
    def _TASK_KWARGS(self):
        return utils.extract_func_kwargs(
            func=tasks.interpolation.InterpolationTask,
            kwargs=self._PARAMS,
        )

    def on_train_epoch_end(self, trainer, DiffEq, *args, **kwargs):
        """"""
        d4_loss, d6_loss = self.interpolation_task(trainer, DiffEq, *args, **kwargs)

        self.log("interpolation_eval_d4_loss", d4_loss)
        self.log("interpolation_eval_d6_loss", d6_loss)
