
# -- import packages: -----
import autodevice
import lightning


from .. import utils
from .. import tasks


class InterpolationCallback(utils.ABCParse):
    def __init__(
        self,
        adata,
        time_key="Time point",
        use_key="X_pca",
        t0=2,
        n_samples=10_000,
        lineage_key="clone_idx",
        device=autodevice.AutoDevice(),
    ):
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
