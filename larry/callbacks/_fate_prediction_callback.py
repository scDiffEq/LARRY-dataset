
# -- import packages: ----------------------------------------------------------
import ABCParse
import autodevice
import lightning
import pathlib
import pandas as pd

# -- import local dependencies: ------------------------------------------------
from .. import tasks

metrics = tasks.fate_prediction.metrics


# -- Callback: -----------------------------------------------------------------
class FatePredictionCallback(lightning.Callback, ABCParse.ABCParse):
    def __init__(
        self, adata, N: int = 2000, device=autodevice.AutoDevice(), *args, **kwargs
    ):
        self.__parse__(locals(), public=[None])

    def __repr__(self) -> str:
        return "FatePredictionCallback()"

    @property
    def F_obs(self):
        F_obs = tasks.fate_prediction.F_obs
        F_obs.index = F_obs.index.astype(str)
        return F_obs

    @property
    def t0_idx(self):
        return self.F_obs.index

    @property
    def _BASE_PATH(self):
        # assumes the first logger is the .CSVLogger
        base_path = pathlib.Path(self._log_dir).joinpath("fate_prediction_metrics")
        if not base_path.exists():
            base_path.mkdir()
        return base_path

    @property
    def _CKPT_METRICS_PATH(self):
        ckpt_metrics_path = self._BASE_PATH.joinpath(self._ckpt_name)
        if not ckpt_metrics_path.exists():
            ckpt_metrics_path.mkdir()
        return ckpt_metrics_path

    @property
    def _F_hat_unfiltered_csv_path(self):
        return self._CKPT_METRICS_PATH.joinpath("F_hat.unfiltered.csv")

    @property
    def _F_hat_processed_csv_path(self):
        return self._CKPT_METRICS_PATH.joinpath("F_hat.processed.csv")

    @property
    def _ACC_CSV_PATH(self):
        return self._CKPT_METRICS_PATH.joinpath("accuracy.csv")

    @property
    def _NEG_CE_CSV_PATH(self):
        return self._CKPT_METRICS_PATH.joinpath("neg_cross_entropy.csv")

    @property
    def _NM_CORR_CSV_PATH(self):
        return self._CKPT_METRICS_PATH.joinpath("neu_mon_corr.csv")

    def compute_fate_bias(self, pl_module):

        fate_bias = tasks.fate_prediction.FateBias()
        F_hat = fate_bias(
            adata=self._adata,
            DiffEq=pl_module.to(self._device),
            t0_idx=self.t0_idx,
            N=self._N,
        )
        F_hat.index = F_hat.index.astype(str)
        return F_hat

    def process_F_hat(self, F_hat: pd.DataFrame, index=None):

        print(f"Saving `F_hat` [ unfilt. ] to: {self._F_hat_unfiltered_csv_path}")
        F_hat.to_csv(self._F_hat_unfiltered_csv_path)

        F_hat = F_hat.drop("Undifferentiated", axis=1)
        F_hat = F_hat.div(F_hat.sum(1), axis=0)
        F_hat = F_hat.fillna(0)
        print(f"Saving `F_hat` [ processed ] to: {self._F_hat_processed_csv_path}")
        F_hat.to_csv(self._F_hat_processed_csv_path)

        if not index is None:
            F_hat.index = index
        else:
            F_hat.index = F_hat.index.astype(str)
        return F_hat

    def _accuracy(self, F_hat):

        accuracy_df = metrics.multi_idx_accuracy(self.F_obs, F_hat)
        print(f"Saving `accuracy_df` to: {self._ACC_CSV_PATH}")
        accuracy_df.to_csv(self._ACC_CSV_PATH)

    def _negative_cross_entropy(self, F_hat):

        neg_cross_entropy_df = metrics.multi_idx_negative_cross_entropy(
            self.F_obs, F_hat
        )
        print(f"Saving `neg_cross_entropy_df` to: {self._NEG_CE_CSV_PATH}")
        neg_cross_entropy_df.to_csv(self._NEG_CE_CSV_PATH)

    def _neu_mon_corr(self, F_hat):

        neu_mon_corr_df = pd.DataFrame(
            metrics.neutrophil_monocyte_correlation(self.F_obs, F_hat)
        )
        print(f"Saving `neu_mon_corr_df` to: {self._NM_CORR_CSV_PATH}")
        neu_mon_corr_df.to_csv(self._NM_CORR_CSV_PATH)

    def compute_metrics(self, F_hat):

        self._accuracy(F_hat)
        self._negative_cross_entropy(F_hat)
        self._neu_mon_corr(F_hat)

    def forward(self, pl_module, ckpt_name, log_dir):

        """ """

        self.__update__(locals())

        self.F_hat_unfiltered = self.compute_fate_bias(pl_module)
        self.F_hat = self.process_F_hat(self.F_hat_unfiltered)
        self.compute_metrics(self.F_hat)

    def on_train_end(self, trainer, pl_module):

        ckpt_name = f"on_train_end.epoch_{pl_module.current_epoch}"

        log_dir = self._log_dir = pl_module.loggers[0].log_dir

        self.forward(pl_module, ckpt_name=ckpt_name, log_dir=log_dir)

        
# -- import packages: ----------------------------------------------------------
# import autodevice
# import lightning
# import ABCParse

# -- import local dependencies: ------------------------------------------------
# from .. import utils
# from .. import tasks


# # -- Callback class: -----------------------------------------------------------
# class FatePredictionCallback(lightning.Callback, ABCParse.ABCParse):
#     def __init__(self, model):
        
#         self.__parse__(locals())
        
#         self._parse_model(model)

#     def _parse_model(self, model):
    
#         adata = model.adata
#         kNN_Graph = model.kNN_Graph
#         PCA = model.reducer.PCA
