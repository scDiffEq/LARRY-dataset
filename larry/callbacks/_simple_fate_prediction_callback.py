import ABCParse
import autodevice
import lightning
import pathlib
import pandas as pd

from .. import tasks

class SimpleFatePredictionCallback(lightning.Callback, ABCParse.ABCParse):
    def __init__(self, adata, *args, **kwargs):
        self.__parse__(locals())

    @property
    def F_obs(self):
        return tasks.fate_prediction.F_obs

    @property
    def t0_idx(self):
        return self.F_obs.index

    def _compute_fate_bias(self, pl_module):

        fate_bias = tasks.fate_prediction.FateBias()
        F_hat = fate_bias(
            adata=self._adata,
            DiffEq=pl_module.to(autodevice.AutoDevice()),
            t0_idx=self.t0_idx,
        )
        F_hat.index = F_hat.index.astype(str)
        return F_hat
    
    def compute_metrics(self, pl_module, F_hat):
        
        F_obs = self.F_obs.copy()
        
        F_hat = F_hat.loc[F_obs.index]
        
        BASE_PATH = pathlib.Path(pl_module.loggers[0].log_dir)
        
        F_hat_csv_path = BASE_PATH.joinpath("F_hat.csv")
        print(f"Saving `F_hat` to: {F_hat_csv_path}")
        F_hat.to_csv(F_hat_csv_path)
        
        F_obs_csv_path = BASE_PATH.joinpath("F_obs.csv")
        print(f"Saving `F_obs` to: {F_obs_csv_path}")
        F_obs.to_csv(F_obs_csv_path)
        
        accuracy_df = tasks.fate_prediction.metrics.multi_idx_accuracy(F_obs, F_hat)
        
        accuracy_csv_path = BASE_PATH.joinpath("accuracy.csv")
        print(f"Saving `accuracy_df` to: {accuracy_csv_path}")
        accuracy_df.to_csv(accuracy_csv_path)
        
        neg_cross_entropy_df = (
            tasks.fate_prediction.metrics.multi_idx_negative_cross_entropy(F_obs, F_hat)
        )
        neg_cross_entropy_csv_path = BASE_PATH.joinpath("neg_cross_entropy.csv")
        print(f"Saving `neg_cross_entropy_df` to: {neg_cross_entropy_csv_path}")
        neg_cross_entropy_df.to_csv(neg_cross_entropy_csv_path)
        
        neu_mon_corr_df = tasks.fate_prediction.metrics.neutrophil_monocyte_correlation(F_obs, F_hat)
        neu_mon_corr_csv_path = BASE_PATH.joinpath("neu_mon_corr.csv")
        print(f"Saving `neu_mon_corr_df` to: {neu_mon_corr_csv_path}")
        pd.DataFrame(neu_mon_corr_df).to_csv(neu_mon_corr_csv_path)
        

    def on_train_end(self, trainer, pl_module):
        
        F_hat = self._compute_fate_bias(pl_module)
        
        self.compute_metrics(pl_module, F_hat)
        
        
