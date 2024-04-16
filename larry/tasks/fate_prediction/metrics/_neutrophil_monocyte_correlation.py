
# -- import packages: ---------------
import ABCParse
from typing import List
import scipy.stats
import pandas as pd

# -- import local dependencies: -----
from .._in_vitro_fate_prediction_subsets import InVitroFatePredictionSubsets

# -- operational class: -------------
class NeuMonCorr(ABCParse.ABCParse):
    def __init__(
        self,
        fate: str = "Neutrophil",
        possible_fates: List[str] = ["Neutrophil", "Monocyte"],
        corr_func=scipy.stats.pearsonr,
        *args,
        **kwargs,
    ):
        """"""
        self.__parse__(locals())
        
        assert fate in possible_fates
        self._indices = InVitroFatePredictionSubsets().indices
        
    @property
    def train_early_nm_idx(self):
        return self._indices['unique_train']['N/M-early']
        
    @property
    def test_early_nm_idx(self):
        return self._indices['unique_test']['N/M-early']
        
    @property
    def train_nm_idx(self):
        return self._indices['unique_train']['N/M']
        
    @property
    def test_nm_idx(self):
        return self._indices['unique_test']['N/M']
        
    @property
    def n_fates(self):
        return len(self._possible_fates)
    
    @property
    def na_fill(self):
        """What nan values should be filled with, effectively substituting a +1 rule"""
        return 1 / self.n_fates
        
    def row_normalize(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Sum normalize a pandas DataFrame across rows.

        Parameters
        ----------
        df: pd.DataFrame
            Input (raw, un-normalized) df

        Returns
        -------
        df_norm: pd.DataFrame
            Row-sum normalized df
        """
        return df.div(df.sum(axis=1), axis=0)
    
    @property
    def true_fate_bias(self):
        if not hasattr(self, "_true_fate_bias"):
            self._true_fate_bias = self.row_normalize(self._F_obs[self._possible_fates]).fillna(self.na_fill)[self._fate]
        return self._true_fate_bias
        
    @property
    def pred_fate_bias(self):
        if not hasattr(self, "_pred_fate_bias"):
            self._pred_fate_bias = self.row_normalize(self._F_hat[self._possible_fates]).fillna(self.na_fill)[self._fate]
        return self._pred_fate_bias
    
    def correlation(self, index):
        
        F_obs = self.true_fate_bias
        F_hat = self.pred_fate_bias
        
        return self._corr_func(F_obs.loc[index], F_hat.loc[index])
    
    @property
    def _ATTRS(self):
        return [key for key in self.__dir__() if "nm_idx" in key]
    
    def _to_dict(self, result):
        return {
            "correlation": result.correlation,
            "statistic": result.statistic,
            "pvalue": result.pvalue,
            "CI_low": result.confidence_interval().low,
            "CI_high": result.confidence_interval().high,
        }

    def forward(self, attr):
        ix = getattr(self, attr)
        return self._to_dict(self.correlation(ix))
    
    def _to_frame(self, Results):
        return pd.DataFrame(Results).T
    
    def __call__(
        self,
        F_obs: pd.DataFrame,
        F_hat: pd.DataFrame,
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        
        """
        Args:
            F_obs (pd.DataFrame): true fate bias matrix.

            F_hat (pd.DataFrame): true fate bias matrix.
            
        Returns:
            (pd.DataFrame)
        """

        self.__parse__(locals())
        
        Results = {}
        for attr in self._ATTRS:
            attr_key = attr.split("_idx")[0]
            Results[attr_key] = self.forward(attr)
        return self._to_frame(Results)
    
    

# -- API-facing function: --------------------------------------------            
def neutrophil_monocyte_correlation(
    F_obs: pd.DataFrame,
    F_hat: pd.DataFrame,
    *args,
    **kwargs,
) -> pd.DataFrame:
    """
    Args:
        F_obs (pd.DataFrame): true fate bias matrix.
        
        F_hat (pd.DataFrame): true fate bias matrix.
    """
    nm_corr = NeuMonCorr()
    return nm_corr(F_obs, F_hat)


# import ABCParse
# from typing import List
# import scipy.stats

# from .... import utils

# class NeuMonCorr(ABCParse.ABCParse):
#     def __init__(
#         self,
#         key: str = "Neutrophil",
#         reference: List[str] = ["Neutrophil", "Monocyte"],
#         corr_func=scipy.stats.pearsonr,
#         *args,
#         **kwargs,
#     ):
#         self.__parse__(locals())

#     def forward(self, df):
#         return utils.row_norm_df(df[self._reference]).fillna(
#             1 / len(self._reference)
#         )[self._key]

#     @property
#     def correlation(self):
#         return self._corr_func(self._F_obs_, self._F_hat_)

#     def __call__(self, F_obs, F_hat, *args, **kwargs):

#         self.__parse__(locals())

#         self._F_obs_ = self.forward(self._F_obs)
#         self._F_hat_ = self.forward(self._F_hat)

#         return self.correlation
    
# def neutrophil_monocyte_correlation(F_obs, F_hat):
#     nm_corr = NeuMonCorr()
#     return nm_corr(F_obs, F_hat)
