import ABCParse
from typing import List
import scipy.stats

from .... import utils

class NeuMonCorr(ABCParse.ABCParse):
    def __init__(
        self,
        key: str = "Neutrophil",
        reference: List[str] = ["Neutrophil", "Monocyte"],
        corr_func=scipy.stats.pearsonr,
        *args,
        **kwargs,
    ):
        self.__parse__(locals())

    def forward(self, df):
        return utils.row_norm_df(df[self._reference]).fillna(
            1 / len(self._reference)
        )[self._key]

    @property
    def correlation(self):
        return self._corr_func(self._F_obs_, self._F_hat_)

    def __call__(self, F_obs, F_hat, *args, **kwargs):

        self.__parse__(locals())

        self._F_obs_ = self.forward(self._F_obs)
        self._F_hat_ = self.forward(self._F_hat)

        return self.correlation
    
def neutrophil_monocyte_correlation(F_obs, F_hat):
    nm_corr = NeuMonCorr()
    return nm_corr(F_obs, F_hat)