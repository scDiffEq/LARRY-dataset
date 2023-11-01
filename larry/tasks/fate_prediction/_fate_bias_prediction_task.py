

# -- import packages: ---------------------------------------------------------
import autodevice
import torch
import ABCParse
import pandas as pd
import anndata


# -- import local dependencies: -----------------------------------------------
from ... import utils
from ._fate_bias_matrix_generation import FateBias
from ._observed_fate_bias import F_obs
from .metrics import MultiClassPrecisionRecall

# -- Operational class: -------------------------------------------------------
class FateBiasPredictionTask(ABCParse.ABCParse):
    _F_HAT_CONFIGURED = False
    _F_OBS_CONFIGURED = False
    def __init__(
        self,
        F_hat: pd.DataFrame = None,
        F_obs: pd.DataFrame = None,
        adata: anndata.AnnData = None,
        DiffEq = None,
        kNN_Graph = None,
    ):
        """
        F_hat: pd.DataFrame = None,
        F_obs: pd.DataFrame = None,
        adata: anndata.AnnData = None,
        DiffEq = None,
        kNN_Graph = None,
        """
        self.__parse__(locals(), public = [None])
        
        
    def _configure_F_hat(self, ):
        
        if self._F_hat is None:
            self._F_hat = FateBias(
                self._DiffEq,
                self._kNN_Graph,
            )(self._adata)
        self._F_HAT_CONFIGURED = True
        
    @property
    def F_hat(self):
        if not self._F_HAT_CONFIGURED:
            self._configure_F_hat()
        return self._F_hat
    
    def _configure_F_obs(self, ):
        if self._F_obs is None:
            self._F_obs = F_obs(self._adata).df
        self._F_OBS_CONFIGURED = True
        
    @property
    def F_obs(self):
        if not self._F_OBS_CONFIGURED:
            self._configure_F_obs()
        return self._F_obs
    
    
    @property
    def METRIC_precision_recall(self):
        multiclass_pr = MultiClassPrecisionRecall()
        return = multiclass_pr(self.F_obs, self.F_hat)
    
    
        