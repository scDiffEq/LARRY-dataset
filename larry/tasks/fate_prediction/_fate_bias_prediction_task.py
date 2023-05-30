
# -- import local dependencies: -----
from ... import utils


# -- import packages: ----
import autodevice
import torch


NoneType = type(None)

class FateBiasPredictionTask(utils.ABCParse):
    _F_HAT_CONFIGURED = False
    _F_OBS_CONFIGURED = False
    def __init__(
        self,
        F_hat = None,
        F_obs = None,
        adata = None,
        DiffEq = None,
        kNN_Graph = None,
    ):
        self.__parse__(locals(), public = [None])
        
        
    def _configure_F_hat(self, ):
        
        if isinstance(self._F_hat, NoneType):
            self._F_hat = larry.tasks.fate_prediction.FateBias(
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
        if isinstance(self._F_obs, NoneType):
            self._F_obs = larry.tasks.fate_prediction.F_obs(self._adata).df
        self._F_OBS_CONFIGURED = True
        
    @property
    def F_obs(self):
        if not self._F_OBS_CONFIGURED:
            self._configure_F_obs()
        return self._F_obs
    
    
    @property
    def METRIC_precision_recall(self):
        multiclass_pr = larry.tasks.fate_prediction.metrics.MultiClassPrecisionRecall()
        return = multiclass_pr(self.F_obs, self.F_hat)
    
    
        