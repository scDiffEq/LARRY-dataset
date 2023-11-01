
# -- import packages: ----------------------------------------------------------
import autodevice
import lightning
import ABCParse

# -- import local dependencies: ------------------------------------------------
from .. import utils
from .. import tasks


# -- Callback class: -----------------------------------------------------------
class FatePredictionCallback(lightning.Callback, ABCParse.ABCParse):
    def __init__(self, model):
        
        self.__parse__(locals())
        
        self._parse_model(model)

    def _parse_model(self, model):
    
        adata = model.adata
        kNN_Graph = model.kNN_Graph
        PCA = model.reducer.PCA