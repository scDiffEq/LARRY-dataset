
# -- import packages: ----------------------------------------------------------
import autodevice
import lightning


# -- import local dependencies: ------------------------------------------------
from .. import utils
from .. import tasks


# -- Callback class: -----------------------------------------------------------
class FatePredictionCallback(utils.ABCParse):
    def __init__(self, adata):
        ...
        