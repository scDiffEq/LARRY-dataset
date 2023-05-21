
from ... import utils

import autodevice
import torch


class FatePredictionData(utils.ABCParse):
    def __init__(self, adata):
        
        self.__parse__(locals())