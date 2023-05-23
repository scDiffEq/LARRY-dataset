
from ... import utils

import autodevice
import torch

from typing import Union
NoneType = type(None)

def fetch_data(
    adata: a.AnnData,
    use_key: str = "X_pca",
    torch: bool = True,
    groupby: Union[NoneType, str] = None,
    device: str = autodevice.AutoDevice(),
)

class FatePredictionData(utils.ABCParse):
    def __init__(
        self,
        adata,
        time_key = "Time point",
        N = 2000,
        groupby: Union[NoneType, str] = None,
        device: str = autodevice.AutoDevice(),
    ):
        
        self.__parse__(locals(), public = ['adata'])
        self.df = adata.obs.copy()
        
    @property
    def t(self):
        self._t = sorted(self.adata.obs[self._time_key].unique().tolist())
        return torch.Tensor(self._t).to(self._device)
        
    @property
    def t0_idx(self):
        return self.df.loc[self.df[self._time_key] == self._t[0]].index
    
    @property
    def t0_adata(self):
        return self.adata[self.t0_idx]
        
    @property
    def X0(self):
        return fetch_data(
            adata = self.t0_adata,
            use_key=self._use_key,
            torch = True,
            groupby = None,
            device: str = autodevice.AutoDevice(),
        )[None, :, :].expand(self._N, -1, -1)
        