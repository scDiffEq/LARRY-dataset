
from ... import utils
from . import _fate_prediction_utils as fate_utils

import pandas as pd
import autodevice
import torch
import time
import ABCParse
import adata_query

# from ._fate_prediction_data import FatePredictionData
from ._observed_fate_bias import F_obs

from typing import List, Optional, Union


### from when the `FatePredictionData` class was in it's own module script: 
### MOVE THIS TO _fate_bias_matrix_generation.py ?? and make one class. would make sense.
### then you have the task class which encompasses the modules for accuracy quantitation.
def is_notebook() -> bool:
    """
    https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook#:~:text=def%20is_notebook()%20%2D%3E,standard%20Python%20interpreter
    """
    
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter
    

class FatePredictionData(ABCParse.ABCParse):
    def __init__(
        self,
        adata,
        time_key = "Time point",
        use_key = "X_pca",
        N = 2000,
        t0_idx: Optional[Union[List, pd.Index]] = None,
        groupby: Optional[Union[None, str]] = None,
        device: str = autodevice.AutoDevice(),
    ):
        
        self.__parse__(locals(), public = ['adata'])
        self.df = adata.obs.copy()
        
    @property
    def _t(self):
        if not hasattr(self, "_t_list"):
            self._t_list = sorted(self.adata.obs[self._time_key].unique().tolist())
        return self._t_list
        
    @property
    def t(self):
        return torch.Tensor(self._t).to(self._device)
        
    @property
    def t0_idx(self):
        if not hasattr(self, "_t0_idx"):
            self._t0_idx = self.df.loc[self.df[self._time_key] == self._t[0]].index
        return self._t0_idx
    
    @property
    def t0_adata(self):
        if not hasattr(self, "_t0_adata"):
            self._t0_adata = self.adata[self.t0_idx]
        return self._t0_adata
        
    @property
    def X0(self):
        if not hasattr(self, "_X0"):
            
            self._X0 = adata_query.fetch(
                adata = self.t0_adata,
                key = self._use_key,
                torch = True,
                groupby = None,
                device = autodevice.AutoDevice(),
            )[:, None, :].expand(-1, self._N, -1)
            
#             self._X0 = fate_utils.fetch_data(
#                 adata = self.t0_adata,
#                 use_key=self._use_key,
#                 torch = True,
#                 groupby = None,
#                 device = autodevice.AutoDevice(),
#             )[:, None, :].expand(-1, self._N, -1)

        return self._X0


# -- main class: ----------------------------------------------------------------------      
class FateBias(ABCParse.ABCParse):
    def __init__(self, device=autodevice.AutoDevice()):
        
        self.__parse__(locals())
        
    @property
    def kNN(self):
        if not hasattr(self, "_kNN"):
            import scdiffeq as sdq
            self._kNN = sdq.tl.kNN(self._adata)
        return self._kNN
    
    @property
    def _MODEL_TYPE(self):
        return str(self._DiffEq)
    
    def _fate_bias_matrix(self):
        value_counts = {
            key: val[self._obs_key].value_counts() for key, val in self.F_hat.items()
        }
        # don't divide by N: integer counts are sometimes useful
        return pd.DataFrame(value_counts).T.fillna(0) #  / self._N

    def _predict_state(self, X0, t0_idx):
        """returns only the FINAL state"""
        return self._DiffEq(X0, t=self.data.t)[-1].detach().cpu().numpy()
    
    def umap_transform(self, umap_model, t0_idx=[], random=0):

        if random > 0:
            t0_idx = np.random.choice(self.t0_idx, random)

        self.X_umap = {}

        for i in tqdm.notebook.tqdm(range(len(t0_idx)), desc="umap projection"):
            self.X_umap[t0_idx[i]] = np.stack(
                [
                    umap_model.transform(x)
                    for x in sdq.tl.DataFormat(self.X_hat[t0_idx[i]]).to_numpy()
                ]
            )
    @property
    def _kNN_BASIS(self):
        return self.kNN.X_use.shape[1]
            
    def _configure_dimension_reduction(self)->bool:
        self._dimension_reduce = self.data.X0.shape[-1] > self._kNN_BASIS
    
    @property
    def _DIMENSION_REDUCE(self)->bool:
        if not hasattr(self, "_dimension_reduce"):
            self._configure_dimension_reduction()
        return self._dimension_reduce
    
    def _configure_t0_idx(self):
        
        if hasattr(self, "_t0_idx"):
            self._t0_idx_config = self._t0_idx
        else:
            try:
                self.F_obs = F_obs
            except:
                self.F_obs = F_obs(self.adata)
                self._t0_idx_config = self.F_obs().index
            finally:
                self._t0_idx_config = self.data.t0_idx
        
    @property
    def t0_idx(self):
        if not hasattr(self, "_t0_idx_config"):
            self._configure_t0_idx()
        return self._t0_idx_config
        
    def __call__(
        self,
        adata,
        DiffEq,
        obs_key: str = "Cell type annotation",
        use_key: str = "X_pca",
        time_key: str = "Time point",
        t0_idx = None,
        N=2000,
        kNN = None,
        PCA = None,
        notebook: bool = is_notebook(),
    ):
        """
        
        Generates F from simulations of input data and a model
        
        DiffEq and Graph are required.
        """

        self.__update__(locals())
        
        self.F_hat = {}
        
        self.data = FatePredictionData(
            self._adata,
            t0_idx = self._t0_idx,
            time_key = self._time_key,
            use_key = self._use_key,
            N = self._N,
            device = self._device,
        )
        
        if notebook:
            import tqdm.notebook
            TQDM = tqdm.notebook.tqdm
        else:
            import tqdm
            TQDM = tqdm.tqdm
            
        
        for i in TQDM(range(len(self.t0_idx)), desc="EVALUATION"):

            X_hat = self._predict_state(X0=self.data.X0[i].to(self._device), t0_idx=self.t0_idx[i])
            
            if self._DIMENSION_REDUCE:
                X_hat = self.PCA.transform(X_hat)
                
            self.F_hat[self.t0_idx[i]] = self.kNN.aggregate(
                X_hat, obs_key=self._obs_key, max_only=True
            )
            
        return self._fate_bias_matrix()
        

    def __repr__(self):
        return f"Fate prediction evaluation | evaluating: {self._MODEL_TYPE}"

