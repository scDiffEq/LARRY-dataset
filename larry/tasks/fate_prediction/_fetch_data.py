import licorice_font as lf
import pandas as pd
import scipy.sparse
import numpy as np
import anndata as a
import autodevice
import torch


from ... import utils


from typing import Union, Dict, List
NoneType = type(None)


class DataFormat(utils.ABCParse):
    def __init__(self, data: Union[torch.Tensor, np.ndarray]):
        self.__parse__(locals(), public=["data"])

    @property
    def device_type(self):
        if hasattr(self.data, "device"):
            return self.data.device.type
        return "cpu"

    @property
    def is_ArrayView(self):
        return isinstance(self.data, a._core.views.ArrayView)

    @property
    def is_numpy_array(self):
        return isinstance(self.data, np.ndarray)

    @property
    def is_torch_Tensor(self):
        return isinstance(self.data, torch.Tensor)

    @property
    def on_cpu(self):
        return self.device_type == "cpu"

    @property
    def on_gpu(self):
        return self.device_type in ["cuda", "mps"]

    def to_numpy(self):
        if self.is_torch_Tensor:
            if self.on_gpu:
                return self.data.detach().cpu().numpy()
            return self.data.numpy()
        elif self.is_ArrayView:
            return self.data.toarray()
        return self.data

    def to_torch(self, device=autodevice.AutoDevice()):
        self.__update__(locals())
        if self.is_torch_Tensor:
            return self.data.to(device)
        elif self.is_ArrayView:
            self.data = self.data.toarray()
        return torch.Tensor(self.data).to(device)


class X_use(ABCParse):
    _idx = None
    """The matrix you want from AnnData"""

    def __init__(self, adata: a.AnnData, use_key: str = "X_pca"):
        self.__parse__(locals(), public=[None])

    @property
    def _DTYPES(self):
        """array of data types over which to enumerate"""
        return np.array([NoneType, np.ndarray, scipy.sparse.spmatrix, torch.Tensor])

    def __is__(self, X, dtype) -> bool:
        """Boolean gate to check data type"""
        return isinstance(X, dtype)

    @property
    def grouped(self):
        return self._adata.obs.copy().groupby(self._groupby)

    @property
    def dtype(self):
        """returns data type of X"""
        return self._DTYPES[[self.__is__(self.X, dtype) for dtype in self._DTYPES]][0]

    @property
    def _LAYER_KEYS(self):
        """returns keys of adata.layers"""
        return list(self._adata.layers)

    @property
    def _OBSM_KEYS(self):
        """returns keys of adata.obsm"""
        return self._adata.obsm_keys()

    def __repr__(self):
        use_key = lf.font_format(self._use_key, ["BOLD"])
        return f"X_use: {use_key}"

    def _fetch_X(self, idx=None):

        if isinstance(idx, NoneType):
            adata = self._adata
        else:
            adata = self._adata[idx]

        if self._use_key in self._OBSM_KEYS:
            return adata.obsm[self._use_key]

        if self._use_key in self._LAYER_KEYS:
            return adata.layers[self._use_key]

        if self._use_key == "X":
            return adata.X

    @property
    def X(self):
        return self._fetch_X(self._idx)

    def _ensure__dense_array(self):
        if self.dtype == scipy.sparse.spmatrix:
            return self.X.toarray()
        return self.X

    def forward(self, df, device=None) -> Union[torch.Tensor, np.ndarray]:

        if not isinstance(df, NoneType):
            self._idx = df.index

        data_format = DataFormat(data=self._ensure__dense_array())

        if self._torch:
            return data_format.to_torch(device=device)
        return data_format.to_numpy()

    def groupby_forward(self) -> Dict:

        X_dict = self.grouped.apply(self.forward, device="cpu").to_dict()
        if self._torch:
            return {key: val.to(self._device) for key, val in X_dict.items()}
        return X_dict

    def __call__(
        self,
        torch: bool = False,
        groupby: Union[NoneType, str] = None,
        device: str = autodevice.AutoDevice(),
    ):

        self.__update__(locals(), private=["torch", "groupby", "device"])

        if not isinstance(groupby, NoneType):
            return self.groupby_forward()

        return self.forward(df=None, device=self._device)


def fetch_data(
    adata: a.AnnData,
    use_key: str = "X_pca",
    torch: bool = True,
    groupby: Union[NoneType, str] = None,
    device: str = autodevice.AutoDevice(),
) -> Union[Dict[torch.Tensor, np.ndarray], torch.Tensor, np.ndarray]:

    xu = X_use(adata=adata, use_key=use_key)
    return xu(torch=torch, groupby=groupby, device=device)
