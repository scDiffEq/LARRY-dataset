from torch_adata import AnnDataset
from anndata import AnnData
from typing import List, Union

NoneType = type(None)


def to_AnnDataset(
    adata: AnnData,
    use_key: str = "X_pca",
    groupby: str = "Time point",
    obs_keys: List[str] = None,
    attr_names: {"obs": List[str], "aux": List[str]} = {"obs": [], "aux": []},
    one_hot: List[bool] = False,
    aux_keys: Union[List[str], NoneType] = None,
    silent: bool = False,
    sampling_weight_key: Union[str, NoneType] = None,
):

    return AnnDataset(**locals())
