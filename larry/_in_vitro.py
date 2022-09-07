
from ._fetch._fetch_data_from_github import _fetch_data_from_github
from ._preprocess._Yeo2021_preprocessing_recipe import _Yeo2021_preprocessing_recipe

def _in_vitro(return_obj=False, silent=True, **kwargs):
    adata = _fetch_data_from_github("in_vitro", silent=silent)
    adata = _Yeo2021_preprocessing_recipe(adata, return_obj=False, **kwargs)
    print("\n", adata)
    return adata