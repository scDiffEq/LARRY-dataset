
from ._fetch._fetch_data_from_github import _fetch_data_from_github
from ._preprocess._Yeo2021_preprocessing_recipe import _Yeo2021_preprocessing_recipe

def _in_vitro(data_dir="./", return_obj=False, silent=True, fetch_kwargs={}, pp_kwargs={}):
    
    adata = _fetch_data_from_github("in_vitro", silent=silent, destination_dir=data_dir, **fetch_kwargs)
    adata = _Yeo2021_preprocessing_recipe(adata, data_dir, return_obj=False, **pp_kwargs)
    print("\n", adata)
    return adata