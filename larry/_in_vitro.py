
from ._fetch._fetch_data_from_github import _fetch_data_from_github
from ._preprocess._Yeo2021_preprocessing_recipe import _Yeo2021_preprocessing_recipe

from ._preprocess._add_data_from_supp_files import _add_data_from_supp_files
from ._preprocess._build_kNN import _build_annoy_adata

def _in_vitro(destination_dir="./", return_obj=False, silent=True, fetch_kwargs={}, pp_kwargs={}):
    
    adata = _fetch_data_from_github(
        "in_vitro",
        silent=silent,
        destination_dir=destination_dir,
        **fetch_kwargs,
    )
    adata = _Yeo2021_preprocessing_recipe(
        adata,
        destination_dir,
        return_obj=False,
        **pp_kwargs,
    )
    adata.uns['data_dir'] = destination_dir
    _add_data_from_supp_files(adata)
    _build_annoy_adata(adata)
    
    print("\n", adata)
    return adata