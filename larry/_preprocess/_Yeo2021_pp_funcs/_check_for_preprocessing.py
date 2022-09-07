import os
import anndata as a


from ._mk_pp_h5ad_path import _mk_pp_h5ad_path

def _check_for_preprocessing(self, adata):

    self.pp_complete = False
    pp_h5ad_path = _mk_pp_h5ad_path(adata)
    if os.path.exists(pp_h5ad_path):
        print("Preprocessing performed previously. Loading...", end="")
        adata = a.read_h5ad(pp_h5ad_path)
        print("done.")
        self.pp_complete = True

    self.adata = adata
    self.pp_h5ad_path = pp_h5ad_path