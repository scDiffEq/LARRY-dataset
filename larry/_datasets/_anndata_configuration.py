
# -- import packages: ----------------------------------------------------------
from pathlib import Path
import pandas as pd
import numpy as np
import anndata
import scipy
import scipy.io


# -- import local dependencies: ------------------------------------------------
from .._utils import AutoParseBase
from ._dataset_utils import Messages

# -- supporting functions and classes: -----------------------------------------
def _read_format_df(path, **kwargs):
    df = pd.read_csv(path, **kwargs)
    df.index = df.index.astype(str)
    return df

class ExpressionMatrix(AutoParseBase):
    """Container for parsing the expression matrix."""

    def __init__(self, mtx_path: str):

        self.__parse__(locals())

    def read_mtx(self)->None:
        self.mtx = scipy.io.mmread(self.mtx_path)

    @property
    def npz_path(self)->Path:
        return Path(str(self.mtx_path).replace("mtx.gz", "npz"))

    @property
    def npz_exists(self)->bool:
        return self.npz_path.exists()

    def save_mtx_as_npz(self)->None:
        if not self.npz_exists:
            scipy.sparse.save_npz(self.npz_path, self.mtx)

    def __call__(self) -> scipy.sparse.csr_matrix:

        """
        Load the expression matrix.

        Parameters:
        -----------
        None

        Returns:
        --------
        X
            type: scipy.sparse.csr_matrix

        Given a path to a .mtx file, initializes the ExpressionMatrix class.
        Looks for the .npz file. If it's not there, the .mtx file is parsed
        and saved as .npz.
        """
        if not self.npz_exists:
            messages.mtx_to_npz()
            self.read_mtx()
            self.save_mtx_as_npz()

        return scipy.sparse.csr_matrix(scipy.sparse.load_npz(self.npz_path))


class AnnDataConfiguration(AutoParseBase):
    """Construct AnnData from constituent components."""

    def __init__(self, X_path, obs_path, var_path, X_clone_path, silent=False):

        """Pass paths of adata components."""

        self.__parse__(locals(), public=[None])
        self._expr_mtx = ExpressionMatrix(mtx_path=self._X_path)

    @property
    def X(self)->scipy.sparse.csr_matrix:
        """returns expression matrix"""
        if not hasattr(self, "_X"):
            self._X = self._expr_mtx()
        return self._X

    @property
    def var(self)->pd.DataFrame:
        """returns var pd.DataFrame"""
        if not hasattr(self, "_var_df"):
            self._var_df = _read_format_df(
                self._var_path, header=None, names=["gene_ids"]
            )
        return self._var_df

    @property
    def obs(self)->pd.DataFrame:
        """returns obs pd.DataFrame"""
        if not hasattr(self, "_obs_df"):
            self._obs_df = _read_format_df(self._obs_path, sep="\t")
        return self._obs_df

    @property
    def X_clone(self)->scipy.sparse.csr_matrix:
        """returns cell x clonal barcode matrix"""
        if not hasattr(self, "_X_clone"):
            self._X_clone = scipy.sparse.csr_matrix(scipy.io.mmread(self._X_clone_path))
            # make sure it is oriented correctly
            if self._X_clone.shape[0] != self.X.shape[0]:
                self._X_clone = self._X_clone.T
        return self._X_clone

    @property
    def adata(self):
        """returns formatted AnnData object"""
        if not hasattr(self, "_adata"):
            self._adata = anndata.AnnData(
                X=self.X,
                dtype=self.X.dtype,
                obs=self.obs,
                var=self.var,
                obsm={"X_clone": self.X_clone},
            )
        return self._adata
