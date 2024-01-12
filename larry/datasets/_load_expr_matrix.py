
# -- import packages: ---------------------------------------------------------
import ABCParse
import pathlib
import scipy


# -- set typing: --------------------------------------------------------------
from typing import Union


# -- operational class: -------------------------------------------------------
class RetypedExpressionMatrix(ABCParse.ABCParse):
    def __init__(self, mtx_path: Union[str, pathlib.Path]):
        """"""

        self.__parse__(locals())

    @property
    def mtx_path(self):
        return pathlib.Path(self._mtx_path)

    @property
    def npz_path(self):
        return pathlib.Path(str(expr_mat.mtx_path).replace("mtx.gz", "npz"))

    @property
    def mtx(self):
        """"""
        return scipy.io.mmread(self.mtx_path)

    @property
    def _load_msg(self) -> str:
        """"""
        return "Loading expression matrix (`.mtx.gz`) and saving as `.npz` for much faster future loading."

    def save_mtx_as_npz(self):
        """"""
        if not self.npz_path.exists():
            self._INFO(self._load_msg)
            scipy.sparse.save_npz(self.npz_path, self.mtx)

    def load_npz(self) -> scipy.sparse.coo_matrix:
        """Load expression matrix into memory.

        Args:
            None

        Returns:
            scipy.sparse.coo_matrix
        """
        return scipy.sparse.load_npz(self.npz_path)


# -- API-facing function: -----------------------------------------------------
def load_expr_matrix(mtx_path: Union[str, pathlib.Path]):
    """Load the expression matrix as .npz (downloaded as .mtx)

    Args:
        mtx_path (Union[str, pathlib.Path])

    Returns:
        scipy.sparse.coo_matrix

    To-Do:
        Check whether coo is disadvantageous to csr or csc.
    """

    expr_mat = RetypedExpressionMatrix(mtx_path=mtx_path)
    expr_mat.save_mtx_as_npz()

    return expr_mat.load_npz()
