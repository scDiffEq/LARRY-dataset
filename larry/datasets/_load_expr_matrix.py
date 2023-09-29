
import pathlib
import scipy

import ABCParse
from .. import utils

from typing import Union


class ExpressionMatrix(ABCParse.ABCParse):
    def __init__(self, mtx_path: Union[str, pathlib.Path]):

        self.__parse__(locals())

    def read_mtx(self):
        self.mtx = scipy.io.mmread(self.mtx_path)

    @property
    def npz_path(self):
        return pathlib.Path(str(self.mtx_path).replace("mtx.gz", "npz"))

    @property
    def npz_exists(self):
        return self.npz_path.exists()

    def save_mtx_as_npz(self):
        if not self.npz_exists:
            scipy.sparse.save_npz(self.npz_path, self.mtx)


def load_expr_matrix(path: Union[str, pathlib.Path]):
    
    expr_mat = ExpressionMatrix(mtx_path = path)

    if not expr_mat.npz_exists:
        print("Loading expression matrix as `.mtx`. This is slow. Will save as `.npz` for much faster future loading.")
        expr_mat.read_mtx()
        expr_mat.save_mtx_as_npz()
        
    mtx = scipy.sparse.load_npz(expr_mat.npz_path)

    return 
