import scipy
import numpy as np


def _sparse_variance(X, tmp, axis=0):
    """Calculate variance across the specified axis of a sparse matrix"""

    mean_gene = X.mean(axis=axis).A.squeeze()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene**2


def _sparse_row_wise_multiplication(X, scalar):
    """Multiply each row of sparse matrix by a scalar"""

    nrow = X.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(scalar)

    return w * X


class _Sparse:

    """Useful sparse functions."""

    def __init__(self, X):

        self.X = X
        self.tmp = self.X.copy()

    def variance(self, axis=0):
        """Calculate variance across the specified axis of a sparse matrix"""
        return _sparse_variance(self.X, self.tmp, axis)

    def row_wise_multiplication(self, scalar):

        """Multiply each row of sparse matrix by a scalar"""

        return _sparse_row_wise_multiplication(self.X, scalar)

    def mean_center(self, column_means=False):
        """Mean-center columns of a sparse matrix"""
        if not column_means:
            column_means = self.X.mean(axis=0)
        return self.X - column_means

    def normalize_variance(self, col_stdev=False):
        """variance-normalize columns of a sparse matrix"""
        if not col_stdev:
            col_stdev = np.sqrt(_sparse_variance(self.X, self.tmp, axis=0))
        return _sparse_row_wise_multiplication(self.X.T, scala1 / col_stdev).T

    def z_score(self, gene_mean=False, gene_stdev=False):

        """"""

        if not gene_mean:
            gene_mean = self.X.mean(0)
        if not gene_stdev:
            gene_stdev = np.sqrt(_sparse_variance(self.X, self.tmp))
        return _sparse_row_wise_multiplication((self.X - gene_mean).T, 1 / gene_stdev).T