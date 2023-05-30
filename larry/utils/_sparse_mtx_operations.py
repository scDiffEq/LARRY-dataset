
import scipy
import numpy as np

def sparse_var(E, axis=0):
    """calculate variance across the specified axis of a sparse matrix."""

    mean_gene = E.mean(axis=axis).A.squeeze()
    tmp = E.copy()
    tmp.data **= 2
    return tmp.mean(axis=axis).A.squeeze() - mean_gene ** 2

def sparse_rowwise_multiply(E, a):
    ''' multiply each row of sparse matrix by a scalar '''

    nrow = E.shape[0]
    w = scipy.sparse.lil_matrix((nrow, nrow))
    w.setdiag(a)
    return w * E

def mean_center(E, column_means=None):
    """mean-center columns of a sparse matrix"""

    if column_means is None:
        column_means = E.mean(axis=0)
    return E - column_means

def normalize_variance(E, column_stdevs=None):
    ''' variance-normalize columns of a sparse matrix '''

    if column_stdevs is None:
        column_stdevs = np.sqrt(sparse_var(E, axis=0))
    return sparse_rowwise_multiply(E.T, 1 / column_stdevs).T

def sparse_zscore(E, gene_mean=None, gene_stdev=None):
    """z-score normalize each column of a sparse matrix"""
    if gene_mean is None:
        gene_mean = E.mean(0)
    if gene_stdev is None:
        gene_stdev = np.sqrt(sparse_var(E))
    return sparse_rowwise_multiply((E - gene_mean).T, 1/gene_stdev).T
