
import pydk
import numpy as np


# -----------------------------------------------------------------------------
def _get_boolean_on_idx(full_length, idx):

    """Starting from an array of zeros, set a boolean value True (i.e., =1) if it's in the idx passed"""

    _tmp = np.zeros(full_length)
    _tmp[idx] = 1
    
    return _tmp.astype(bool)


def _annotate_train_test(adata, key="Well", train=[0, 1], test=[0, 2]):
    
    """Annotate test/train spit using a column of adata.obs[key]"""

    df = adata.obs.copy()

    train_idx = df.loc[df[key].isin(train)].index.astype(int)
    test_idx  = df.loc[df[key].isin(test)].index.astype(int)

    df["train"] = _get_boolean_on_idx(len(df), train_idx)
    df["test"]  = _get_boolean_on_idx(len(df), test_idx)
    
    df['cell_idx'] = df['cell_idx'].astype(int)
    df = df.sort_values("cell_idx").reset_index(drop=True)
    
    adata.obs = df

def _annotate_unique_test_train_lineages(adata, lineage_key="clone_idx", silent=True):

    """
    Annotate cells that belong to lineages only in the training set and only in the test set.
    
    Parameters:
    -----------
    adata
    
    lineage_key
        default: "clone_idx"
        type: str
    
    Returns:
    --------
    None, adata modified in-place.
    
    Notes:
    ------
    (1) Everything is non-unique in the test/train when you include well 0 because
        this is our common ground. Here, we start from wells 1 and 2.
    """

    df = adata.obs.copy()

    well_01 = df.loc[df["Well"] == 1].dropna()[lineage_key].unique().astype(float)
    well_02 = df.loc[df["Well"] == 2].dropna()[lineage_key].unique().astype(float)
    
    well_unique = pydk.unique(well_01, well_02, ["train", "test"])

    n_train = well_unique["train"].shape[0]
    n_test  = well_unique["test"].shape[0]

    if not silent:
        print("Unique training lineages: {}".format(n_train))
        print("Unique test lineages: {}".format(n_test))

    train_idx = df.loc[df[lineage_key].isin(well_unique["train"])].index.astype(int)
    test_idx  = df.loc[df[lineage_key].isin(well_unique["test"])].index.astype(int)
    
    adata.obs["unique_train_lineage"] = _get_boolean_on_idx(len(adata), train_idx)
    adata.obs["unique_test_lineage"]  = _get_boolean_on_idx(len(adata), test_idx)
    
    
# -----------------------------------------------------------------------------
def _annotate_fate_test_train(adata):
    
    _annotate_train_test(adata, key="Well", train=[0, 1], test=[0, 2])
    _annotate_unique_test_train_lineages(adata, lineage_key="clone_idx", silent=True)