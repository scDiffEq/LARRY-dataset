
# import packages #
# --------------- #
from torch.utils.data import random_split


def _split_training_data(train_dataset, train_val_split):

    """

    Parameters:
    -----------
    train_dataset

    train_val_split

    Returns:
    --------
    train, val
    """

    n_train_init = train_dataset.__len__()
    n_train = int(n_train_init * train_val_split)
    n_val = int(n_train_init - n_train)

    return random_split(train_dataset, [n_train, n_val])