

# -----------------------------------------------------------------------------
import numpy as np
from licorice_font import font_format


# -----------------------------------------------------------------------------
from ..._preprocess._annotate_test_train import _annotate_test_train


# -----------------------------------------------------------------------------
def _options():

    opt1 = font_format("fate_prediction", ["BOLD", "BLUE"])
    opt2 = font_format("timepoint_recovery", ["BOLD", "BLUE"])
    print("Choose: {} or {}".format(opt1, opt2))
    

# -----------------------------------------------------------------------------
class AnnDataTaskAnnotation:
    """Annotate adata according to ask with test/train split."""

    def __init__(self, adata, train_key="train", test_key="test", silent=True):

        self.adata = adata
        self._train_key = train_key
        self._test_key = test_key
        self._silent = silent

    def timepoint_recovery(
        self, train_time=[2, 6], test_time=[2, 4], time_key="Time point"
    ):

        """
        Notes:
        ------
        (1) We could make the these data vectors of type, pd.Categorical(), which
            is the standard AnnData convention. However, by not doing so, this
            allows us to slice the adata object directly.
        """

        adata = self.adata.copy()

        test, train = np.zeros(len(adata)), np.zeros(len(adata))
        test_idx = np.where(adata.obs[time_key].isin(test_time))[0]
        train_idx = np.where(adata.obs[time_key].isin(train_time))[0]
        test[test_idx], train[train_idx] = 1, 1

        adata.obs[self._train_key] = train.astype(bool)  # **1
        adata.obs[self._test_key] = test.astype(bool)  # **1

        return adata

    def fate_prediction(self):

        adata = self.adata.copy()
        _annotate_fate_test_train(adata)
        return adata

    def run(self, task, train_time=[2, 6], test_time=[2, 4], time_key="Time point"):

        if task == "fate_prediction":
            adata = self.fate_prediction()

        elif task == "timepoint_recovery":
            adata = self.timepoint_recovery(train_time, test_time, time_key)

        elif not task in ["fate_prediction", "timepoint_recovery"]:
            _options()

        adata.obs.index = adata.obs.index.astype(str)

        return adata
    
    
def _annotate_task_groups(adata, task, train_key, test_key, train_time, test_time, silent):
    
    task_test_train = AnnDataTaskAnnotation(
            adata=adata,
            train_key=train_key,
            test_key=test_key,
            silent=silent,
        )
    return task_test_train.run(
                task=task,
                train_time=train_time,
                test_time=test_time,
                time_key=time_key,
            )
