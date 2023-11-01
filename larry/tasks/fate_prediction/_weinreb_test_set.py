
import ABCParse
import pathlib
import numpy as np
import pandas as pd
import anndata

from typing import Optional


class WeinrebTestSet(ABCParse.ABCParse):
    """335 cell test set used in Weinreb, 2020 and Yeo, 2021"""

    def __init__(self, *args, **kwargs):
        
        """
        Parameters
        ----------
        base_path: str
        
        """

        self.__parse__(locals())
        
    @property
    def _BASE_PATH(self):
        
        self._sub_path = "assets/files/Figure5_Weinreb2020"
        return list(pathlib.Path(__file__).parents)[3].joinpath(self._sub_path)

    @property
    def _TIME_POINTS(self):
        if not hasattr(self, "_time_points"):
            PATH = self._BASE_PATH.joinpath("timepoints.npy")
            self._time_points = np.load(PATH)
        return self._time_points

    @property
    def _NM_MASK(self):
        if not hasattr(self, "_nm_mask"):
            PATH = self._BASE_PATH.joinpath("neutrophil_monocyte_trajectory_mask.npy")
            self._nm_mask = np.load(PATH)
        return self._nm_mask

    @property
    def _CLONAL_FATE_MATRIX(self):
        if not hasattr(self, "_clonal_fate_matrix"):
            PATH = self._BASE_PATH.parent.parent.parent.joinpath(
                "fate_bias_files/clonal_fate_matrix.npy"
            )
            self._clonal_fate_matrix = np.load(PATH)
        return self._clonal_fate_matrix

    @property
    def _PBA_predictions(self):
        if not hasattr(self, "_pba_preds"):
            PATH = self._BASE_PATH.joinpath("PBA_predictions.npy")
            self._pba_preds = np.load(PATH)
        return self._pba_preds

    @property
    def _FateID_predictions(self):
        if not hasattr(self, "_fateid_preds"):
            PATH = self._BASE_PATH.joinpath("FateID_predictions.npy")
            self._fateid_preds = np.load(PATH)
        return self._fateid_preds

    @property
    def _WOT_predictions(self):
        if not hasattr(self, "_wot_preds"):
            PATH = self._BASE_PATH.joinpath("WOT_predictions.npy")
            self._wot_preds = np.load(PATH)
        return self._wot_preds

    @property
    def _EARLY_CELLS(self):
        if not hasattr(self, "_early_cells"):
            PATH = self._BASE_PATH.joinpath("early_cells.npy")
            self._early_cells = np.load(PATH)
        return self._early_cells

    @property
    def _HELDOUT_MASK(self):
        if not hasattr(self, "_heldout_mask"):
            PATH = self._BASE_PATH.joinpath("heldout_mask.npy")
            self._heldout_mask = np.load(PATH)
        return self._heldout_mask

    @property
    def _SMOOTHED_GROUNDTRUTH(self):
        if not hasattr(self, "_smoothed_groundtruth"):
            PATH = self._BASE_PATH.joinpath("smoothed_groundtruth_from_heldout.npy")
            self._smoothed_groundtruth = np.load(PATH)
        return self._smoothed_groundtruth

    @property
    def _SMOOTHED_GROUNDTRUTH_SUBSET(self):
        if not hasattr(self, "_smoothed_groundtruth_subset"):
            self._smoothed_groundtruth_subset = self._SMOOTHED_GROUNDTRUTH[
                np.all(
                    [
                        self._EARLY_CELLS[self._NM_MASK[self._TIME_POINTS == 2]],
                        ~self._HELDOUT_MASK,
                    ],
                    axis=0,
                )
            ]
        return self._smoothed_groundtruth_subset

    @property
    def _NM_COUNTS(self):
        return self._CLONAL_FATE_MATRIX[:, 5:7].sum(1)

    @property
    def _HAS_FATE_MASK(self):
        return np.all(
            [self._NM_COUNTS > 0, self._NM_MASK[self._TIME_POINTS == 2]], axis=0
        )

    @property
    def _NM_PERCENT(self):
        return self._CLONAL_FATE_MATRIX[
            self._HAS_FATE_MASK, 5
        ] / self._CLONAL_FATE_MATRIX[self._HAS_FATE_MASK, 5:7].sum(1)

    @property
    def _HELDOUT_MASK_EARLY_CELLS(self):
        return (~self._HELDOUT_MASK)[
            self._EARLY_CELLS[self._NM_MASK[self._TIME_POINTS == 2]]
        ]

    def _subset_preds(self, method_prediction):
        return method_prediction[
            self._EARLY_CELLS[self._NM_MASK[self._TIME_POINTS == 2]]
        ][self._HELDOUT_MASK_EARLY_CELLS]

    @property
    def _PREDICTIONS(self):
        if not hasattr(self, "_predictions"):
            predictions = {
                "smoothed_ground_truth": self._SMOOTHED_GROUNDTRUTH_SUBSET,
            }
            for attr_name in self.__dir__():
                if "_predictions" in attr_name:
                    preds = getattr(self, attr_name)
                    predictions[attr_name.strip("_").split("_")[0]] = self._subset_preds(preds)
            self._predictions = predictions
        return self._predictions

    @property
    def _CLONAL_DATA_NO_HELDOUT(self):
        return self._NM_PERCENT[
            np.all(
                [
                    self._EARLY_CELLS[self._NM_MASK[self._TIME_POINTS == 2]],
                    ~self._HELDOUT_MASK,
                ],
                axis=0,
            )[self._HAS_FATE_MASK[self._NM_MASK[self._TIME_POINTS == 2]]]
        ]

    @property
    def _HAS_FATE_MASK_EARLY_HELDOUT(self):
        return self._HAS_FATE_MASK[
            np.all([self._NM_MASK[self._TIME_POINTS == 2], self._EARLY_CELLS], axis=0)
        ][self._HELDOUT_MASK_EARLY_CELLS]

    @property
    def _PREDICTIONS_FILTERED(self):
        if not hasattr(self, "_predictions_filtered"):
            self._predictions_filtered = {
                method: pred[self._HAS_FATE_MASK_EARLY_HELDOUT]
                for method, pred in self._PREDICTIONS.items()
            }
        return self._predictions_filtered

    @property
    def _D2_IDX(self):
        if hasattr(self, "_adata"):
            return self._adata.obs.loc[self._adata.obs["Time point"] == 2].index

    @property
    def _TEST_SET_IDX(self):
        if not self._D2_IDX is None:
            return self._D2_IDX[
                np.all(
                    [
                        self._NM_MASK[self._TIME_POINTS == 2],
                        self._EARLY_CELLS,
                    ],
                    axis=0,
                )
            ][self._HELDOUT_MASK_EARLY_CELLS][self._HAS_FATE_MASK_EARLY_HELDOUT]

    @property
    def df(self):
        if not hasattr(self, "_df"):
            COLS = [
                'clonal',
                'smoothed_ground_truth',
                'PBA',
                'FateID',
                'WOT',
            ]
            pred_dict = self._PREDICTIONS_FILTERED
            pred_dict["clonal"] = self._CLONAL_DATA_NO_HELDOUT
            self._df = pd.DataFrame(pred_dict, index=self._TEST_SET_IDX)[COLS]
            
        return self._df

    def __call__(self, adata: Optional[anndata.AnnData] = None, *args, **kwargs):
        
        """
        Parameters
        ----------
        adata: Optional[anndata.AnnData], default = None
        
        Returns
        -------
        test_set_df: pd.DataFrame
        """

        self.__update__(locals())

        return self.df

    
def weinreb_test_set(adata: Optional[anndata.AnnData] = None) -> pd.DataFrame:
    """
    Parameters
    ----------
    adata: Optional[anndata.AnnData], default = None

    Returns
    -------
    test_set_df: pd.DataFrame
    """
    
    test_set = WeinrebTestSet()

    return test_set(adata = adata)
