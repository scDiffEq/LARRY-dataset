
import ABCParse
import pickle
import pathlib
import pandas as pd
from typing import Dict


class InVitroFatePredictionSubsets(ABCParse.ABCParse):
    """
    Dictionary of train/test split subsets (indices and tables) for the
    LARRY in vitro fate prediction task.
    """

    def __init__(self, disable: bool = False, *args, **kwargs):
        
        file = pathlib.Path(__file__)
        fname = "FatePredictionSubsets.pkl"
        path = file.parent.joinpath(fname)
        
        self.__parse__(locals(), public = [None])
        
        if not self._disable:
            self._f = self._load()

    def _load(self) -> Dict:
        try:
            return pickle.load(self._path.open("rb"))
        except:
            return pd.read_pickle(self._path)

    @property
    def indices(self) -> Dict:
        return self._f["Indices"]

    @property
    def tables(self) -> Dict:
        return self._f["Tables"]