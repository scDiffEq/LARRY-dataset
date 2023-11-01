import pathlib
import pandas as pd


class _F_obs:
    def __init__(self):
        """"""

    @property
    def _PATH(self):
        return pathlib.Path(__file__).parent.joinpath("F_obs.csv")

    @property
    def df(self):

        df = pd.read_csv(self._PATH, index_col=0)
        df.index = df.index.astype(str)

        return df

    def __call__(self):
        return self.df
