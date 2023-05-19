

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sklearn
import os


from ._abc_parse import ABCParse
from ._multi_fate_color_map import mk_multifate_cmap

NoneType = type(None)

class ConfusionMatrix(ABCParse):
    def __init__(self):

        self.__parse__(locals())
        if not os.path.exists("label_cmap.csv"):
            self.label_cmap = mk_multifate_cmap()
        else:
            self.label_cmap = pd.read_csv("label_cmap.csv")

    @property
    def labels(self):
        if isinstance(self._labels, NoneType):
            return self.F_obs.sum(0).sort_values(ascending=False).index.tolist()
        return self._labels

    @property
    def y_true(self):
        if isinstance(self.F_obs, pd.DataFrame):
            return self.F_obs.idxmax(1).values
        return self.F_obs

    @property
    def y_pred(self):
        if isinstance(self.F_hat, pd.DataFrame):
            return self.F_hat.idxmax(1).values
        return self.F_hat

    @property
    def conf_mtx(self):
        return sklearn.metrics.confusion_matrix(
            y_true=self.y_true,
            y_pred=self.y_pred,
            labels=self.labels,
        )

    @property
    def _ANNOT(self):
        _ANNOT = self.conf_mtx.astype(str)
        _ANNOT[np.where(_ANNOT.astype(str) == "0")] = ""
        return _ANNOT

    @property
    def _VMAX(self):
        if isinstance(self._vmax, NoneType):
            return np.round(self.conf_mtx.max(), -2)
        return self._vmax

    def __heatmap__(self):
        
        cg = sns.clustermap(
            self.conf_mtx,
            figsize=self._figsize,
            cmap=self._cmap,
            annot=self._ANNOT,
            xticklabels=self.labels,
            yticklabels=self.labels,
            # need to fix here.... not sure why this doesn't work....?
            # needs re-ordering of some sort....
            row_colors=self.label_cmap["color"].to_numpy(),
            col_colors=self.label_cmap["color"].to_numpy(),
            row_cluster=False,
            col_cluster=False,
            fmt="",
            vmin=0,
            vmax=self._VMAX,
            cbar_pos=(0.92, 0.30, 0.02, 0.4),
            annot_kws={"size": 6},
        )
        cg.ax_heatmap.tick_params(axis="both", which="both", labelsize=6)

    def __call__(
        self,
        F_obs,
        F_hat,
        labels=None,
        title=None,
        ax=None,
        cmap="Blues",
        vmax=None,
        figsize=(4, 4),
        save=False,
    ):

        self.__parse__(
            locals(),
            private=["labels", "ax", "cmap", "figsize", "vmax"],
        )
        self.__heatmap__()

        if save:
            plt.savefig(save)

