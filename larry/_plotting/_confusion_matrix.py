

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sklearn
import os


from .. import utils

from ._multi_fate_in_vitro_color_map import multi_fate_in_vitro_cmap

NoneType = type(None)


class ConfusionMatrix(utils.ABCParse):
    def __init__(self):

        self.__parse__(locals())
        
        self.multi_fate_in_vitro_cmap = multi_fate_in_vitro_cmap()

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
    
    def _format_max_val(self):
        max_val = self.conf_mtx.max()
        
        if max_val >= 100:
            divisor, rounder = 100, -2
        else:
            divisor, rounder = 10, -1
        
        diff = divisor - (max_val % divisor)
        return np.round(max_val + diff, rounder)

    @property
    def _VMAX(self):
        if isinstance(self._vmax, NoneType):
            return self._format_max_val()
        return self._vmax
    
    @property
    def _CMAP(self):
        cmap_labels = []
        _df = self.multi_fate_in_vitro_cmap.copy()
        return [
            _df.loc[_df["label"] == label]["color"].values[0] for label in self.labels
        ]

    def __heatmap__(self):
        
        cg = sns.clustermap(
            self.conf_mtx,
            figsize=self._figsize,
            cmap=self._cmap,
            annot=self._ANNOT,
            xticklabels=self.labels,
            yticklabels=self.labels,
            row_colors=self._CMAP,
            col_colors=self._CMAP,
            row_cluster=False,
            col_cluster=False,
            fmt="",
            vmin=0,
            vmax=self._VMAX,
            cbar_pos=(0.92, 0.30, 0.02, 0.4),
            annot_kws={"size": 6},
        )
        cg.ax_heatmap.tick_params(axis="both", which="both", labelsize=6)
        cg.ax_heatmap.set_title(
            self._title,
            fontsize = self._title_fontsize,
            y = 1.04,
        )

    def __call__(
        self,
        F_obs,
        F_hat,
        labels=None,
        title=None,
        title_fontsize = 10,
        ax=None,
        cmap="Blues",
        vmax=None,
        figsize=(4, 4),
        save=False,
        *args,
        **kwargs,
    ):

        self.__parse__(
            locals(),
            public = ["F_obs", "F_hat"],
        )
        self.__heatmap__()

        if save:
            plt.savefig(save)
            

def confusion_matrix(
    F_obs,
    F_hat,
    labels=None,
    title=None,
    ax=None,
    cmap="Blues",
    vmax=None,
    figsize=(4, 4),
    save=False,
    *args,
    **kwargs,
):
    confusion_mtx = ConfusionMatrix()
    KWARGS = utils.extract_func_kwargs(confusion_mtx, locals())
    return confusion_mtx(**KWARGS)