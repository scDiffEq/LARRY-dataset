
import pandas as pd
import ABCParse

from ._in_vitro_color_map import InVitroColorMap
from .. import utils


class MultiFateInVitroColorMap(ABCParse.ABCParse):
    def __init__(self, save=False, fname="label_cmap.csv"):

        self.__parse__(locals(), public=[None])
        self._CMAP = InVitroColorMap()._dict
        self.__build_cmap__()
        if self._save:
            self._to_csv()

    @property
    def Baso(self):
        return {
            "Baso": "#0a9396",
            "Baso_Eos": "#29a1a4",
            "Baso_Eos_Neutrophil": "#48aeb1",
            "Baso_Erythroid": "#67bcbe",
            "Baso_Mast": "#85c9cb",
            "Baso_Meg": "#a4d7d8",
            "Baso_Monocyte": "#c2e4e5",
            "Baso_Neutrophil": "#e1f2f2",
        }

    @property
    def Ccr_7(self):
        return {
            "Ccr_7": "#ae2012",
            "Ccr7_DC_Lymphoid": "#ae2012",
        }

    @property
    def Eos(self):
        return {
            "Eos": "#005F73",
            "Eos_Baso": "#80AFB9",
            "Eos_Neutrophil": "#C0D7DC",
        }

    @property
    def Erythroid(self) -> dict:
        return {
            "Erythroid": "#E9D8A6",
            "Erythroid_Mast": "#EFE2BD",
            "Erythroid_Meg": "#FAF6E9",
        }

    @property
    def Lymphoid(self) -> dict:
        return {
            "Lymphoid": "#EE9B00",
            "Lymphoid_Monocyte": "#F3B440",
            "Lymphoid_pDC": "#F7CD80",
        }

    @property
    def Mast(self) -> dict:
        return {
            "Mast": "#94D2BD",
            "Mast_Baso": "#AFDECE",
            "Mast_Erythroid": "#BDE4D6",
            "Mast_Meg": "#CAE9DE",
            "Mast_Neutrophil": "#E5F4EF",
        }

    @property
    def Meg(self) -> dict:
        return {
            "Meg": "#9B2226",
            "Meg_Baso": "#A83E42",
            "Meg_Erythroid": "#B45A5D",
            "Meg_Monocyte": "#CD9193",
            "Meg_Neutrophil": "#E6C8C9",
        }

    @property
    def Monocyte(self) -> dict:
        return {
            "Monocyte": "#F08700",
            "Monocyte_Lymphoid": "#F29620",
            "Monocyte_Mast": "#F4A540",
            "Monocyte_Neutrophil": "#F8C380",
        }

    @property
    def Neutrophil(self) -> dict:
        return {
            "Neutrophil": "#023047",
            "Neutrophil_Baso": "#224A5E",
            "Neutrophil_Ccr7_DC": "#426475",
            "Neutrophil_Monocyte": "#8198A3",
        }

    @property
    def cmap(self):
        return (
            pd.DataFrame.from_dict(self._CMAP, orient="index")
            .reset_index()
            .rename({"index": "label", 0: "color"}, axis=1)
        )

    def __build_cmap__(self):
        for attr in self.__dir__():
            if not (attr.startswith("_")) and (attr != "cmap"):
                self._CMAP.update(getattr(self, attr))

    def _to_csv(self):
        self.cmap.to_csv(self._fname)
        
def multi_fate_in_vitro_cmap():
    return MultiFateInVitroColorMap().cmap