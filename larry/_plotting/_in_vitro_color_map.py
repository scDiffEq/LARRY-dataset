
import pandas as pd


class InVitroColorMap:
    def __init__(self):
        self._dict = {
            "Neutrophil": "#023047",
            "Eos": "#005f73",
            "Baso": "#0a9396",
            "Mast": "#94d2bd",
            "Erythroid": "#e9d8a6",
            "Lymphoid": "#ee9b00",
            "Monocyte": "#F08700",
            "pDC": "#bb3e03",
            "Ccr7_DC": "#ae2012",
            "Meg": "#9b2226",
            "Undifferentiated": "#f0efeb",
            "undiff": "#f0efeb",
        }

    @property
    def _KEYS(self):
        return list(self._dict.keys())

    @property
    def _VALS(self):
        return list(self._dict.values())

    @property
    def df(self):
        return pd.DataFrame({"Cell type annotation": self._KEYS, "color": self._VALS})
    
def in_vitro_cmap():
    return InVitroColorMap().df
