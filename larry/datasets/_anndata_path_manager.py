import os
from pathlib import Path
import licorice_font


class AnnDataPathManager:
    def __init__(self, dataset, download_dir):

        self._dataset = dataset
        self._download_dir = os.path.join(download_dir, f"KleinLabData/{self._dataset}")

    def _path_constructor(self, obj_specifier):
        fpath = f"adata.Weinreb2020.{self._dataset}.{obj_specifier}.h5ad"
        return Path(os.path.join(self._download_dir, fpath))

    @property
    def raw(self):
        return self._path_constructor("raw")

    @property
    def gene_filtered(self):
        return self._path_constructor("gene_filtered")

    @property
    def timepoint_recovery(self):
        return self._path_constructor("task_01.timepoint_recovery")

    @property
    def fate_prediction(self):
        return self._path_constructor("task_02.fate_prediction")

    @property
    def uniform(self):
        return self._path_constructor("uniform")

    def __repr__(self):

        description = "{}:\n".format(
            licorice_font.font_format("AnnData Objects", ["BOLD"])
        )
        description += "----------------"

        attrs = {
            attr: getattr(self, attr)
            for attr in self.__dir__()
            if not attr.startswith("_")
        }

        for k, v in attrs.items():
            if v.exists():
                exists = "[{:^20}]".format(
                    licorice_font.font_format(str(v.exists()), ["BOLD", "BLUE"])
                )
            else:
                exists = "[{:^15}]".format(
                    licorice_font.font_format(str(v.exists()), ["BOLD"])
                )
            description += f"\n- {exists} | adata.Weinreb2020.{self._dataset}.{k}.h5ad"

        return description