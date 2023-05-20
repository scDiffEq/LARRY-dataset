import os
from pathlib import Path
from typing import Dict
import wget
import licorice_font

from .._utils import AutoParseBase
from ._directory_manager import mkdir


# -- base classes: -------------------------------------------------------------
class FormattedURLFile(AutoParseBase):
    def __init__(self, file, download_dir, download_structure, base_url):

        self.__parse__(locals(), public=[None])
    
    @property
    def _fname(self):
        return os.path.basename(self._file)

    @property
    def URL(self):
        return os.path.join(self._base_url, self._fname)

    @property
    def _download_path(self):
        return Path(os.path.join(self._download_dir, self._download_structure)).absolute()

    @property
    def fpath(self):
        return Path(os.path.join(self._download_path, self._fname)).absolute()

    @property
    def downloaded(self):
        return self.fpath.exists()
        

class URLPathInterface:
    
    def __init__(self, download_path=os.getcwd()):
        
        self._download_path = download_path
        self._base_url = "https://kleintools.hms.harvard.edu/paper_websites/state_fate2020/"

    def _format(self, file):
        formatted = FormattedURLFile(
            file=f"stateFate_{self._dataset}_{file}",
            download_dir=self._download_path,
            download_structure=self._download_structure,
            base_url=self._base_url,
            )
        self._download_dir = formatted._download_path
        return formatted

    @property
    def gene_names(self):
        return self._format(file = "gene_names.txt.gz")

    @property
    def clone_matrix(self):
        return self._format(file = "clone_matrix.mtx.gz")

    @property
    def normed_counts(self):
        return self._format(file = "normed_counts.mtx.gz")

    @property
    def metadata(self):
        return self._format(file = "metadata.txt.gz")

    @property
    def __attributes__(self)->Dict:

        AttrDict = {}
        for attr in self.__dir__():
            if not attr.startswith("_") and attr != "download":
                AttrDict[attr] = getattr(self, attr)

        return AttrDict

    def download(self, silent = False, bar = False):

        for fname, file in self.__attributes__.items():
            if not file.downloaded:
                path = str(file.fpath)
                mkdir(path = path, silent = silent)
                wget.download(url=file.URL, out= path, bar=bar)

    @property
    def _attrs(self):
        return {
            attr: getattr(self, attr)
            for attr in self.__dir__()
            if not attr.startswith("_")
        }

    @property
    def _header(self):
        header_items = [
            "{:<13}".format("Component"),
            "{:^10}".format("Downloaded"),
            "Filepath",
        ]
        download_header = (
            " | ".join(
                [licorice_font.font_format(item, ["BOLD"]) for item in header_items]
            )
            + "\n"
            + "".join(["-"] * 38)
            + "\n"
        )
        return download_header

    @property
    def _description(self):
        header = self._header
        for k, v in self._attrs.items():
            if not k == "download":
                if v.downloaded:
                    color = "GREEN"
                else:
                    color = "RED"
                header += "{:<13} |{:^25}| {}\n".format(
                    k,
                    licorice_font.font_format(str(v.downloaded), ["BOLD", color]),
                    v.fpath,
                )
        return header
    
    def __repr__(self):
        return self._description


# -- API-facing sub-classes: -------------------------------------------------
class inVitroURLPaths(URLPathInterface):
    _download_structure = "KleinLabData/in_vitro"
    _dataset = "inVitro"
    def __init__(self, download_path=os.getcwd()):
        super(inVitroURLPaths, self).__init__(download_path=download_path)


class inVivoURLPaths(URLPathInterface):
    _download_structure = "KleinLabData/in_vivo"
    _dataset = "inVivo"
    def __init__(self, download_path=os.getcwd()):
        super(inVivoURLPaths, self).__init__(download_path=download_path)


class CytokinePerturbationURLPaths(URLPathInterface):
    _download_structure = "KleinLabData/cytokine_perturbation"
    _dataset = "cytokinePerturbation"
    def __init__(self, download_path=os.getcwd()):
        super(CytokinePerturbationURLPaths, self).__init__(download_path=download_path)
        