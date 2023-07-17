import requests
import anndata
import pathlib
import os

from ..utils import InfoMessage

class FigshareDownloader:
    """Download AnnData stored in Figshare"""

    def __init__(
        self, block_size: int = 1024, notebook: bool = True, silent: bool 
= False
    ):
        """
        # 1 Kibibyte
        """
        self._block_size = block_size
        self._notebook = notebook
        self._INFO = InfoMessage()

    @property
    def _RESPONSE(self):
        if not hasattr(self, "_response"):
            self._response = requests.get(self._url, stream=True)
            self._response.raise_for_status()
        return self._response

    @property
    def _TOTAL(self):
        """"""
        return int(self._RESPONSE.headers.get("content-length", 0))

    @property
    def _TQDM(self):
        if not hasattr(self, "_tqdm"):
            if self._notebook:
                import tqdm.notebook

                tqdm = tqdm.notebook.tqdm
            else:
                from tqdm import tqdm

            self._tqdm = tqdm
        return self._tqdm

    @property
    def _PROGRESS_BAR(self):
        if not hasattr(self, "_progress_bar"):
            self._progress_bar = self._TQDM(
                total=self._TOTAL, unit="iB", unit_scale=True
            )
        return self._progress_bar

    @property
    def _PASSES_DOWNLOAD_QC_CHECK(self):
        return self._PROGRESS_BAR.n == self._TOTAL

    @property
    def fpath(self):
        return pathlib.Path(self._fpath).absolute()

    @property
    def _FILE_DIR(self):
        return self.fpath.parent

    def _configure_file_destination(self):
        if not self._FILE_DIR.exists():
            self._FILE_DIR.mkdir()

    @property
    def _FILE_DOWNLOADED(self):
        return os.path.exists(self._fpath)

    def _read_adata(self):
        if not hasattr(self, "_adata"):
            self._adata = anndata.read_h5ad(self._fpath)
            self._adata.uns["source_url"] = self._url
        return self._adata

    def __call__(self, url, fpath, return_adata=True):
        """


        Parameters
        ----------

        Returns
        -------
        """

        self._url = url
        self._fpath = fpath

        self._configure_file_destination()

        if not self._FILE_DOWNLOADED:
            self._INFO(f"Downloading to: {fpath}")
            with open(self.fpath, "wb") as file:
                for data in self._RESPONSE.iter_content(self._block_size):
                    self._PROGRESS_BAR.update(len(data))
                    file.write(data)
            self._PROGRESS_BAR.close()
            self._PASSES_DOWNLOAD_QC_CHECK

        if return_adata:
            return self._read_adata()
        return self.fpath
