
# -- import local dependencies: ------------------------------------------------
from ._figshare_downloader import FigshareDownloader


# -- import packages: ----------------------------------------------------------
import anndata
import pathlib


# -- set typing: ---------------------------------------------------------------
from typing import Optional, Union


# -- API-facing function: ------------------------------------------------------
def LARRY_task_one(
    fpath: Optional[Union[str, pathlib.Path]] = None,
    url: Optional[Union[str, pathlib.Path]] = None,
) -> anndata.AnnData:
    """Downloader for Timepoint Recovery Task AnnData"""

    if fpath is None:
        fpath = (
            "scdiffeq_data/adata.Weinreb2020.in_vitro.task_01.timepoint_recovery.h5ad"
        )

    if url is None:
        url = "https://figshare.com/ndownloader/files/41576409"

    figshare_downloader = FigshareDownloader()
    return figshare_downloader(url=url, fpath=fpath)