
# -- import local dependencies: ------------------------------------------------
from ._figshare_downloader import FigshareDownloader


# -- import packages: ----------------------------------------------------------
import anndata
import pathlib


# -- set typing: ---------------------------------------------------------------
from typing import Optional, Union


# -- API-facing function: ------------------------------------------------------
def LARRY_timepoint_interpoltion(
    fpath: Optional[Union[str, pathlib.Path]] = None,
    url: Optional[Union[str, pathlib.Path]] = None,
) -> anndata.AnnData:
    """
    Downloader for Timepoint Recovery Task AnnData
    https://figshare.com/articles/dataset/LARRY_LT-scRNA-seq_dataset_formatted_for_timepoint_recovery/23692095
    """

    if fpath is None:
        fpath = (
            "scdiffeq_data/adata.Weinreb2020.in_vitro.task_01.timepoint_recovery.h5ad"
        )

    if url is None:
        url = "https://figshare.com/ndownloader/files/41576409"

    figshare_downloader = FigshareDownloader()
    return figshare_downloader(url=url, fpath=fpath)


def LARRY_fate_prediction(
    fpath: Optional[Union[str, pathlib.Path]] = None,
    url: Optional[Union[str, pathlib.Path]] = None,
) -> anndata.AnnData:
    """
    Downloader for Fate Prediction Task AnnData
    https://figshare.com/articles/dataset/LARRY_LT-scRNA-seq_dataset_formatted_for_fate_prediction/23699991
    """

    if fpath is None:
        fpath = (
            "scdiffeq_data/adata.Weinreb2020.in_vitro.task_02.fate_prediction.h5ad"
        )

    if url is None:
        url = "https://figshare.com/ndownloader/files/41591139"

    figshare_downloader = FigshareDownloader()
    return figshare_downloader(url=url, fpath=fpath)


def LARRY_spliced_unspliced_reads(
    fpath: Optional[Union[str, pathlib.Path]] = None,
    url: Optional[Union[str, pathlib.Path]] = None,
) -> anndata.AnnData:
    """
    Downloader for AnnData object containing spliced / unspliced reads.
    """

    if fpath is None:
        fpath = (
            "scdiffeq_data/adata.Weinreb2020.in_vitro.spliced_unspliced_reads.h5ad"
        )

    if url is None:
        url = "https://figshare.com/ndownloader/files/37028569"

    figshare_downloader = FigshareDownloader()
    return figshare_downloader(url=url, fpath=fpath)