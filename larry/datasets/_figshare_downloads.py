
from ._figshare_downloader import FigshareDownloader

def LARRY_task_one(
    
fpath="scdiffeq_data/adata.Weinreb2020.in_vitro.task_01.timepoint_recovery.h5ad",
    url="https://figshare.com/ndownloader/files/41576409",
):
    """Downloader for Timepoint Recovery Task AnnData"""

    figshare_downloader = FigshareDownloader()
    return figshare_downloader(url=url, fpath=fpath)
