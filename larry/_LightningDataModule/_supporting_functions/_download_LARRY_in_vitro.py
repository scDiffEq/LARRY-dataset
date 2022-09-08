
# import packages #
# --------------- #
import subprocess
import pydk
import os

## TO-DO: adjust this to be a general cloud download shortcut if accessible to the user

def _download_LARRY_in_vitro(
    gsutil_uri="scdiffeq-data/Weinreb2020/adata.fate_scores_added.h5ad",
    dest="./scdiffeq_data",
    force=False,
):
    h5ad_path = os.path.join(dest, os.path.basename(gsutil_uri))

    if not os.path.exists(h5ad_path) or force:
        gsutil_command = "gsutil -m cp -r gs://{} {}".format(gsutil_uri, dest).split(
            " "
        )
        pydk.mkdir_flex(dest)
        _ = subprocess.run(gsutil_command)

    return h5ad_path