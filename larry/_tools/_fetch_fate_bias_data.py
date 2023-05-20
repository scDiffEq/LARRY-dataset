
import os
import glob
import numpy as np
import pandas as pd


def fetch_fate_bias_data():
    
    data = {}
    
    dpath = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "fate_bias_files")
    
    data["NM_traj"] = np.load(os.path.join(dpath, "neutrophil_monocyte_trajectory_mask.npy"))
    data["early_cells"] = np.load(os.path.join(dpath, "early_cells.npy"))
    data["fate_names"] = pd.read_csv(os.path.join(dpath, "FINAL_fate_names.txt"), header=None)[0].values.tolist()
    data["clonal_fate_matrix"] = np.load(os.path.join(dpath, "clonal_fate_matrix.npy"))
    data["timepoints"] = np.load(os.path.join(dpath, "timepoints.npy"))
    data['fate_df'] = fate_df = pd.DataFrame(data["clonal_fate_matrix"], columns=data["fate_names"])
    
    NM_cols = ["Neutrophil", "Monocyte"]
    data['NM_fate_df'] = fate_df[fate_df[NM_cols].sum(1) > 0][NM_cols]
    
    return data
