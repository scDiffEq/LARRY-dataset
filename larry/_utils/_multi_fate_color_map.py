
import pandas as pd

MultiFateColorMap = {}
MultiFateColorMap["Baso"] = {
    "Baso": "#0a9396",
    "Baso_Eos": "#29a1a4",
    "Baso_Eos_Neutrophil": "#48aeb1",
    "Baso_Erythroid": "#67bcbe",
    "Baso_Mast": "#85c9cb",
    "Baso_Meg": "#a4d7d8",
    "Baso_Monocyte": "#c2e4e5",
    "Baso_Neutrophil": "#e1f2f2",
}
MultiFateColorMap["Ccr_7"] = {
    "Ccr7_DC_Lymphoid": "#ae2012",
}
MultiFateColorMap["Eos"] = {
    "Eos": "#005F73",
    "Eos_Baso": "#80AFB9",
    "Eos_Neutrophil": "#C0D7DC",
}
MultiFateColorMap["Erythroid"] = {
    "Erythroid": "#E9D8A6",
    "Erythroid_Mast": "#EFE2BD",
    "Erythroid_Meg": "#FAF6E9",
}
MultiFateColorMap["Lymphoid"] = {
    "Lymphoid": "#EE9B00",
    "Lymphoid_Monocyte": "#F3B440",
    "Lymphoid_pDC": "#F7CD80",
}
MultiFateColorMap["Mast"] = {
    "Mast": "#94D2BD",
    "Mast_Baso": "#AFDECE",
    "Mast_Erythroid": "#BDE4D6",
    "Mast_Meg": "#CAE9DE",
    "Mast_Neutrophil": "#E5F4EF",
}
MultiFateColorMap["Meg"] = {
    "Meg": "#9B2226",
    "Meg_Baso": "#A83E42",
    "Meg_Erythroid": "#B45A5D",
    "Meg_Monocyte": "#CD9193",
    "Meg_Neutrophil": "#E6C8C9",
}
MultiFateColorMap["Monocyte"] = {
    "Monocyte": "#F08700",
    "Monocyte_Lymphoid": "#F29620",
    "Monocyte_Mast": "#F4A540",
    "Monocyte_Neutrophil": "#F8C380",
}
MultiFateColorMap["Neutrophil"] = {
    "Neutrophil": "#023047",
    "Neutrophil_Baso": "#224A5E",
    "Neutrophil_Ccr7_DC": "#426475",
    "Neutrophil_Monocyte": "#8198A3",
}

def mk_multifate_cmap():
    
    # -- make a uniform dict: --------------
    multifate_cmap = {}
    for val in MultiFateColorMap.values():
        multifate_cmap.update(val)
    multifate_cmap["None"] = "lightgrey"
    # --------------------------------------

    label_cmap = (
        pd.DataFrame.from_dict(multifate_cmap, orient="index")
        .reset_index()
        .rename({"index": "label", 0: "color"}, axis=1)
    )

    label_cmap.to_csv("label_cmap.csv")
    
    return label_cmap