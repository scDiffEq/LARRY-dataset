
from licorice_font import font_format

    
from ._time_occupance import time_occupance


def fated_idx(
    adata,
    fate_time=[4, 6],
    exclude_fate: tuple = ("Cell type annotation", ["undiff"]),
    lineage_key: str = "clone_idx",
    time_key: str = "Time point",
) -> list:
    """
    generate a list of lineages that are seen at d2 as well as d4 and d6

    Notes:
    ------
    (1) This function returns the indices of lineages NOT cells.
    """
    
    if "time_occupance" in adata.uns_keys():
        t_occ = adata.uns["time_occupance"]
    else:
        t_occ = time_occupance(
            adata,
            fate_time=fate_time,
            exclude_fate=exclude_fate,
            lineage_key=lineage_key,
            time_key=time_key,
            return_df=True,
        )
    
    return t_occ[t_occ["2"] & t_occ[["4", "6"]].any(axis=1)].index.tolist()


def annotate_fated(
    adata,
    lineage_key="clone_idx",
    time_key="Time point",
    t0=2,
    fate_time=[4, 6],
    key_added="fate_observed",
    t0_key_added="t0_fated",
    exclude_fate: tuple = ("Cell type annotation", ["undiff"]),
) -> None:

    """We use this function to denote lineages with cells at d2
    and one or more cells in [d4, d6]
    
    Updates adata.obs with two columns -> adata.obs[['fate_observed', 't0_fated']]
    """

    df = adata.obs.copy()

    f_idx = fated_idx(
        adata,
        fate_time=fate_time,
        exclude_fate=exclude_fate,
        lineage_key=lineage_key,
        time_key=time_key,
    )
    
    info = font_format("INFO", ['PURPLE'])
    df[key_added] = df[lineage_key].isin(f_idx)
    msg = f"- [ {info} ] | Fated cells annotated at: adata.obs['{key_added}']"
    print(msg)
    df[t0_key_added] = (df[time_key] == t0) & df[key_added]
    msg = f"- [ {info} ] | Fated cells (t=t0) annotated at: adata.obs['{t0_key_added}']"
    print(msg)

    adata.obs = df
