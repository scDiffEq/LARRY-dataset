

def _count_report_features(adata):

    adata.uns["n_hv"] = n_hv = adata.var["highly_variable"].sum()
    adata.uns["n_corr_cell_cycle"] = n_corr_cell_cycle = adata.var[
        "corr_cell_cycle"
    ].sum()
    adata.uns["n_mito"] = n_mito = adata.var["gene_name"].str.startswith("mt").sum()
    adata.uns["n_pass"] = n_pass = adata.var["pass_filter"].sum()
    adata.uns["n_total"] = n_total = adata.var.shape[0]

    messages = [
        "{:<5} genes annotated as correlated with cell cycle".format(n_corr_cell_cycle),
        "{:<5} genes annotated as mitochondrial".format(n_mito),
        "{:<5} genes annotated as highly-variable".format(n_hv),
        "\n{:<5} genes / {} total genes passed filters".format(n_pass, n_total),
    ]

    for msg in messages:
        print(msg)