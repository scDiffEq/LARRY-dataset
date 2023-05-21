
def cell_cycle_genes(genes_added = []):
    std = ["Ube2c", "Hmgb2", "Hmgn2", "Tuba1b", "Ccnb1", "Tubb5", "Top2a", "Tubb4b"]
    cc_genes = std + genes_added
    return cc_genes