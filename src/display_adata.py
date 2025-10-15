def pp_adata(adata):
    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of genes: {adata.n_vars}")
    print()
    print('-- Display:')
    display(adata)
    print()
    print('-- Display adata.obs.head:')
    display(adata.obs.head())
    print()
    print('-- Display adata.var.head:')
    display(adata.var.head())
    print()
    print('-- Nbr values in obs cols:')
    for c in adata.obs.columns:
       print(f"There are {len(set(adata.obs[c]))} values available for field: {c}")
    print()
    return
