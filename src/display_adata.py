"""
Misc utils to display info about adata
"""

from IPython.display import display

import pandas as pd
import scanpy as sc
# noinspection PyUnresolvedReferences
import hdf5plugin


def reset_df_index(df: pd.DataFrame, restart: int = 1) -> pd.DataFrame:
    """Convenience function"""
    df = df.reset_index(drop=True)
    if restart != 0:
        df.index += restart
    return df


def pp_adata(adata):
    if isinstance(adata, str):
        adata = sc.read_h5ad(adata, backed="r")

    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells: {adata.n_obs:,d}")
    print(f"Number of genes: {adata.n_vars:,d}")
    print()

    print('-- Display:')
    display(adata)
    print()

    print('-- Display adata.obs.head:')
    display(adata.obs.head())
    print()

    if 'X_hvg' in adata.obsm:
        print("-- adata.obsm['X_hvg']")
        # noinspection PyUnresolvedReferences
        print(f"shape = {adata.obsm['X_hvg'].shape}")
        print()

    print('-- Display adata.var.head:')
    display(adata.var.head())
    print()

    print('-- Nbr values in obs cols:')
    df = pd.DataFrame.from_records([(c, len(set(adata.obs[c]))) for c in adata.obs.columns],
                                   columns=["obs.col", "nUniq-vals"])
    print(reset_df_index(df).to_markdown(intfmt=",", floatfmt=",.0f"))
    print()

    if adata.isbacked:
        adata.file.close()

    return


def pp_train_test_splits(train_file, test_file):

    print("Reading data ...")
    adata_train = sc.read_h5ad(train_file)
    adata_test = sc.read_h5ad(test_file)
    print()

    print("##  Training Data info:")
    print()
    print("File:", train_file)
    print()
    pp_adata(adata_train)

    print("##  Test Data info:")
    print()
    print("File:", test_file)
    print()
    pp_adata(adata_test)

    print()
    print("##  Splits")
    print()

    trn_cnts = adata_train.obs['target_gene'].value_counts()
    trn_cnts.name = "train"

    tst_cnts = adata_test.obs['target_gene'].value_counts()
    tst_cnts.name = "test"

    df = pd.concat([trn_cnts, tst_cnts], axis=1)
    df['total_count'] = df.train + df.test
    df['train_pct'] = df.train / df.total_count
    df['test_pct'] = df.test / df.total_count

    print(df['total_count  train_pct  test_pct'.split()].to_markdown(floatfmt=['', ',.0f', '.2f', '.2f']))

    return adata_train, adata_test
