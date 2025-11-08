"""
Misc utils to display info about adata
"""

from IPython.display import display

import pandas as pd
import scanpy as sc
# noinspection PyUnresolvedReferences
import hdf5plugin


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


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

    # ---
    def pp_obs_fld(fld):
        if fld in adata.obs:
            print(f"-- Values in adata.obs['{fld}']:")
            uniq = adata.obs[fld].unique()
            print(f"    nbr Unique values = {len(uniq):,d}")
            print(f"    first few values = {uniq[:10].tolist()}")
            print()
        return
    # ---

    print('-- Display adata.obs.head:')
    display(adata.obs.head())
    print()

    pp_obs_fld('tissue_ontology_term_id')
    pp_obs_fld('perturbation_name')
    pp_obs_fld('target_gene')

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


# ======================================================================================================
#   Main
# ======================================================================================================

# To run
# ------
#
# [Python]$ python -m test_loader loaders tmp/split2.toml
#
#

if __name__ == "__main__":

    import argparse
    from datetime import datetime

    _argparser = argparse.ArgumentParser(
        description='Display AnnData objects.',
    )

    _subparsers = _argparser.add_subparsers(dest='subcmd',
                                            title='Available commands',
                                            )
    # Make the sub-commands required
    _subparsers.required = True

    # ... pp ADATA_PATH
    _sub_cmd_parser = _subparsers.add_parser('pp',
                                             help="DIspay summary of AnnData.")
    _sub_cmd_parser.add_argument('adata_path', type=str,
                                 help="Path to AnnData h5ad file.")

    # ... loaders TOML_PATH
    _sub_cmd_parser = _subparsers.add_parser('loaders',
                                             help="Test data loaders for TOML.")
    _sub_cmd_parser.add_argument("-l", '--log_info', action="store_true",
                                 help="Log INFO messages.")
    _sub_cmd_parser.add_argument('toml', type=str,
                                 help="Path to TOML file.")

    # ...

    _args = _argparser.parse_args()
    # .................................................................................................

    start_time_ = datetime.now()

    print("---------------------------------------------------------------------")

    if _args.subcmd == 'pp':

        pp_adata(_args.adata_path)

    else:

        raise NotImplementedError(f"Command not implemented: {_args.subcmd}")

    # /

    print('\nTotal Run time =', datetime.now() - start_time_)
    print()
