"""
Misc utils to display info about adata
"""

from IPython.display import display

from collections import defaultdict
import os.path
from pathlib import Path
import pprint

import anndata as an
import pandas as pd
import scanpy as sc
# noinspection PyUnresolvedReferences
import hdf5plugin

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from cell_load.data_modules.perturbation_dataloader import PerturbationDataModule


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


def plot_distribution(vec: np.ndarray, title: str):
    """
    Plots the distribution of vector `vec`
    """

    sns.histplot(vec, kde=True, bins=30, color='darkblue', edgecolor='black')
    plt.title(f"Distribution of {title}")
    plt.xlabel("Value")
    plt.ylabel("Density/Count")
    plt.show()
    return


def pp_underlined_hdg(hdg, overline=False, linechar="-", file=None):
    if overline:
        print(linechar * len(hdg), file=file)
    print(hdg, linechar * len(hdg), sep="\n", file=file)
    print(file=file)
    return


def reset_df_index(df: pd.DataFrame, restart: int = 1) -> pd.DataFrame:
    """Convenience function"""
    df = df.reset_index(drop=True)
    if restart != 0:
        df.index += restart
    return df


def pp_adata(adata: an.AnnData | str,
             cell_type_key: str = "tissue_ontology_term_id",
             pert_col: str = "target_gene",
             control_pert: str = "non-targeting",
             batch_col: str = "batch",
             ):

    if isinstance(adata, str):
        adata = sc.read_h5ad(adata, backed="r")

    print(f"Dataset shape: {adata.shape}")
    print(f"Number of cells  (.n_obs): {adata.n_obs:,d}")
    print(f"Number of genes (.n_vars): {adata.n_vars:,d}")
    print()

    print('-- Display:')
    display(adata)
    print()

    # ---
    def pp_obs_fld(fld, fldname: str = None, chkvalue: str = None):
        # noinspection PyPackages
        if fld in adata.obs:
            if fldname is not None:
                fldname = f"{fldname}: "
            else:
                fldname = ""

            print(f"-- Values in {fldname}adata.obs['{fld}']:")

            value_counts = adata.obs[fld].value_counts()

            print(f"    nbr Unique values = {len(value_counts):,d}")
            print(f"    first few values:", end="")
            print("", *value_counts[:10].to_markdown(floatfmt=',.0f').splitlines(), sep="\n\t")
            if chkvalue is not None and chkvalue not in value_counts.index[:10]:
                print(f"\t{chkvalue} = {value_counts[chkvalue]:,d}")
            print()
        return

    def top_vals_str(col, maxlen=50):
        top_vals = col.value_counts().head().index.to_list()
        s = ", ".join([str(x) for x in top_vals])
        if len(s) > maxlen:
            s = s[:maxlen - 3] + "..."
        return s
    # ---

    print('-- Display adata.obs.head:')
    display(adata.obs.head())
    print()

    pp_obs_fld(cell_type_key, "cell_type_key")
    pp_obs_fld(pert_col, "pert_col", control_pert)
    pp_obs_fld('perturbation_name')
    pp_obs_fld(batch_col, "batch_col")

    print("-- adata.obsm['X_hvg']")
    if 'X_hvg' in adata.obsm:
        # noinspection PyUnresolvedReferences
        print(f"shape = {adata.obsm['X_hvg'].shape}")
    else:
        print("... not found.")
    print()

    print('-- Display adata.var.head: shape =', adata.var.shape)
    display(adata.var.head())
    print()

    print('-- Unique values in obs cols:')
    df = adata.obs.nunique().to_frame(name='nUniq-vals')
    df["Top-vals"] = adata.obs.apply(lambda col: top_vals_str(col))
    df.index.name = "obs.col"
    print(df.to_markdown(intfmt=",", floatfmt=",.0f"))
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


def describe_toml(toml_path: str,
                  cell_type_key: str = "cell_type",
                  pert_col: str = "gene",
                  control_pert: str = "non-targeting",
                  batch_col: str = "gem_group",
                  ):
    """
    Prints summary of AnnData objects as described in TOML file.
    """

    print()
    pp_underlined_hdg("TOML Description", linechar="=")

    print(f"toml_path = '{toml_path}'")
    print(f"cell_type_key = '{cell_type_key}'")
    print(f"pert_col = '{pert_col}'")
    print(f"control_pert = '{control_pert}'")
    print(f"batch_col = '{batch_col}'")

    print("\n")

    # Arbitrary values, won't be used
    batch_size = 4
    cell_sentence_len = 16

    pdm = PerturbationDataModule(toml_path,
                                 batch_size=batch_size,
                                 cell_sentence_len=cell_sentence_len,
                                 cell_type_key=cell_type_key,
                                 pert_col=pert_col,
                                 batch_col=batch_col,
                                 control_pert=control_pert,
                                 embed_key=None,
                                 output_space="gene",
                                 basal_mapping_strategy="random",
                                 num_workers=1,
                                 cache_perturbation_control_pairs=False,
                                 barcode=False,
                                 )
    pdm.setup()
    config = pdm.config

    print("\nPerturbationDataModule created.\n")

    for dataset in config.get_all_datasets():
        dataset_path = config.datasets[dataset]

        # noinspection PyProtectedMember
        files = pdm._find_dataset_files(Path(dataset_path))

        pp_underlined_hdg(f"Dataset {dataset} ... {os.path.basename(dataset_path)}", overline=True)

        print("Fewshot perturbations:\n    ", end="")
        if fewshot_perts := config.get_fewshot_celltypes(dataset):
            print("\n    ", end="")
            pprint.pp(fewshot_perts)
        else:
            print("None.")
        print()

        print("Zeroshot perturbations:\n    ", end="")
        if zeroshot_perts := config.get_zeroshot_celltypes(dataset):
            print("\n    ", end="")
            pprint.pp(zeroshot_perts)
        else:
            print("None.")
        print()

        for fname, fpath in files.items():

            print()
            pp_underlined_hdg(f"AnnData file {fname} ... {os.path.basename(fpath)}")

            adata = sc.read_h5ad(fpath, backed="r")

            pp_adata(adata, cell_type_key=cell_type_key, pert_col=pert_col, control_pert=control_pert,
                     batch_col=batch_col)

            print()
            pp_underlined_hdg("Perturb-Cell sets")

            print("Perturb-Cell sets defined by the compound key:", f"{cell_type_key=}, {pert_col=}")
            print()

            # Get groups with actual observations (no empty groups)
            grouped = adata.obs.groupby([cell_type_key, pert_col], sort=False, observed=True)
            pert_cell_counts = grouped.size()

            # Number of perts per cell-type
            cell_type_n_perts = (adata.obs.groupby(cell_type_key, observed=True)[pert_col].nunique()
                                 .to_frame(name="Total_n_perts"))

            print(f"Nbr. pert-cell groups = {len(pert_cell_counts):,d}")
            print("Smallest pert-cell group:", pert_cell_counts.idxmin(), "=", format(pert_cell_counts.min(), ",d"))
            print("Largest pert-cell group:", pert_cell_counts.idxmax(), "=", format(pert_cell_counts.max(), ",d"))
            print()
            print(f"Mean pert-cell count   = {pert_cell_counts.mean():,.1f}")
            print(f"Median pert-cell count = {pert_cell_counts.median():,.1f}")
            print()

            print("Control cells:")
            ctl_mask = pert_cell_counts.index.get_level_values(pert_col) == control_pert

            print(f"    Nbr control cell groups = {ctl_mask.sum():,d}")
            print("    Control cell groups:\n\t", end="")
            counts_df = pert_cell_counts[ctl_mask].sort_index().to_frame("count")
            print(*counts_df.to_markdown(floatfmt=",.0f").splitlines(), sep="\n\t")
            print()

            for nmax in [10, 20, 50]:
                print(f"Nbr pert-cell groups with cell count < {nmax} = {(pert_cell_counts < nmax).sum():5,d}")

            print()

            if fewshot_perts or zeroshot_perts:
                if zeroshot_perts:
                    print()
                    pp_underlined_hdg("Held-out Perturbations: Zeroshot")

                    zeroshot_splits = defaultdict(list)
                    for ct, split in zeroshot_perts.items():
                        zeroshot_splits[split].append(ct)

                    for split, cts in zeroshot_splits.items():
                        mask = pert_cell_counts.index.get_level_values(cell_type_key).isin(cts)
                        counts_df = pert_cell_counts[mask].sort_index().to_frame(name='count')

                        print(f"Zeroshot {split} pert-cell sets = {len(counts_df):,d}")
                        print("", *counts_df.to_markdown(floatfmt=",.0f").splitlines(),
                              sep="\n\t")
                        print()

                if fewshot_perts:
                    print()
                    pp_underlined_hdg("Held-out Perturbations: Fewshot")

                    split_cell_types = defaultdict(list)
                    for ct, split_dict in fewshot_perts.items():
                        for split in split_dict.keys():
                            split_cell_types[split].append(ct)

                    split_cell_type_nperts = defaultdict(dict)

                    fewshot_splits = defaultdict(list)

                    for ct, split_dict in fewshot_perts.items():
                        if ct in pert_cell_counts.index.get_level_values(cell_type_key):
                            for split, perts in split_dict.items():
                                ct_split_perts = [(ct, p) for p in perts if (ct, p) in pert_cell_counts]
                                split_cell_type_nperts[split][ct] = len(ct_split_perts)
                                fewshot_splits[split].extend(ct_split_perts)

                    for split, ct_perts in fewshot_splits.items():
                        counts_df = pert_cell_counts.loc[ct_perts].sort_index().to_frame(name='count')

                        print(f"Fewshot {split} pert-cell sets = {len(counts_df):,d}")
                        print("", *counts_df.to_markdown(floatfmt=",.0f").splitlines(),
                              sep="\n\t")
                        print()

                    print("Nbr Few-shot pert-sets per cell-type and split:")
                    print("    Notes: Control perturbation is counted in 'train'.")
                    print("           Zeroshot counts are not included in this table.")
                    df = cell_type_n_perts.copy().reset_index()
                    df["train"] = df["Total_n_perts"]
                    for split in split_cell_type_nperts:
                        df[split] = df[cell_type_key].apply(lambda x: split_cell_type_nperts[split].get(x, 0))
                        df["train"] = df["train"] - df[split]
                    print("", *df.to_markdown(index=False, intfmt=",d").splitlines(),
                          sep="\n\t")
                    print()

            print("\n--- *** ---\n")

    print("========= ***** =========")
    return


# ======================================================================================================
#   Main
# ======================================================================================================

# To run
# ------
#
# [Python]$ python -m display_adata ppadata --ct "cell_type" --pert "cytokine" --control "PBS" --batch "donor"   \
#       /Users/smohan/Home/Projects/pertai/Data/scg-llm/parse1m_adata_hvg.h5ad
#
# [Python]$ python -m display_adata pptoml --ct "cell_type" --pert "cytokine" --control "PBS" --batch "donor"   \
#       /Users/smohan/Home/Projects/pertai/Data/scg-llm/parse1m_state.toml
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

    # ... ppadata [...] ADATA_PATH
    _sub_cmd_parser = _subparsers.add_parser('ppadata',
                                             help="Dispay summary of AnnData.")
    _sub_cmd_parser.add_argument('--ct', type=str, default=None,
                                 help="Column containing the Cell type.")
    _sub_cmd_parser.add_argument('--pert', type=str, default=None,
                                 help="Column containing the Perturbation label.")
    _sub_cmd_parser.add_argument('--control', type=str, default=None,
                                 help="Perturbation label value for Control cells.")
    _sub_cmd_parser.add_argument('--batch', type=str, default=None,
                                 help="Column containing the Batch label.")
    _sub_cmd_parser.add_argument('adata_path', type=str,
                                 help="Path to AnnData h5ad file.")

    # ... pptoml [...] TOML_PATH
    _sub_cmd_parser = _subparsers.add_parser('pptoml',
                                             help="Dispay summary of AnnData in TOML file.")
    _sub_cmd_parser.add_argument('--ct', type=str, default=None,
                                 help="Column containing the Cell type.")
    _sub_cmd_parser.add_argument('--pert', type=str, default=None,
                                 help="Column containing the Perturbation label.")
    _sub_cmd_parser.add_argument('--control', type=str, default=None,
                                 help="Perturbation label value for Control cells.")
    _sub_cmd_parser.add_argument('--batch', type=str, default=None,
                                 help="Column containing the Batch label.")
    _sub_cmd_parser.add_argument('toml_path', type=str,
                                 help="Path to TOML file.")

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

    if _args.subcmd == 'ppadata':

        pp_adata(_args.adata_path,
                 cell_type_key=_args.ct,
                 pert_col=_args.pert,
                 control_pert=_args.control,
                 batch_col=_args.batch
                 )

    elif _args.subcmd == 'pptoml':

        describe_toml(_args.toml_path,
                      cell_type_key=_args.ct,
                      pert_col=_args.pert,
                      control_pert=_args.control,
                      batch_col=_args.batch
                      )

    else:

        raise NotImplementedError(f"Command not implemented: {_args.subcmd}")

    # /

    print('\nTotal Run time =', datetime.now() - start_time_)
    print()
