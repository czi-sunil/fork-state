"""
For testing the cell-load PerturbationDataModule
"""

from collections import defaultdict
import logging
import os

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

from sklearn.model_selection import train_test_split

import torch

from cell_load.data_modules.perturbation_dataloader import PerturbationDataModule


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


def get_data_module(toml_path: str = "../../../Data/Arc/vcc_curated/test.toml",
                    cell_set_sz=4,
                    batch_sz=2,
                    cache_perturbation_control_pairs=False):

    print("Loading toml file:", toml_path)
    print()

    pdm = PerturbationDataModule(toml_path,
                                 batch_size=batch_sz,
                                 pert_col="target_gene",
                                 batch_col="batch",
                                 cell_type_key="tissue_ontology_term_id",
                                 control_pert="non-targeting",
                                 embed_key=None,
                                 output_space="gene",
                                 basal_mapping_strategy="random",
                                 cell_sentence_len=cell_set_sz,
                                 num_workers=1,
                                 cache_perturbation_control_pairs=cache_perturbation_control_pairs,
                                 barcode=False,
                                 )

    return pdm


def create_anndata(cell_type: str = "CT_1", n_perts=8, n_cells_per_pert=10, n_genes=10):

    step = 1 / n_genes
    n_cells = n_perts * n_cells_per_pert
    counts = np.arange(1.0, n_cells + 1, step, dtype=np.float32).reshape(n_cells, n_genes)
    adata = ad.AnnData(csr_matrix(counts))

    adata.obs_names = [f"Cell_{i+1:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i+1:d}" for i in range(adata.n_vars)]

    # cell_type
    adata.obs["tissue_ontology_term_id"] = pd.Categorical(np.repeat(cell_type, adata.n_obs))
    # perturbation
    pert_labels = ["non-targeting"] + [f"pert_{i+1:d}" for i in range(n_perts - 1)]
    adata.obs["target_gene"] = pd.Categorical(np.repeat(pert_labels, n_cells_per_pert))
    # batch ... all the same
    adata.obs["batch"] = pd.Categorical(np.repeat("batch", adata.n_obs))

    # HVG's ... just take the top 5
    num_hvgs = min(5, n_genes // 2)

    # ... copy and pre-process for `sc.pp.highly_variable_genes`. This will create all the metadata.
    adata2 = adata.copy()
    sc.pp.normalize_total(adata2)
    sc.pp.log1p(adata2)
    sc.pp.highly_variable_genes(adata2, n_top_genes=num_hvgs)

    adata.var = adata2.var
    # Mark only the first 5 as HVG
    adata.var['highly_variable'][:num_hvgs] = True
    adata.var['highly_variable'][num_hvgs:] = False

    adata.obsm["X_hvg"] = adata[:, adata2.var.highly_variable].X.toarray()

    return adata


def split_adata_stratified(adata, test_pct=0.3):
    """
    Every perturbation, including 'non-targeting', is split based on `test_pct`.
    """

    pert_col = "target_gene"

    train_df, test_df = train_test_split(adata.obs,
                                         test_size=test_pct,
                                         stratify=adata.obs[pert_col],
                                         random_state=0,
                                         )
    test_cells_mask = adata.obs.index.isin(test_df.index)

    adata_test = adata[test_cells_mask, :].copy()
    adata_trng = adata[np.logical_not(test_cells_mask), :].copy()

    return adata_trng, adata_test


# noinspection PyTypeChecker
def create_dummy_splits(outdir: str):
    """
    Creates 2 cell-types. Data written to `outdir`, ONLY IF it does not already exist.

    The first cell-type gets a stratified-split into trng/test.
    For validation: the first 2 non-ctrl perturbations are selected FROM THE TEST DATA.

    The 2nd ct is allocated completely to trng.

    Note: No validation split!
    """

    assert not os.path.exists(outdir), f"Dir already exists: {outdir}"

    # ---
    def write_split(x_adata, x_split, x_cell_type):
        out_path = f"{outdir}/{x_cell_type}_{x_split}.h5ad"

        print("Writing:", out_path)
        print(f"    Dataset shape: {x_adata.shape}")
        print(f"    nbr Perturb labels (incl. 'non-targeting') = {x_adata.obs['target_gene'].nunique()}")

        x_adata.write_h5ad(out_path)

        return os.path.abspath(out_path)
    # ---

    pert_col = "target_gene"

    all_perts = dict()
    test_perts = dict()
    split_file_paths = defaultdict(dict)
    split_cts = ["CT_1"]

    outdir_created = False

    # --- Create AnnData, Splits, and write

    for ct in ["CT_1", "CT_2"]:
        adata = create_anndata(ct)

        all_perts[ct] = set(adata.obs[pert_col]) - {"non-targeting"}

        if not outdir_created:
            os.mkdir(outdir)
            outdir_created = True

        split_file_paths[ct]["all"] = write_split(adata, "all", ct)

        # Split only CT_1
        if ct not in split_cts:
            continue

        adata_trng, adata_test = split_adata_stratified(adata)

        test_perts[ct] = set(adata_test.obs[pert_col]) - {"non-targeting"}

        split_file_paths[ct]["trng"] = write_split(adata_trng, "trng", ct)
        split_file_paths[ct]["test"] = write_split(adata_test, "test", ct)

    # --- Write the TOML

    toml_path = f"{outdir}/split1.toml"
    with open(toml_path, "w") as f:
        print("[datasets]", file=f)
        for ct, ct_files in split_file_paths.items():
            if ct in split_cts:
                write_splits = ["trng", "test"]
            else:
                write_splits = ["all"]

            for split in write_splits:
                print(f'{ct}_{split} = "{ct_files[split]}"', file=f)

        print(file=f)

        print("[training]", file=f)
        for ct, ct_files in split_file_paths.items():
            if ct in split_cts:
                split = "trng"
            else:
                split = "all"
            print(f'{ct}_{split} = "train"', file=f)

        print(file=f)

        print("[zeroshot]", file=f)
        print(file=f)
        print("[fewshot]", file=f)
        print(file=f)
        for ct, test_perts in test_perts.items():
            print(f'[fewshot."{ct}_test.{ct}"]', file=f)

            val_perts = sorted(test_perts - {"non-targeting"})[:2]
            print("val =", val_perts, file=f)

            # print("test =", sorted(test_perts), file=f)
            print("test = ", '"_rest_"', file=f)

            print(file=f)

        print("TOML written:", toml_path)

    return


def summarize_loader(name, dl):
    print("------------------------------------")
    print("Processing data loader:", name)
    print()

    data = defaultdict(lambda: defaultdict(list))

    # -- Run thru the DataLoader and collect all the data

    n_batches = 0
    tot_n_vals = 0
    for batch in dl:
        n_batches += 1
        tot_n_vals += batch["pert_cell_emb"].shape[0]

        cell_types = np.asarray(batch["cell_type"])
        uniq_cell_types = set(batch["cell_type"])
        for k, v in batch.items():
            if k == "cell_type":
                continue
            if isinstance(v, list):
                v = np.asarray(v)
            for ct in uniq_cell_types:
                data[ct][k].append(v[cell_types == ct])

    # -- Determine unique values

    uniq_vals = defaultdict(lambda: defaultdict(set))

    pert_cell_embs = dict()
    pert_names = dict()

    for ct, ct_vdict in data.items():
        for k, v in ct_vdict.items():
            for vv in v:
                if isinstance(vv, torch.Tensor):
                    vv = vv.numpy()

                if k.endswith("cell_emb"):
                    uniq_vals[ct][k].update(vv[:, 0].tolist())
                elif isinstance(vv, np.ndarray) and len(vv.shape) == 2:
                    vals = [", ".join(str(x) for x in row) for row in vv]
                    try:
                        uniq_vals[ct][k].update(vals)
                    except TypeError as e:
                        print(f"Error in {k=}, {type(vv)=}, {vv.shape=}, {vv.dtype=}")
                        raise e
                else:
                    uniq_vals[ct][k].update(vv)

        pert_cell_embs[ct] = np.concatenate([x[:, 0] for x in ct_vdict["pert_cell_emb"]])
        pert_names[ct] = np.concatenate(ct_vdict["pert_name"])

    # -- pp summary

    print()
    print("nbr Batches =", n_batches)
    print("nbr samples =", tot_n_vals)
    print("Cell types  =", ", ".join(sorted(uniq_vals.keys())))
    print_flds = ("pert_cell_emb, ctrl_cell_emb, pert_name, "
                  "batch, batch_name, pert_cell_barcode, ctrl_cell_barcode").split(", ")
    for k in print_flds:
        for ct in sorted(uniq_vals.keys()):
            vals = sorted(uniq_vals[ct][k])
            print()
            print(f"{ct}/{k}: nbr uniq values = {len(vals):3d}")
            if not vals:
                continue
            try:
                print(f"vals: {vals[0]}", *vals[1:], sep="\n      ")
            except IndexError as e:
                print(f"Error in {ct}/{k=}")
                raise e

            if k == "pert_cell_emb":
                print("Nbr uniq vals per `pert_name`:")
                for pert_name in sorted(uniq_vals[ct]["pert_name"]):
                    pert_name_embs = pert_cell_embs[ct][pert_names[ct] == pert_name]
                    print(f"    {pert_name} = {len(np.unique(pert_name_embs))}")

    print()
    return


def test_all_loaders(toml_path: str = "tmp/split2.toml", log_info=False):

    if log_info:
        logging.basicConfig(level=logging.INFO)

    print()
    pdm = get_data_module(toml_path)
    pdm.setup("fit")

    print("\n")

    summarize_loader("Validation", pdm.val_dataloader())
    summarize_loader("Test", pdm.test_dataloader())
    summarize_loader("Training", pdm.train_dataloader())

    return


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
        description='Test cell-load PerturbationDataModule Data-Loaders.',
    )

    _subparsers = _argparser.add_subparsers(dest='subcmd',
                                            title='Available commands',
                                            )
    # Make the sub-commands required
    _subparsers.required = True

    # ... create OUTDIR
    _sub_cmd_parser = _subparsers.add_parser('create',
                                             help="Create and write data to OUTDIR.")
    _sub_cmd_parser.add_argument('outdir', type=str,
                                 help="Path to dir, which is created and data written there.")

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

    start_time = datetime.now()

    print("---------------------------------------------------------------------")

    if _args.subcmd == 'create':

        create_dummy_splits(_args.outdir)

    elif _args.subcmd == 'loaders':

        test_all_loaders(_args.toml, log_info=_args.log_info)

    else:

        raise NotImplementedError(f"Command not implemented: {_args.subcmd}")

    # /

    print('\nTotal Run time =', datetime.now() - start_time)
    print()
