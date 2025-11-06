"""
[Sunil] New file, for pre-processing data that has already been split
into Trainining and Test, using simple stratified splits on ALL perturbation labels.
"""

import argparse as ap
import logging
from pathlib import Path

import anndata as ad
import scanpy as sc
import hdf5plugin


# ------------------------------------------------------------------------------
# Globals
# ------------------------------------------------------------------------------


logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.INFO)


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------


def add_arguments_preprocess_splits(parser: ap.ArgumentParser):
    """Add arguments for the preprocess_train subcommand."""
    parser.add_argument(
        "--train_split",
        type=str,
        required=True,
        help="Path to the Training split AnnData file (.h5ad).",
    )
    parser.add_argument(
        "--test_split",
        type=str,
        required=True,
        help="Path to the Test split AnnData file (.h5ad).",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="Path to output dir where preprocessed AnnData files (.h5ad) will be written. Default is input dir.",
    )
    parser.add_argument(
        "--num_hvgs",
        type=int,
        default=2000,
        help="Number of highly variable genes to select (default = 2000).",
    )
    parser.add_argument(
        "--target_sum",
        type=int,
        default=1e4,
        help="Target library size to scale each cell to (default = 1e4).",
    )


def run_tx_preprocess_splits(train_split: str, test_split: str, output_dir: str, num_hvgs: int, target_sum: int):
    """
    Preprocess data that has already been split into Trainining and Test using simple stratified splits
    on ALL perturbation labels.
        - Scale every cell to `target_sum` (different from taking median count over the whole dataset)
        - Log-1p normalize
        - Get HVG's from Training-split (the larger split)
        - Add HVG info to .obsm['X_hvg']
        - Project adata.X to just the HVG's (to make data smaller)

    Processed data will look like state._cli._tx._preprocess_train.run_tx_preprocess_train.
    """
    # sanity check
    assert num_hvgs > 10

    # -- Get the Training split, use it to determine the HVGs, and process it

    adata = get_scaled_normalized_data(train_split, target_sum)

    logger.info(f"Finding top {num_hvgs} highly variable genes in Training split")
    sc.pp.highly_variable_genes(adata, n_top_genes=num_hvgs)

    # This is a boolean mask, as a pd.Series
    hvgs_mask = adata.var.highly_variable
    all_gene_names = adata.var["gene_symbol"]
    num_cells_trng = adata.n_obs

    logger.info("Storing highly variable genes in .obsm['X_hvg']")
    adata.obsm["X_hvg"] = adata[:, hvgs_mask].X.toarray()

    logger.info("Projecting adata.X to only the HVGs")
    adata = adata[:, adata.var.highly_variable]

    write_preprocessed_adata(adata, train_split, output_dir)

    # -- Process the Test split

    adata = get_scaled_normalized_data(test_split, target_sum)

    if adata.n_obs > num_cells_trng:
        print(f"! Warning: nbr cells in Test-split ({adata.n_obs}) > Training-split ({num_cells_trng}).")
        print(" ... Note that the Training split is used to determine the HVGs.")
        print(flush=True)

    assert adata.var["gene_symbol"].equals(all_gene_names), \
        f"Test split gene-names in `adata.var['gene_symbol']` are different from Training split!"

    adata = post_process_adata(adata, hvgs_mask)

    write_preprocessed_adata(adata, test_split, output_dir)

    return


def get_scaled_normalized_data(adata_path: str, target_sum: int) -> ad.AnnData:

    # sanity check
    assert target_sum > 10

    logger.info(f"Loading AnnData from {adata_path}")
    adata = ad.read_h5ad(adata_path)

    logger.info("Normalizing total counts per cell")
    sc.pp.normalize_total(adata, target_sum=target_sum)

    logger.info("Applying log1p transformation")
    sc.pp.log1p(adata)

    return adata


def post_process_adata(adata, hvgs_mask):
    logger.info("Post-procssing data to HVGs")
    adata.obsm["X_hvg"] = adata[:, hvgs_mask].X.toarray()

    adata = adata[:, hvgs_mask]

    return adata


def write_preprocessed_adata(adata, adata_path: str, outdir: str | None):
    adata_path = Path(adata_path)
    output_filename = adata_path.name.replace(adata_path.suffix, "_preprocessed" + adata_path.suffix)

    if outdir is not None:
        output_path = Path(outdir).joinpath(output_filename)
    else:
        output_path = adata_path.with_name(output_filename)

    logger.info(f"Saving preprocessed data to {output_path}")
    adata.write_h5ad(output_path, compression=hdf5plugin.FILTERS["zstd"])

    return
