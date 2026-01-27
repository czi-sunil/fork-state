"""
Test cell-load for use with pseudo-bulked Training and Validation data,
but Test data is not pseudo-bliked.
"""

from datetime import datetime
import logging
import os

from test_loader import create_anndata, summarize_loader, get_data_module


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


def create_dummy_zeroshot(outdir: str):
    """
    Creates 2 anndata files:
        Trng file: Training and Validation data
        Test file: Only Test data

    TOML File: Zero-shot entry for Test data.
    """

    zershot_toml_path = f"{outdir}/zeroshot.toml"

    assert not os.path.exists(zershot_toml_path), f"Zero-shot TOML already exists: {zershot_toml_path}"

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print("Output dir created:", outdir)

    # ---
    def write_split(x_adata, x_split, x_cell_type):
        out_file_basename = f"{x_cell_type}_{x_split}.h5ad"
        out_path = f"{outdir}/{out_file_basename}"

        print("Writing:", out_path)
        print(f"    Dataset shape: {x_adata.shape}")
        print(f"    nbr Perturb labels (incl. 'non-targeting') = {x_adata.obs['target_gene'].nunique()}")

        x_adata.write_h5ad(out_path)

        # return os.path.abspath(out_path)
        return out_file_basename
    # ---

    # --- Create the AnnData files and write

    ct_trng = "CT_1"
    ct_test = "CT_1"    # Test same cell type as Trng

    adata_trng = create_anndata(ct_trng)
    adata_test = create_anndata(ct_test, n_perts=4, n_cells_per_pert=5, base_count=100)

    # trng_perts = set(adata_trng) - {"non-targeting"}
    # test_perts = set(adata_test) - {"non-targeting"}

    trng_file = os.path.basename(write_split(adata_trng, "trng", ct_trng))
    test_file = os.path.basename(write_split(adata_test, "test", ct_test))

    # --- Write a basic README file
    with open(f"{outdir}/README_zeroshot.txt", "w") as f:
        print("Test dir and files created by `Arc/state/src/test_pseudobulk.create_dummy_zeroshot()`", file=f)
        print("for testing zero-shot files with Arc/cell-load", file=f)
        print("Files:", file=f)
        print(f"    {trng_file:<30}   Training and Validation data", file=f)
        print(f"    {test_file:<30}   Test data", file=f)
        print(f"    {os.path.basename(zershot_toml_path):<30}   TOML file", file=f)

    # --- Write the TOML

    with open(zershot_toml_path, "w") as f:
        print("[datasets]", file=f)
        trng_name = f'{ct_trng}_trng'
        test_name = f'{ct_test}_test'
        print(f'{trng_name} = "${{toml_dir}}/{trng_file}"', file=f)
        print(f'{test_name} = "${{toml_dir}}/{test_file}"', file=f)

        print(file=f)

        print("# All cell types in these datasets automatically go into training (excl. zero-/few-shot overrides)",
              file=f)
        print("[training]", file=f)
        print(f'{trng_name} = "train"', file=f)

        print(file=f)

        print("# Zeroshot specifications - entire cell types go to val or test", file=f)
        print("[zeroshot]", file=f)
        print(f'"{test_name}.{ct_test}" = "test"', file=f)

        print(file=f)

        print("# Fewshot specifications - explicit perturbation lists", file=f)
        print("[fewshot]", file=f)
        print(file=f)
        print(f'[fewshot."{trng_name}.{ct_trng}"]', file=f)
        print("val = ['pert_1', 'pert_2']", file=f)

    print("TOML written:", zershot_toml_path)

    return


def test_all_loaders(toml_path: str = "tmp/zeroshot.toml", log_info=False):

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
# [src]$ python -m test_pseudobulk create ./tmp
# [src]$ python -m test_pseudobulk loaders tmp/zeroshot.toml
#
#

if __name__ == "__main__":

    import argparse

    # --- Defaults

    DEFAULT_N_BATCHES = 5
    DEFAULT_NUM_WORKERS = 1
    # ---

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

    # ... timeloader [-n N_BATCHES]
    _sub_cmd_parser = _subparsers.add_parser('timeloader',
                                             help="Time Test of data loader.")
    _sub_cmd_parser.add_argument("-n", '--n_batches', type=int, default=DEFAULT_N_BATCHES,
                                 help=f"Nbr of batches (default {DEFAULT_N_BATCHES}).")
    _sub_cmd_parser.add_argument("-w", '--num_workers', type=int, default=DEFAULT_NUM_WORKERS,
                                 help=f"Nbr of workers (default {DEFAULT_NUM_WORKERS}).")
    _sub_cmd_parser.add_argument('dataset', type=str,
                                 help="Name of the dataset ('xaira' or 'vcc') or path to toml file.")

    # ...

    _args = _argparser.parse_args()
    # .................................................................................................

    start_time_ = datetime.now()

    print("---------------------------------------------------------------------")

    if _args.subcmd == 'create':

        create_dummy_zeroshot(_args.outdir)

    elif _args.subcmd == 'loaders':

        test_all_loaders(_args.toml, log_info=_args.log_info)

    else:

        raise NotImplementedError(f"Command not implemented: {_args.subcmd}")

    # /

    print('\nTotal Run time =', datetime.now() - start_time_)
    print()
