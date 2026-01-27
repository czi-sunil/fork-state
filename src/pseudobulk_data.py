"""
Test cell-load for use with pseudo-bulked Training and Validation data,
but Test data is not pseudo-bliked.
"""

from collections import defaultdict
from datetime import datetime
import logging
import os
from pathlib import Path
import shutil

import numpy as np
import anndata as ad
from anndata import AnnData

from cell_load.config import ExperimentConfig
from cell_load.data_modules import PerturbationDataModule
from cell_load.utils.data_utils import GlobalH5MetadataCache
from cell_load.utils.data_utils import H5MetadataCache

from test_loader import create_anndata, summarize_loader, get_data_module


# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------


def convert_to_pseudobulk(toml_path: str,
                          outdir: str = None,
                          n_cells_pb: int = 20,
                          skip_low_counts: bool = True,
                          cell_type_key: str = "cell_line",
                          pert_col: str = "gene",
                          control_pert: str = "non-targeting",
                          batch_col: str = "gem_group",
                          ):
    """
    Convert datasets described in TOML file to pseudo-bulked dataset.

    :param toml_path:
    :param outdir: Dir where to store the pseudo-bulked dataset.
        Default is subdir of dir containing toml file
    :param n_cells_pb: Nbr cells to average for each pseudo-bulked cell.
    :param skip_low_counts: pert-cell sets with count below `n_cells_pb` are not pseudo-bulked (left as is).
    :param cell_type_key:
    :param pert_col:
    :param control_pert:
    :param batch_col:
    :return:

    from pseudobulk_data import convert_to_pseudobulk
    toml_path = '~/Home/Projects/pertai/Data/Arc/Replogle-Nadig-Preprint/test.toml'
    test_data_info = convert_to_pseudobulk(toml_path)
    """

    toml_path = Path(toml_path).expanduser()

    if outdir is None:
        outdir = f"{toml_path.parent}/{toml_path.stem}_pbulk_n{n_cells_pb}"

    print()
    if os.path.isdir(outdir):
        print("Output dir already exists:", outdir)
        print("Clearing old contents and re-creating the directory...")
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    print("Output dir created:", outdir)
    print()

    # Arbitrary values, won't be used
    batch_size = 4
    cell_sentence_len = 16

    pdm = PerturbationDataModule(str(toml_path),
                                 batch_size=batch_size,
                                 pert_col=pert_col,
                                 batch_col=batch_col,
                                 cell_type_key=cell_type_key,
                                 control_pert=control_pert,
                                 embed_key=None,
                                 output_space="gene",
                                 basal_mapping_strategy="random",
                                 cell_sentence_len=cell_sentence_len,
                                 num_workers=1,
                                 cache_perturbation_control_pairs=False,
                                 barcode=False,
                                 )
    pdm.setup()

    # Close the adata files that pdm opened. We will reopen them later (o/w causes OS error).
    close_opened_adatas(pdm)

    print("\nPerturbationDataModule created.\n")

    config = pdm.config

    # For each file:
    #   IF there is a Test split THEN copy it (along with control) to a separate file
    #   Pseudo-bulk each (Cell-type, Perturbation) set (the adata.obsm['X_hvg'] data)

    # Files where test extracted data was written to
    test_data_info = []

    # Dataset-name => [pbulk_file, ...]
    pbulk_files = defaultdict(list)

    rng = np.random.default_rng(12345)

    for dataset_name in config.get_all_datasets():
        dataset_path = Path(config.datasets[dataset_name])
        # noinspection PyProtectedMember
        files = pdm._find_dataset_files(dataset_path)

        # stem, full-path
        for fname, fpath in files.items():

            print("Reading adata:", fpath)
            # Don't use `backed="r"` as later `write_h5ad` causes IndexError: Fancy indexing out of range for (0-5693)
            adata = ad.read_h5ad(fpath)
            print()

            if adata is None:
                print("Dataset not found:", fpath)
                continue

            # Extract and write the test data, if any
            if test_extract_info := extract_save_test_split(adata, fpath, fname,
                                                            config,
                                                            cell_type_key, pert_col, control_pert,
                                                            outdir):
                test_data_info.append(test_extract_info)

            print("Pseudo-bulking file:", fpath, "    ...")

            # Get the metadata cache created by PerturbationDataModule
            cache = GlobalH5MetadataCache().get_cache(str(fpath),
                                                      pdm.pert_col,
                                                      pdm.cell_type_key,
                                                      pdm.control_pert,
                                                      pdm.batch_col)

            # Process each cell type in this file
            for ct_idx, ct in enumerate(cache.cell_type_categories):
                ct_mask = cache.cell_type_codes == ct_idx
                n_cells = np.sum(ct_mask)

                if n_cells == 0:
                    continue

                # Pseudo-bulk each perturbation set in this cell-type
                for pert_code, pert_label in enumerate(cache.pert_categories):

                    pert_indices = np.where(np.logical_and(ct_mask, cache.pert_codes == pert_code))[0]

                    pseudobulk_pert_set(adata, pert_indices, n_cells_pb, skip_low_counts, rng)

            print()

            # Remove test data
            if test_extract_info:
                adata = remove_test_pert_cells(adata, test_extract_info, cache)

            # Write the pseudo-bulked data
            pb_fpath = f"{outdir}/{fpath.stem}_pbulk.h5ad"

            pbulk_files[dataset_name].append(pb_fpath)

            print()
            print("Writing pseudo-bulked data to:", pb_fpath)

            adata.write_h5ad(pb_fpath)

            print()
            print("Pseudo-bulked data written:", pb_fpath)
            print()

            if adata.isbacked:
                adata.file.close()

    # Now write the Pseudo-bulk TOML file

    write_pseudobulk_toml(outdir, toml_path, config, pbulk_files, test_data_info,
                          n_cells_pb, skip_low_counts, cell_type_key, pert_col, control_pert, batch_col)

    return test_data_info


def close_opened_adatas(pdm: PerturbationDataModule):
    # pdm has opened all the files, so get pointers to them
    ds_adata = dict()
    for ds in (pdm.test_datasets, pdm.val_datasets, pdm.train_datasets):
        for subset_ds in ds:
            # noinspection PyUnresolvedReferences
            ds_adata[subset_ds.dataset.h5_path] = subset_ds.dataset.h5_file

    # Now close these files
    for ds_path, ds_file in ds_adata.items():
        ds_file.close()

    return


def pseudobulk_pert_set(adata: AnnData,
                        pert_indices: np.ndarray,
                        n_cells_pb: int,
                        skip_low_counts: bool,
                        rng: np.random.Generator,
                        ):
    """
    `adata.obsm["X_hvg"]` is modified:
        Replaces `adata.obsm["X_hvg"][pert_indices]` with equal nbr of pseudo-bulked rows.

    :param adata:
    :param pert_indices:
    :param n_cells_pb:
    :param skip_low_counts:
    :param rng:
    :return:
    """

    n_cells = len(pert_indices)
    if n_cells == 0:
        return

    if skip_low_counts and n_cells < n_cells_pb:
        return

    pb_cells = []

    # Replace only if n_cells is small
    replace = n_cells < (1.5 * n_cells_pb)

    for _ in range(n_cells):
        pb_idxs = rng.choice(pert_indices, n_cells_pb, replace=replace)

        pb_row = adata.obsm["X_hvg"][pb_idxs].mean(axis=0)
        pb_cells.append(pb_row)

    pb_data = np.vstack(pb_cells)

    adata.obsm["X_hvg"][pert_indices] = pb_data

    return


def extract_save_test_split(adata: AnnData,
                            adata_path: Path,
                            dataset_name: str,
                            config: ExperimentConfig,
                            cell_type_key: str,
                            pert_col: str,
                            control_pert: str,
                            outdir: str,
                            ) -> dict | None:
    """
    Only extract the Test split, and write to another file.

    :return: For the Test data extracted:
        test_extract_info:
        {
            'zeroshot' => {dataset-name: [cell-type, ...]}
            'fewshot' =>  {dataset-name: { cell-type => [pert, ...] } }
                ... there can be more than one cell-type here.
            'dataset_path' => {dataset_name: path to new dataset}
        }
        or
        None, if no data extracted.
    """

    test_split_name = "test"

    # [1] Set up dicts indexed on split, for this dataset_name

    # split => [ct, ...]
    zeroshot_splits = defaultdict(list)
    # split => [(ct, pert), ...]
    fewshot_splits = defaultdict(list)

    if zeroshot_perts := config.get_zeroshot_celltypes(dataset_name):
        # zeroshot_perts: Cell-type => Split
        for ct, split in zeroshot_perts.items():
            zeroshot_splits[split].append(ct)

    if fewshot_perts := config.get_fewshot_celltypes(dataset_name):
        # fewshot_perts: Cell-type => Split => [ pert, ...]
        for ct, split_dict in fewshot_perts.items():
            for split, perts in split_dict.items():
                fewshot_splits[split].extend([(ct, p) for p in perts])
                # For fewsoht, also add control
                fewshot_splits[split].append((ct, control_pert))

    # Any `test` split in this dataset_name?
    if not zeroshot_splits[test_split_name] and not fewshot_splits[test_split_name]:
        return None

    print()
    print("Extracting test split from:", adata_path)

    test_adata_path = f"{outdir}/{adata_path.stem}_{test_split_name}{adata_path.suffix}"

    test_extract_info = dict()

    # adata may not have all (or any) of the required test data, since
    # the test data is defined by `dataset_name`, which could map to a group of AnnData files.

    test_adata = None

    # [2] test_adata := Zero-shot cell-types that are in this file

    if zs_cts := zeroshot_splits[test_split_name]:
        test_adata = adata[adata.obs[cell_type_key].isin(zs_cts)]

        if test_adata.shape[0] == 0:
            test_adata = None
            print(f"... No Zero-shot cell-types in this file.")
        else:
            test_zs_cts = test_adata.obs[cell_type_key].unique().tolist()
            test_extract_info['zeroshot'] = {dataset_name: test_zs_cts}
            print(f"... nbr Zero-shot cell-types extracted = {len(test_zs_cts)}")
            print("   ", test_zs_cts)

        print()

    # [3] Extract Few-shot perts that are in this file, append to test_adata

    if ct_perts := fewshot_splits[test_split_name]:

        grouped = adata.obs.groupby([cell_type_key, pert_col], sort=False, observed=True)

        test_indices = []
        # ct => [pert, ...]
        test_ct_perts = defaultdict(list)
        nex = 0

        for ct_pert_ in ct_perts:
            # Use `get` since (ct, pert) may not be in this file
            if (indices := grouped.indices.get(ct_pert_)) is not None:
                test_indices.append(indices)
                test_ct_perts[ct_pert_[0]].append(ct_pert_[1])
                nex += 1

        if test_indices:
            idxs = np.concatenate(test_indices)
            # sort the indices so order is same (shouldn't really matter)
            idxs.sort()
            fewshot_adata = adata[idxs]

            # Build the retval
            test_extract_info["fewshot"] = {dataset_name: test_ct_perts}

            print(f"... nbr Few-shot pert-cell sets extracted = {nex:,d}")
            print(f"    nbr Few-shot pert-cell sets not found = {len(ct_perts) - nex:,d}")
            print()

            if test_adata is not None:
                test_adata = ad.concat([test_adata, fewshot_adata], merge="same")
            else:
                test_adata = fewshot_adata

    if test_adata is None:
        print("... no test data found")
        return None

    test_adata.write_h5ad(test_adata_path)

    test_extract_info["dataset_path"] = {dataset_name: test_adata_path}

    print("Test split written to:", test_adata_path)
    print(f"    Dataset shape: {test_adata.shape}")
    print()

    return test_extract_info


def remove_test_pert_cells(adata: AnnData,
                           test_extract_info: dict,
                           cache: H5MetadataCache,
                           ) -> AnnData:

    # ---
    def get_ct_code(ct) -> int | None:
        idxs = np.nonzero(cache.cell_type_categories == ct)[0]
        if idxs.shape[0] > 0:
            return idxs[0]
        return None

    def get_pert_code(pert_label) -> int | None:
        idxs = np.nonzero(cache.pert_categories == pert_label)[0]
        if idxs.shape[0] > 0:
            return idxs[0]
        return None
    # ---

    if test_extract_info is None or not test_extract_info:
        return adata

    print("Removing test cells from:", Path(cache.h5_path).name)

    test_mask = np.zeros_like(cache.pert_codes, dtype=bool)

    # 'zeroshot' => {dataset-name: [cell-type, ...]}
    if zs_info := test_extract_info.get("zeroshot"):
        zs_cell_types = next(iter(zs_info.values()))
        for zs_ct in zs_cell_types:
            if (ct_idx := get_ct_code(zs_ct)) is not None:
                test_mask = np.logical_or(test_mask, cache.cell_type_codes == ct_idx)

    # 'fewshot' =>  {dataset-name: { cell-type => [pert, ...] } }
    if fs_info := test_extract_info.get("fewshot"):
        fs_ct_perts = next(iter(fs_info.values()))
        for fs_ct, ct_perts in fs_ct_perts.items():

            if (ct_idx := get_ct_code(fs_ct)) is not None:
                ct_mask = cache.cell_type_codes == ct_idx

                for pert in ct_perts:
                    if (pert_idx := get_pert_code(pert)) is not None:
                        # Don't remove Few-shot control cells
                        if pert_idx != cache.control_pert_code:
                            pert_mask = np.logical_and(ct_mask, cache.pert_codes == pert_idx)
                            test_mask = np.logical_or(test_mask, pert_mask)

    if test_mask.sum() > 0:
        print(f"Nbr cells removed = {test_mask.sum():,d}")
        return adata[np.logical_not(test_mask)]
    else:
        print("No test cells found")

    return adata


# noinspection PyUnusedLocal
def write_pseudobulk_toml(outdir: str,
                          toml_path: Path,
                          config: ExperimentConfig,
                          pbulk_files: dict[str, list[str]],
                          test_data_info: list[dict],
                          n_cells_pb: int,
                          skip_low_counts: bool,
                          cell_type_key: str,
                          pert_col: str,
                          control_pert: str,
                          batch_col: str,
                          ):

    pb_toml_path = f"{outdir}/{toml_path.stem}_pbulk.toml"

    print("Writing pseudobulk toml...")

    # ---
    def format_dataset_paths(fpaths: str | list[str]) -> str:
        # Assume all paths are files from outdir
        dir_pfx = "${toml_dir}/"

        if isinstance(fpaths, str) or len(fpaths) == 1:
            if isinstance(fpaths, str):
                fp = fpaths
            else:
                fp = fpaths[0]
            return f"{dir_pfx}{Path(fp).name}"

        fnames = [Path(fp).name for fp in fpaths]
        combo_path = f"{dir_pfx}{{{','.join(fnames)}}}"
        return combo_path
    # ---

    with open(pb_toml_path, "w") as f:
        print("# Pseudobulk version of:", toml_path, file=f)
        print("# Args:", file=f)
        print(f"#    {n_cells_pb =}", file=f)
        print(f"#    {skip_low_counts =}", file=f)
        print(f"#    cell_type_key = '{cell_type_key}'", file=f)
        print(f"#    pert_col = '{pert_col}'", file=f)
        print(f"#    control_pert = '{control_pert}'", file=f)
        print(f"#    batch_col = '{batch_col}'", file=f)

        print("[datasets]", file=f)
        for dsname, ds_files in pbulk_files.items():
            print(f'{dsname} = "{format_dataset_paths(ds_files)}"', file=f)

        # Specify the zs-test data
        for n, test_extract_info in enumerate(test_data_info, start=1):
            ds_dict = test_extract_info["dataset_path"]
            ds_name, ds_path = next(iter(ds_dict.items()))
            print(f'{ds_name}_test_{n} = "{format_dataset_paths(ds_path)}"', file=f)

        print(file=f)

        print("[training]", file=f)
        for dsname in pbulk_files.keys():
            print(f'{dsname} = "train"', file=f)

        print(file=f)

        print("# Zeroshot specifications - entire cell types go to val or test", file=f)
        print("[zeroshot]", file=f)
        # Specify the zs-test data
        for n, test_extract_info in enumerate(test_data_info, start=1):
            ds_dict = test_extract_info["dataset_path"]
            ds_name = next(iter(ds_dict.keys()))

            zs_cts = set()
            if zs_info := test_extract_info.get("zeroshot"):
                zs_cts = set(zs_info[ds_name])

            if fs_info := test_extract_info.get("fewshot"):
                zs_cts = zs_cts | set(fs_info[ds_name].keys())

            ds_name += f"_test_{n}"
            for ct in sorted(zs_cts):
                print(f'"{ds_name}.{ct}" = "test"', file=f)

    print("Pseudobulk TOML written to:", pb_toml_path)

    return pb_toml_path


# -----------------------------------------------------------------------------
# Functions: Testing Zero-shot
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


def test_zeroshot_loaders(toml_path: str = "tmp/zeroshot.toml", log_info=False):

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
# [src]$ python -m pseudobulk_data test_create ./tmp
# [src]$ python -m pseudobulk_data test_loaders tmp/zeroshot.toml
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

    # ... test_create OUTDIR
    _sub_cmd_parser = _subparsers.add_parser('test_create',
                                             help="Create and write data to OUTDIR.")
    _sub_cmd_parser.add_argument('outdir', type=str,
                                 help="Path to dir, which is created and data written there.")

    # ... test_loaders TOML_PATH
    _sub_cmd_parser = _subparsers.add_parser('test_loaders',
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

    if _args.subcmd == 'test_create':

        create_dummy_zeroshot(_args.outdir)

    elif _args.subcmd == 'test_loaders':

        test_zeroshot_loaders(_args.toml, log_info=_args.log_info)

    else:

        raise NotImplementedError(f"Command not implemented: {_args.subcmd}")

    # /

    print('\nTotal Run time =', datetime.now() - start_time_)
    print()
