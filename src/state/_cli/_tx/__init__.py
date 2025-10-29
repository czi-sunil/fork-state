import argparse as ap

from ._infer import add_arguments_infer, run_tx_infer
from ._predict import add_arguments_predict, run_tx_predict
from ._preprocess_infer import add_arguments_preprocess_infer, run_tx_preprocess_infer
from ._preprocess_train import add_arguments_preprocess_train, run_tx_preprocess_train
from ._train import add_arguments_train, run_tx_train

__all__ = [
    "run_tx_train",
    "run_tx_predict",
    "run_tx_infer",
    "run_tx_preprocess_train",
    "run_tx_preprocess_infer",
    "add_arguments_tx",
    "add_arguments_tx2",
]


def add_arguments_tx(parser: ap.ArgumentParser):
    """"""
    subparsers = parser.add_subparsers(required=True, dest="subcommand")
    add_arguments_train(subparsers.add_parser("train", add_help=False))
    add_arguments_predict(subparsers.add_parser("predict"))
    add_arguments_infer(subparsers.add_parser("infer"))
    add_arguments_preprocess_train(subparsers.add_parser("preprocess_train"))
    add_arguments_preprocess_infer(subparsers.add_parser("preprocess_infer"))


# [Sunil] Added func:
def add_arguments_tx2(parser: ap.ArgumentParser):
    """
    Copy of the "tx command",
       # that has a required argument: CONFIG_FILE.
    Remaining args are the same.
    """
    subparsers = parser.add_subparsers(required=True, dest="subcommand")

    # Subcommand: tx2 train
    train_parser = subparsers.add_parser("train", add_help=False)

    # required arg
    # train_parser.add_argument('config_file', type=str,
    #                           help="Path to Hydra Config file for this run.")

    add_arguments_train(train_parser)

    return
