"""
Misc utilities for pytorch
"""

import numpy as np

try:
    import torch
    TORCH_FOUND = True
except ModuleNotFoundError:
    TORCH_FOUND = False

if TORCH_FOUND:
    import torch.distributed

    # Pytorch support for Apple MPS?
    try:
        import torch.backends.mps
        PYTORCH_MPS_SUPPORTED = True
    except ModuleNotFoundError:
        PYTORCH_MPS_SUPPORTED = False


# -----------------------------------------------------------------------------
#   Functions
# -----------------------------------------------------------------------------


def pp_arr(msg, arr, print_value=False, indent: str = "", end=""):
    dtype = None
    shape = None
    device = None
    if arr is not None:
        if isinstance(arr, torch.Tensor):
            dtype = arr.type() if isinstance(arr, torch.Tensor) else arr.dtype
            shape = arr.shape
            device = arr.device
        elif isinstance(arr, np.ndarray):
            dtype = arr.dtype
            shape = arr.shape
        else:
            dtype = type(arr[0]) if isinstance(arr, (list, tuple)) else type(arr)
            shape = len(arr) if isinstance(arr, (list, tuple)) else None

    print(f"{indent}{msg}: type = {type(arr)}, {shape = }, {dtype = }", end="")
    if device is not None:
        print(f", {device = }")
    else:
        print()

    if print_value:
        for line in str(arr).splitlines():
            print(indent, line, sep="")

    print(end, end="")
    return


def is_non_scalar(x):
    return isinstance(x, torch.Tensor) or isinstance(x, np.ndarray) or isinstance(x, (list, tuple, dict))


def pp_arr_struct(tnsr_name, arr_struct, print_value=False, indent: str = "", extra_indent: str = "  "):
    pp_arr(tnsr_name, arr_struct, print_value=print_value, indent=indent)
    if isinstance(arr_struct, (list, tuple)):
        print(indent, "- size: ", len(arr_struct), sep="")
        if not is_non_scalar(arr_struct[0]):
            all_types = set(type(x) for x in arr_struct)
            if len(all_types) == 1:
                i = 0
                pp_arr(f"{tnsr_name}[{i}]", arr_struct[i], print_value=print_value, indent=indent + extra_indent)
                return

        for i in range(len(arr_struct)):
            pp_arr(f"{tnsr_name}[{i}]", arr_struct[i], print_value=print_value, indent=indent + extra_indent)

    elif isinstance(arr_struct, dict):
        print(indent, "- Keys: ", ", ".join(arr_struct.keys()), sep="")
        for k, v in arr_struct.items():
            pp_arr_struct(f"{tnsr_name}[{k}]", v, print_value=print_value, indent=indent + extra_indent)
    return
