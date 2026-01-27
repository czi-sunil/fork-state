"""
# [Sunil] Added this file
"""

import json
import yaml


def deep_diff_dicts(d1, d2):
    """
    Recursively calculates the deep difference between two dictionaries.

    Returns a dictionary of differences, including:
    - 'added': keys present in d2 but not d1
    - 'removed': keys present in d1 but not d2
    - 'changed': keys present in both but with different values (nested diff for dicts)
    """
    diff = {'added': {}, 'removed': {}, 'changed': {}}

    # Find added and removed keys
    keys1 = set(d1.keys())
    keys2 = set(d2.keys())

    added_keys = keys2 - keys1
    removed_keys = keys1 - keys2
    common_keys = keys1 & keys2

    for key in added_keys:
        diff['added'][key] = d2[key]
    for key in removed_keys:
        diff['removed'][key] = d1[key]

    # Find changed values
    for key in common_keys:
        v1 = d1[key]
        v2 = d2[key]

        if v1 != v2:
            if isinstance(v1, dict) and isinstance(v2, dict):
                # Recurse for nested dictionaries
                nested_diff = deep_diff_dicts(v1, v2)
                if nested_diff:
                    diff['changed'][key] = nested_diff
            else:
                # Value changed (not a nested dict)
                diff['changed'][key] = {'old': v1, 'new': v2}

    # Filter out empty difference categories
    return {k: v for k, v in diff.items() if v}


def cmp_configs(ref_config, new_config):
    print()
    print("Delta between configs:")
    print("   old:", ref_config)
    print("   new:", new_config)
    print()

    with open(ref_config) as cfg:
        ref_cfg_dict = yaml.safe_load(cfg)

    with open(new_config) as cfg:
        new_cfg_dict = yaml.safe_load(cfg)

    print(json.dumps(deep_diff_dicts(ref_cfg_dict, new_cfg_dict), indent=4))
    print()

    return


# ======================================================================================================
#   Main
# ======================================================================================================

# To run
# ------
#
# [Python]$ python -m state.utils cmpconfig REF_CONFIG.yaml NEW_CONFIG.yaml
#
#

if __name__ == "__main__":

    import argparse
    from datetime import datetime

    _argparser = argparse.ArgumentParser(
        description='Comparing YAML files.',
    )

    _subparsers = _argparser.add_subparsers(dest='subcmd',
                                            title='Available commands',
                                            )
    # Make the sub-commands required
    _subparsers.required = True

    # ... cmpconfig REF_CONFIG.yaml NEW_CONFIG.yaml
    _sub_cmd_parser = _subparsers.add_parser('cmpconfig',
                                             help="Compare two config.yaml files.")
    _sub_cmd_parser.add_argument('config_ref', type=str,
                                 help="Path to reference config.yaml file.")
    _sub_cmd_parser.add_argument('config_new', type=str,
                                 help="Path to new config.yaml file.")

    # ...

    _args = _argparser.parse_args()
    # .................................................................................................

    start_time_ = datetime.now()

    print("---------------------------------------------------------------------")

    if _args.subcmd == 'cmpconfig':

        cmp_configs(_args.config_ref, _args.config_new)

    else:

        raise NotImplementedError(f"Command not implemented: {_args.subcmd}")

    # /

    print('\nTotal Run time =', datetime.now() - start_time_)
    print()
