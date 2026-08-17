import json
import os
from importlib.resources import files

_DEFAULTS = {
    "sequence": {
        "minimum_average_phred_quality": 8,
        "min_sequence_length": 900,
        "max_sequence_length": 1500,
        "sequence_headcrop": 0,
        "sequence_tailcrop": 0,
        "barcode_window_size": 150,
    },
    "amplicones": {
        "barcode_length": 24,
        "spacer_length": 15,
        "linker_length": 15,
    },
    "max_barcode_distance": 5,
}


def _deep_merge(base, overrides):
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[key] = value


def load_config(path=None):
    """Return a config dict, optionally overridden by a user JSON file at `path`."""
    import copy
    config = copy.deepcopy(_DEFAULTS)

    if path and os.path.exists(path):
        with open(path) as f:
            _deep_merge(config, json.load(f))

    # Resolve bundled barcodes FASTA if the user hasn't specified a custom one
    if "barcodes_path" not in config["amplicones"]:
        bundled = files("torchlex.data").joinpath("barcodes_custom2.fasta")
        config["amplicones"]["barcodes_path"] = str(bundled)

    return config
