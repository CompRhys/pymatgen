from __future__ import annotations

import os

from monty.serialization import loadfn

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


def _load_wyckoff_position_relab_dict():
    relab_dict = loadfn(f"{MODULE_DIR}/assets/wyckoff-position-relabelings.json.gz")
    return {
        spg_num: [{int(key): line for key, line in val.items()} for val in vals] for spg_num, vals in relab_dict.items()
    }


AFLOW_PROTOTYPE_LIBRARY = loadfn(f"{MODULE_DIR}/assets/aflow_prototypes.json.gz")
WYCKOFF_MULTIPLICITY_DICT = loadfn(f"{MODULE_DIR}/assets/wyckoff-position-multiplicities.json.gz")
WYCKOFF_POSITION_PARAM_DICT = loadfn(f"{MODULE_DIR}/assets/wyckoff-position-params.json.gz")
WYCKOFF_POSITION_RELAB_DICT = _load_wyckoff_position_relab_dict()
WYCKOFF_POSITION_SPLIT_DICT = loadfn(f"{MODULE_DIR}/assets/wyckoff-position-splits.json.gz")

CRYSTAL_FAMILY_SYMBOLS = {
    "triclinic": "a",
    "monoclinic": "m",
    "orthorhombic": "o",
    "tetragonal": "t",
    "trigonal": "h",
    "hexagonal": "h",
    "cubic": "c",
}

CRYSTAL_LATTICE_PARAMETERS_COUNTS = {
    "a": 6,
    "m": 4,
    "o": 3,
    "t": 2,
    "h": 2,
    "c": 1,
}
