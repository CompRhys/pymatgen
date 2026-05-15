from __future__ import annotations

import os
from typing import TYPE_CHECKING

from monty.serialization import loadfn

if TYPE_CHECKING:
    from collections.abc import Callable

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


class LazyLoad:
    """Lazy proxy for prototype data files."""

    def __init__(self, loader: Callable):
        self._loader = loader
        self._data = None

    @property
    def data(self):
        if self._data is None:
            self._data = self._loader()
        return self._data

    def __bool__(self):
        return bool(self.data)

    def __contains__(self, key):
        return key in self.data

    def __getitem__(self, key):
        return self.data[key]

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return repr(self.data)

    def get(self, *args, **kwargs):
        return self.data.get(*args, **kwargs)

    def items(self):
        return self.data.items()

    def keys(self):
        return self.data.keys()

    def values(self):
        return self.data.values()


def _load_wyckoff_position_relab_dict():
    relab_dict = loadfn(f"{MODULE_DIR}/wyckoff-position-relabelings.json.gz")
    return {
        spg_num: [{int(key): line for key, line in val.items()} for val in vals] for spg_num, vals in relab_dict.items()
    }


AFLOW_PROTOTYPE_LIBRARY = LazyLoad(lambda: loadfn(f"{MODULE_DIR}/assets/aflow_prototypes.json.gz"))
WYCKOFF_MULTIPLICITY_DICT = LazyLoad(lambda: loadfn(f"{MODULE_DIR}/assets/wyckoff-position-multiplicities.json.gz"))
WYCKOFF_POSITION_PARAM_DICT = LazyLoad(lambda: loadfn(f"{MODULE_DIR}/assets/wyckoff-position-params.json.gz"))
WYCKOFF_POSITION_RELAB_DICT = LazyLoad(_load_wyckoff_position_relab_dict)
WYCKOFF_POSITION_SPLIT_DICT = LazyLoad(lambda: loadfn(f"{MODULE_DIR}/assets/wyckoff-position-splits.json.gz"))

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
