from __future__ import annotations

try:
    from pyxtal import pyxtal

    HAS_PYXTAL: bool = True
except ImportError:
    pyxtal = None
    HAS_PYXTAL = False

try:
    import moyopy
    from moyopy.interface import MoyoAdapter

    HAS_MOYOPY: bool = True
except ImportError:
    moyopy = None
    MoyoAdapter = None
    HAS_MOYOPY = False
