"""
This module is intended to match crystal structures against known crystallographic "prototype"
structures.

The module also contains functions for getting the protostructure labels from a
variety of symmetry detection libraries (spglib, moyopy, aflow-sym). The protostructure
label is defined as the canonicalized aflow label with the alphabetically sorted
chemical system appended - `aflow_sym_label:chemsys`.

The utilities for determining the protostructure label are upstreamed from the
`aviary` package (https://github.com/CompRhys/aviary). If using these functions,
please cite the following publications:

Goodall, R. E., Parackal, A. S., Faber, F. A., Armiento, R., & Lee, A. A. (2022).
Rapid discovery of stable materials by coordinate-free coarse graining.
Science advances, 8(30), eabn4117. https://doi.org/10.1126/sciadv.abn4117

Parackal, A. S., Goodall, R. E., Faber, F. A., & Armiento, R. (2024).
Identifying crystal structures beyond known prototypes from x-ray powder diffraction spectra.
Physical Review Materials, 8(10), 103801. https://doi.org/10.1103/PhysRevMaterials.8.103801
"""

from __future__ import annotations

from pymatgen.analysis.prototypes._data import (
    AFLOW_PROTOTYPE_LIBRARY,
    CRYSTAL_FAMILY_SYMBOLS,
    CRYSTAL_LATTICE_PARAMETERS_COUNTS,
    MODULE_DIR,
    WYCKOFF_MULTIPLICITY_DICT,
    WYCKOFF_POSITION_PARAM_DICT,
    WYCKOFF_POSITION_RELAB_DICT,
    WYCKOFF_POSITION_SPLIT_DICT,
)
from pymatgen.analysis.prototypes._optional import HAS_MOYOPY, HAS_PYXTAL
from pymatgen.analysis.prototypes.formula import (
    _get_anonymous_formula_dict,
    get_anonymous_formula_from_prototype_formula,
    get_formula_from_protostructure_label,
    get_prototype_formula_from_composition,
    get_prototype_formula_from_counts,
)
from pymatgen.analysis.prototypes.matcher import AflowPrototypeMatcher, PrototypeDatabaseMatcher
from pymatgen.analysis.prototypes.protostructure import (
    AflowPrototypeLabel,
    ProtostructureLabel,
    _count_from_dict,
    _find_translations,
    _get_all_wyckoffs_substring_and_element_dict,
    count_crystal_dof,
    count_crystal_sites,
    count_distinct_wyckoff_letters,
    count_wyckoff_positions,
    get_protostructure_label,
    get_protostructure_label_from_aflow,
    get_protostructure_label_from_moyopy,
    get_protostructure_label_from_spg_analyzer,
    get_protostructure_label_from_spglib,
    get_protostructures_from_aflow_label_and_composition,
    get_prototype_from_protostructure,
    get_random_structure_for_protostructure,
    parse_aflow_prototype_label,
    parse_protostructure_label,
)
from pymatgen.analysis.prototypes.saps import (
    WyckoffSplitBranch,
    get_saps_prototypes_from_aflow_label,
    validate_wyckoff_split_branch,
)
from pymatgen.analysis.prototypes.wyckoff_utils import (
    RE_ANONYMOUS,
    RE_ELEMENT_NO_SUFFIX,
    RE_SUBST_ONE_PREFIX,
    RE_SUBST_ONE_SUFFIX,
    RE_WYCKOFF,
    RE_WYCKOFF_NO_PREFIX,
    canonicalize_element_wyckoffs,
    count_values_for_wyckoff,
    remove_digits,
    sort_and_score_element_wyckoffs,
    split_alpha_numeric,
)

__all__ = [
    "AFLOW_PROTOTYPE_LIBRARY",
    "CRYSTAL_FAMILY_SYMBOLS",
    "CRYSTAL_LATTICE_PARAMETERS_COUNTS",
    "HAS_MOYOPY",
    "HAS_PYXTAL",
    "MODULE_DIR",
    "RE_ANONYMOUS",
    "RE_ELEMENT_NO_SUFFIX",
    "RE_SUBST_ONE_PREFIX",
    "RE_SUBST_ONE_SUFFIX",
    "RE_WYCKOFF",
    "RE_WYCKOFF_NO_PREFIX",
    "WYCKOFF_MULTIPLICITY_DICT",
    "WYCKOFF_POSITION_PARAM_DICT",
    "WYCKOFF_POSITION_RELAB_DICT",
    "WYCKOFF_POSITION_SPLIT_DICT",
    "AflowPrototypeLabel",
    "AflowPrototypeMatcher",
    "ProtostructureLabel",
    "PrototypeDatabaseMatcher",
    "WyckoffSplitBranch",
    "_count_from_dict",
    "_find_translations",
    "_get_all_wyckoffs_substring_and_element_dict",
    "_get_anonymous_formula_dict",
    "canonicalize_element_wyckoffs",
    "count_crystal_dof",
    "count_crystal_sites",
    "count_distinct_wyckoff_letters",
    "count_values_for_wyckoff",
    "count_wyckoff_positions",
    "get_anonymous_formula_from_prototype_formula",
    "get_formula_from_protostructure_label",
    "get_protostructure_label",
    "get_protostructure_label_from_aflow",
    "get_protostructure_label_from_moyopy",
    "get_protostructure_label_from_spg_analyzer",
    "get_protostructure_label_from_spglib",
    "get_protostructures_from_aflow_label_and_composition",
    "get_prototype_formula_from_composition",
    "get_prototype_formula_from_counts",
    "get_prototype_from_protostructure",
    "get_random_structure_for_protostructure",
    "get_saps_prototypes_from_aflow_label",
    "parse_aflow_prototype_label",
    "parse_protostructure_label",
    "remove_digits",
    "sort_and_score_element_wyckoffs",
    "split_alpha_numeric",
    "validate_wyckoff_split_branch",
]
