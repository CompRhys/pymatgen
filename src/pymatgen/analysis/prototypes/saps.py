from __future__ import annotations

import re
from itertools import combinations
from typing import Literal, NotRequired, TypedDict

from pymatgen.analysis.prototypes._data import WYCKOFF_MULTIPLICITY_DICT, WYCKOFF_POSITION_SPLIT_DICT
from pymatgen.analysis.prototypes.formula import get_prototype_formula_from_counts
from pymatgen.analysis.prototypes.protostructure import (
    AflowPrototypeLabel,
    _count_from_dict,
    parse_aflow_prototype_label,
)
from pymatgen.analysis.prototypes.wyckoff_utils import (
    RE_SUBST_ONE_PREFIX,
    RE_WYCKOFF_NO_PREFIX,
    canonicalize_element_wyckoffs,
    split_alpha_numeric,
)
from pymatgen.util.due import Doi, due


class WyckoffSplitBranch(TypedDict):
    """Serializable subgroup relation for Wyckoff position splitting."""

    parent_space_group: int
    child_space_group: int
    relation_type: Literal["t", "k"]
    index: int
    transformation: list[list[int | float]]
    origin_shift: list[int | float]
    splits: dict[str, list[str]]
    parent_setting: NotRequired[str]
    child_setting: NotRequired[str]
    branch_id: NotRequired[str]
    child_pearson_symbol: NotRequired[str]


def validate_wyckoff_split_branch(branch: WyckoffSplitBranch) -> None:
    """Validate multiplicity conservation for a split branch."""
    parent_sg = str(branch["parent_space_group"])
    child_sg = str(branch["child_space_group"])
    index = int(branch["index"])
    parent_mults = WYCKOFF_MULTIPLICITY_DICT[parent_sg]
    child_mults = WYCKOFF_MULTIPLICITY_DICT[child_sg]

    for parent_letter, child_letters in branch["splits"].items():
        parent_count = parent_mults[parent_letter] * index
        child_count = sum(child_mults[letter] for letter in child_letters)
        if parent_count != child_count:
            raise ValueError(
                f"Invalid split for SG {parent_sg} {parent_letter} -> {child_sg} "
                f"{child_letters}: expected multiplicity {parent_count}, got {child_count}"
            )


@due.dcite(
    Doi("10.1038/s41586-023-06735-9"),
    description="Scaling deep learning for materials discovery.",
)
def get_saps_prototypes_from_aflow_label(
    aflow_label: str,
    max_subgroup_index: int | None = None,
    relation_types: tuple[Literal["t", "k"], ...] = ("t", "k"),
    wyckoff_split_table: dict[str, list[WyckoffSplitBranch]] | None = None,
) -> list[str]:
    """Enumerate SAPS-inspired child prototype labels from an AFLOW prototype label.

    This prototype-level implementation is inspired by the GNoME SAPS workflow,
    but is unlikely to be a faithful reproduction of the Google DeepMind method.
    The reference table describes subgroup relations, and each generated label
    uses the child space group and child Wyckoff letters from a single branch.
    """
    parsed = parse_aflow_prototype_label(aflow_label)
    split_table = wyckoff_split_table if wyckoff_split_table is not None else WYCKOFF_POSITION_SPLIT_DICT
    branches = split_table.get(parsed.space_group, [])
    saps_labels = set()

    for branch in branches:
        if branch["relation_type"] not in relation_types:
            continue
        if max_subgroup_index is not None and int(branch["index"]) > max_subgroup_index:
            continue
        if "child_pearson_symbol" not in branch:
            continue

        validate_wyckoff_split_branch(branch)
        split_wyckoffs = [_split_element_wyckoffs(wyckoffs, branch["splits"]) for wyckoffs in parsed.element_wyckoffs]

        for element_idx, child_wyckoffs in enumerate(split_wyckoffs):
            if len(child_wyckoffs) < 2:
                continue
            for subset_size in range(1, len(child_wyckoffs)):
                for subset in combinations(range(len(child_wyckoffs)), subset_size):
                    child_parts = [list(wyckoffs) for wyckoffs in split_wyckoffs]
                    replacement = [child_parts[element_idx][idx] for idx in subset]
                    child_parts[element_idx] = [
                        letter for idx, letter in enumerate(child_parts[element_idx]) if idx not in subset
                    ]
                    child_parts.insert(element_idx + 1, replacement)
                    if all(child_parts):
                        saps_labels.add(_build_child_aflow_label(branch, child_parts))

    return sorted(saps_labels)


def _split_element_wyckoffs(element_wyckoffs: str, splits: dict[str, list[str]]) -> list[str]:
    """Rewrite one element's parent Wyckoff substring in the child subgroup."""
    sep_el_wyks = split_alpha_numeric(re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, element_wyckoffs))
    child_wyckoffs = []
    for count, parent_letter in zip(sep_el_wyks["numeric"], sep_el_wyks["alpha"], strict=True):
        split_letters = splits[parent_letter]
        for _ in range(int(count)):
            child_wyckoffs.extend(split_letters)
    return child_wyckoffs


def _build_child_aflow_label(branch: WyckoffSplitBranch, child_parts: list[list[str]]) -> str:
    child_sg = str(branch["child_space_group"])
    child_wyckoff_parts = [_collapse_wyckoff_letters(wyckoffs) for wyckoffs in child_parts]
    child_wyckoffs = re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, "_".join(child_wyckoff_parts))
    canonical_wyckoffs = canonicalize_element_wyckoffs(child_wyckoffs, child_sg)
    formula_counts = [
        _count_from_dict([wyckoffs], WYCKOFF_MULTIPLICITY_DICT, child_sg) for wyckoffs in child_wyckoff_parts
    ]
    formula = get_prototype_formula_from_counts(formula_counts)
    return str(
        AflowPrototypeLabel(
            prototype_formula=formula,
            pearson_symbol=branch["child_pearson_symbol"],
            space_group=child_sg,
            element_wyckoffs=tuple(canonical_wyckoffs.split("_")),
        )
    )


def _collapse_wyckoff_letters(wyckoff_letters: list[str]) -> str:
    return "".join(
        f"{count}{letter}" if count > 1 else letter
        for letter, count in sorted({letter: wyckoff_letters.count(letter) for letter in set(wyckoff_letters)}.items())
    )
