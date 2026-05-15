from __future__ import annotations

import re
from collections import defaultdict
from string import ascii_uppercase
from typing import TYPE_CHECKING

from monty.fractions import gcd

from pymatgen.analysis.prototypes.wyckoff_utils import RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, split_alpha_numeric

if TYPE_CHECKING:
    from pymatgen.core import Composition


def get_prototype_formula_from_composition(composition: Composition) -> str:
    """An anonymized formula. Unique species are arranged in alphabetical order
    and assigned ascending alphabets. This format is used in the aflow structure
    prototype labelling scheme.

    Args:
        composition (Composition): Pymatgen Composition to process

    Returns:
        str: anonymized formula where the species are in alphabetical order
    """
    reduced = composition.element_composition
    return get_prototype_formula_from_counts([reduced[key] for key in sorted(reduced, key=str)])


def get_prototype_formula_from_counts(counts) -> str:
    """Get an AFLOW-style anonymous formula from element counts."""
    amounts = list(counts)
    if all(amt == int(amt) for amt in amounts):
        divisor = gcd(*(int(amt) for amt in amounts))
        amounts = [amt / divisor for amt in amounts]

    anon = ""
    for elem, amt in zip(ascii_uppercase, amounts, strict=False):
        if amt == 1:
            amt_str = ""
        elif abs(amt % 1) < 1e-8:
            amt_str = str(int(amt))
        else:
            amt_str = str(amt)
        anon += f"{elem}{amt_str}"
    return anon


def get_anonymous_formula_from_prototype_formula(prototype_formula: str) -> str:
    """Get an anonymous formula from a prototype formula."""
    prototype_formula = re.sub(RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, prototype_formula)
    anom_list = split_alpha_numeric(prototype_formula)

    return "".join(
        f"{el}{num}" if num != 1 else el
        for el, num in zip(anom_list["alpha"], sorted(map(int, anom_list["numeric"])), strict=True)
    )


def get_formula_from_protostructure_label(protostructure_label: str) -> str:
    """Get a formula from a protostructure label."""
    aflow_label, chemsys = protostructure_label.split(":")
    prototype_formula = aflow_label.split("_")[0]
    prototype_formula = re.sub(RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, prototype_formula)
    anom_list = split_alpha_numeric(prototype_formula)

    return "".join(
        f"{el}{num}" if num != 1 else el
        for el, num in zip(chemsys.split("-"), map(int, anom_list["numeric"]), strict=True)
    )


def _get_anonymous_formula_dict(anonymous_formula: str) -> dict:
    """Get a dictionary of element to count from an anonymous formula."""
    result: defaultdict = defaultdict(int)
    element = ""
    count = ""

    for char in anonymous_formula:
        if char.isalpha():
            if element:
                result[element] += int(count) if count else 1
                count = ""
            element = char
        else:
            count += char

    if element:
        result[element] += int(count) if count else 1

    return dict(result)
