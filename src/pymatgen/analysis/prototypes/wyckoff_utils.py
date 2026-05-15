from __future__ import annotations

import re
from itertools import groupby
from string import digits

from pymatgen.analysis.prototypes._data import WYCKOFF_POSITION_RELAB_DICT

remove_digits = str.maketrans("", "", digits)

RE_WYCKOFF_NO_PREFIX = re.compile(r"((?<![0-9])[A-z])")
RE_ELEMENT_NO_SUFFIX = re.compile(r"([A-z](?![0-9]))")
RE_WYCKOFF = re.compile(r"(?<!\d)([a-zA-Z])")
RE_ANONYMOUS = re.compile(r"([A-Z])(?![0-9])")
RE_SUBST_ONE_PREFIX = r"1\g<1>"
RE_SUBST_ONE_SUFFIX = r"\g<1>1"


def split_alpha_numeric(s: str) -> dict[str, list[str]]:
    """Split a string into separate lists of alpha and numeric groups.

    Args:
        s (str): The input string to split.

    Returns:
        dict[str, list[str]]: A dictionary with keys 'alpha' and 'numeric',
                              each containing a list of the respective groups.
    """
    groups = ["".join(g) for _, g in groupby(s, str.isalpha)]
    return {
        "alpha": [g for g in groups if g.isalpha()],
        "numeric": [g for g in groups if g.isnumeric()],
    }


def count_values_for_wyckoff(
    element_wyckoffs: list[str],
    counts: list[str],
    spg_num: str,
    lookup_dict: dict[str, dict[str, int]],
):
    """Count values from a lookup table and scale by wyckoff multiplicities."""
    return sum(
        int(count) * lookup_dict[spg_num][wyckoff_letter]
        for count, wyckoff_letter in zip(counts, element_wyckoffs, strict=True)
    )


def canonicalize_element_wyckoffs(element_wyckoffs: str, spg_num: int | str) -> str:
    """Given an element ordering, canonicalize the associated Wyckoff positions
    based on the alphabetical weight of equivalent choices of origin.

    Args:
        element_wyckoffs (str): wyckoff substring section from aflow_label with the
            wyckoff letters for different elements separated by underscores.
        spg_num (int | str): International space group number.

    Returns:
        str: element_wyckoff string with canonical ordering of the wyckoff letters.
    """
    isopointal_element_wyckoffs = list(
        {element_wyckoffs.translate(str.maketrans(trans)) for trans in WYCKOFF_POSITION_RELAB_DICT[str(spg_num)]}
    )

    scored_element_wyckoffs = [
        sort_and_score_element_wyckoffs(element_wyckoffs) for element_wyckoffs in isopointal_element_wyckoffs
    ]

    return min(scored_element_wyckoffs, key=lambda x: (x[1], x[0]))[0]


def sort_and_score_element_wyckoffs(element_wyckoffs: str) -> tuple[str, int]:
    """Determines the order or Wyckoff positions when canonicalizing AFLOW labels.

    Args:
        element_wyckoffs (str): wyckoff substring section from aflow_label with the
            wyckoff letters for different elements separated by underscores.

    Returns:
        tuple: containing
        - str: sorted Wyckoff position substring for AFLOW-style prototype label
        - int: integer score to rank order when canonicalizing
    """
    score = 0
    sorted_element_wyckoffs = []
    for el_wyks in element_wyckoffs.split("_"):
        wp_counts = split_alpha_numeric(el_wyks)
        sorted_element_wyckoffs.append(
            "".join(
                f"{count}{wyckoff_letter}" if count != "1" else wyckoff_letter
                for count, wyckoff_letter in sorted(
                    zip(wp_counts["numeric"], wp_counts["alpha"], strict=False),
                    key=lambda x: x[1],
                )
            )
        )
        score += sum(0 if wyckoff_letter == "A" else ord(wyckoff_letter) - 96 for wyckoff_letter in wp_counts["alpha"])

    return "_".join(sorted_element_wyckoffs), score
