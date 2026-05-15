# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "pyxtal==1.1.3",
# ]
# ///
"""Generate Wyckoff position data used by pymatgen prototype labels.

Run with:

    uv run dev_scripts/generate_wyckoff_position_data.py

This regenerates the three compressed JSON files in
`src/pymatgen/analysis/prototypes`: Wyckoff multiplicities, positional degrees
of freedom, and equivalent Wyckoff-letter relabelings. The PEP 723 dependency
block above pins the crystallographic data source for reproducible runs.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable

DEFAULT_OUTPUT_DIR = Path("src/pymatgen/analysis/prototypes/assets")
GENERIC_FREE_COORDS = (0.123, 0.234, 0.345)
RE_WYCKOFF_LABEL = re.compile(r"(?P<multiplicity>\d+)(?P<letter>[A-Za-z])")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Wyckoff position data for pymatgen prototypes.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where wyckoff-position-*.json.gz files are written.",
    )
    parser.add_argument(
        "--space-groups",
        nargs="*",
        type=int,
        default=range(1, 231),
        help="Space group numbers to generate. Defaults to all 230 space groups.",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    space_groups = tuple(args.space_groups)

    multiplicities, params = generate_wyckoff_multiplicities_and_params(space_groups)
    relabelings = generate_wyckoff_relabelings(space_groups)

    dump_json_gz(multiplicities, output_dir / "wyckoff-position-multiplicities.json.gz")
    dump_json_gz(params, output_dir / "wyckoff-position-params.json.gz")
    dump_json_gz(relabelings, output_dir / "wyckoff-position-relabelings.json.gz")


def generate_wyckoff_multiplicities_and_params(
    space_groups: Iterable[int] = range(1, 231),
) -> tuple[dict[str, dict[str, int]], dict[str, dict[str, int]]]:
    """Generate Wyckoff multiplicities and positional degrees of freedom."""
    from pyxtal.symmetry import Group

    space_groups = tuple(space_groups)
    multiplicities = {}
    params = {}
    for spg_num in space_groups:
        group = Group(spg_num)
        multiplicities[str(spg_num)] = {}
        params[str(spg_num)] = {}
        for wp in group:
            multiplicity, letter = get_wyckoff_multiplicity_and_letter(wp)
            multiplicities[str(spg_num)][letter] = multiplicity
            params[str(spg_num)][letter] = int(wp.get_dof())
    if 0 in space_groups or space_groups == tuple(range(1, 231)):
        multiplicities["0"] = {}
    return multiplicities, params


def generate_wyckoff_relabelings(space_groups: Iterable[int] = range(1, 231)) -> dict[str, list[dict[str, str]]]:
    """Generate equivalent Wyckoff-letter relabelings for each space group."""
    from pyxtal.symmetry import Group

    relabelings = {}
    for spg_num in space_groups:
        group = Group(spg_num)
        letters = [get_wyckoff_multiplicity_and_letter(wp)[1] for wp in group]
        relabelings[str(spg_num)] = get_equivalent_letter_maps(group, letters)
    return relabelings


def get_equivalent_letter_maps(group: Any, letters: list[str]) -> list[dict[str, str]]:
    """Map each Wyckoff letter to its equivalent labels under origin shifts."""
    mappings = []
    candidate_shifts = {tuple(0.0 for _ in range(group.dim))}
    for wp in group:
        candidate_shifts.add(_get_wp_fixed_translation(wp))

    for shift in candidate_shifts:
        mapping = get_equivalent_letter_map(group, shift)
        if mapping is not None:
            mappings.append(mapping)

    mappings = sorted(
        {tuple(sorted(mapping.items())): mapping for mapping in mappings}.values(),
        key=lambda mapping: tuple(mapping[str(ord(letter))] for letter in sorted(letters)),
    )
    identity = {str(ord(letter)): letter for letter in letters}
    if not mappings:
        mappings.append(identity)
    elif identity not in mappings:
        mappings.insert(0, identity)
    return mappings


def get_equivalent_letter_map(group: Any, shift: tuple[float, ...]) -> dict[str, str] | None:
    mapping = {}
    for wp in group:
        multiplicity, letter = get_wyckoff_multiplicity_and_letter(wp)
        shifted_position = tuple(
            (coord + offset) % 1 for coord, offset in zip(_get_wp_position(wp), shift, strict=True)
        )
        shifted_wp = group.get_wyckoff_position_from_xyz(shifted_position, decimals=4)
        if shifted_wp is None:
            try:
                shifted_multiplicity, shifted_letter = get_equivalent_letter_from_unique_position(group, wp)
            except ValueError:
                return None
        else:
            shifted_multiplicity, shifted_letter = get_wyckoff_multiplicity_and_letter(shifted_wp)
        if shifted_multiplicity != multiplicity:
            return None
        mapping[str(ord(letter))] = shifted_letter
    return mapping


def get_equivalent_letter_from_unique_position(group: Any, wp: Any) -> tuple[int, str]:
    multiplicity, _ = get_wyckoff_multiplicity_and_letter(wp)
    dof = int(wp.get_dof())
    candidates = [
        get_wyckoff_multiplicity_and_letter(candidate)
        for candidate in group
        if get_wyckoff_multiplicity_and_letter(candidate)[0] == multiplicity and int(candidate.get_dof()) == dof
    ]
    if len(candidates) != 1:
        raise ValueError(f"Could not identify unique equivalent position for {get_wyckoff_label(wp)!r}")
    return candidates[0]


def _get_wp_position(wp: Any) -> tuple[float, ...]:
    return tuple(round(float(coord) % 1, 8) for coord in wp.get_position_from_free_xyzs(GENERIC_FREE_COORDS))


def _get_wp_fixed_translation(wp: Any) -> tuple[float, ...]:
    return tuple(round(float(coord) % 1, 8) for coord in wp[0].translation_vector)


def get_wyckoff_multiplicity_and_letter(wp: Any) -> tuple[int, str]:
    label = get_wyckoff_label(wp)
    match = RE_WYCKOFF_LABEL.search(label)
    if match is None:
        raise ValueError(f"Could not parse Wyckoff label from {label!r}")
    return int(match.group("multiplicity")), match.group("letter")


def get_wyckoff_label(wp: Any) -> str:
    if hasattr(wp, "get_label"):
        return str(wp.get_label())
    if hasattr(wp, "multiplicity") and hasattr(wp, "letter"):
        return f"{wp.multiplicity}{wp.letter}"
    return str(wp)


def dump_json_gz(data: dict, path: Path) -> None:
    with (
        open(path, "wb") as raw_file,
        gzip.GzipFile(fileobj=raw_file, mode="wb", mtime=0) as gzip_file,
        io.TextIOWrapper(gzip_file, encoding="utf-8") as file,
    ):
        json.dump(data, file, sort_keys=True, separators=(",", ":"))
        file.write("\n")


if __name__ == "__main__":
    main()
