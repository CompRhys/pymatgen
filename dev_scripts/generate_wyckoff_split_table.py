# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "pymatgen==2025.10.7",
#   "pyxtal==1.1.3",
# ]
# ///
"""Generate subgroup Wyckoff-splitting reference data for SAPS.

Run with:

    uv run dev_scripts/generate_wyckoff_split_table.py

This regenerates `src/pymatgen/analysis/prototypes/wyckoff-position-splits.json.gz`
from PyXtal's subgroup and Wyckoff-splitting tables. Use `--space-groups` to
limit parent groups and `--max-index` to restrict subgroup branches. The PEP
723 dependency block above pins the data-generation dependencies.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import numbers
from collections.abc import Iterable
from pathlib import Path
from typing import Any

DEFAULT_PROTOTYPE_DIR = Path("src/pymatgen/analysis/prototypes/assets")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate prototype subgroup Wyckoff split reference data.")
    parser.add_argument(
        "output",
        nargs="?",
        default=DEFAULT_PROTOTYPE_DIR / "wyckoff-position-splits.json.gz",
        help="Output JSON or JSON.GZ file.",
    )
    parser.add_argument("--space-groups", nargs="*", type=int, help="Parent space groups to export.")
    parser.add_argument("--relation-types", nargs="*", default=("t",), choices=("t", "k"))
    parser.add_argument("--max-index", type=int, help="Maximum subgroup index to include.")
    parser.add_argument("--prototype-dir", type=Path, default=DEFAULT_PROTOTYPE_DIR)
    args = parser.parse_args()

    table = generate_wyckoff_split_table(
        parent_space_groups=args.space_groups or get_aflow_parent_space_groups(args.prototype_dir),
        relation_types=args.relation_types,
        max_index=args.max_index,
    )
    dump_json_gz(table, Path(args.output))


def get_aflow_parent_space_groups(prototype_dir: Path = DEFAULT_PROTOTYPE_DIR) -> list[int]:
    """Get parent space groups represented in the shipped AFLOW prototype library."""
    aflow_prototypes = load_json_gz(prototype_dir / "aflow_prototypes.json.gz")

    spg_nums = {
        int(entry["tags"]["aflow"].split("_")[2]) for entry in aflow_prototypes if entry.get("tags", {}).get("aflow")
    }
    return sorted(spg_nums)


def generate_wyckoff_split_table(
    parent_space_groups: Iterable[int],
    relation_types: Iterable[str] = ("t",),
    max_index: int | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """Generate subgroup Wyckoff split records using PyXtal's crystallographic tables."""
    try:
        from pyxtal.symmetry import Group
    except ImportError as exc:
        raise ImportError("pyxtal is required to generate Wyckoff split reference data") from exc

    table: dict[str, list[dict[str, Any]]] = {}
    for parent_sg in parent_space_groups:
        group = Group(parent_sg)
        branches = []
        for relation_type in relation_types:
            subgroup_data = _get_subgroup_data(group, relation_type)
            subgroup_numbers = subgroup_data.get("subgroup", [])
            for idx, child_sg in enumerate(subgroup_numbers):
                if max_index is not None and _get_subgroup_index(parent_sg, child_sg, subgroup_data, idx) > max_index:
                    continue
                branch = _build_branch(
                    group,
                    parent_sg,
                    int(child_sg),
                    relation_type,
                    subgroup_data,
                    idx,
                )
                if branch["splits"]:
                    branches.append(branch)
        if branches:
            table[str(parent_sg)] = branches
    return table


def _get_subgroup_data(group, relation_type: str) -> dict[str, Any]:
    method_name = "get_max_t_subgroup" if relation_type == "t" else "get_max_k_subgroup"
    method = getattr(group, method_name)
    return method()


def _build_branch(group, parent_sg: int, child_sg: int, relation_type: str, subgroup_data, idx: int) -> dict:
    split_map = {}
    relations = _extract_list(subgroup_data, "relations", idx, default=[])
    for parent_wp, child_wps in zip(group, relations, strict=False):
        parent_letter = _get_wyckoff_label(parent_wp)[-1]
        child_letters = _letters_from_obj(child_wps)
        if child_letters:
            split_map[parent_letter] = child_letters

    return {
        "parent_space_group": parent_sg,
        "child_space_group": child_sg,
        "relation_type": relation_type,
        "index": _get_subgroup_index(parent_sg, child_sg, subgroup_data, idx),
        "branch_id": str(idx),
        "transformation": _extract_list(subgroup_data, "transformation", idx, default=[]),
        "origin_shift": _extract_list(subgroup_data, "origin_shift", idx, default=[]),
        "child_pearson_symbol": _extract_list(subgroup_data, "pearson", idx, default=""),
        "splits": split_map,
    }


def _get_wyckoff_label(wp) -> str:
    if hasattr(wp, "get_label"):
        return str(wp.get_label())
    return str(wp)


def _letters_from_obj(obj) -> list[str]:
    if isinstance(obj, str):
        return [obj[-1]]
    if isinstance(obj, Iterable):
        letters = []
        for item in obj:
            if isinstance(item, str):
                letters.append(item[-1])
            elif hasattr(item, "letter"):
                letters.append(item.letter)
            elif hasattr(item, "get_label"):
                letters.append(item.get_label()[-1])
            else:
                letters.extend(_letters_from_obj(item))
        return letters
    if hasattr(obj, "letter"):
        return [obj.letter]
    if hasattr(obj, "get_label"):
        return [obj.get_label()[-1]]
    return []


def _get_subgroup_index(parent_sg: int, child_sg: int, subgroup_data: dict, idx: int) -> int:
    index = _extract_list(subgroup_data, "index", idx)
    if index:
        return int(index)
    parent_order = _extract_list(subgroup_data, "order", idx, default=0)
    child_order = _extract_list(subgroup_data, "subgroup_order", idx, default=0)
    if parent_order and child_order:
        return int(parent_order / child_order)
    return 1 if parent_sg == child_sg else 2


def _extract_list(data: dict, key: str, idx: int, default=None):
    value = data.get(key, default)
    if isinstance(value, list | tuple) and idx < len(value):
        return value[idx]
    return value


def load_json_gz(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as file:
        return json.load(file)


def dump_json_gz(data: dict, path: Path) -> None:
    with (
        open(path, "wb") as raw_file,
        gzip.GzipFile(fileobj=raw_file, mode="wb", mtime=0) as gzip_file,
        io.TextIOWrapper(gzip_file, encoding="utf-8") as file,
    ):
        json.dump(_jsonable(data), file, sort_keys=True, separators=(",", ":"))
        file.write("\n")


def _jsonable(obj):
    if isinstance(obj, dict):
        return {key: _jsonable(val) for key, val in obj.items()}
    if isinstance(obj, list | tuple):
        return [_jsonable(val) for val in obj]
    if isinstance(obj, numbers.Integral):
        return int(obj)
    if isinstance(obj, numbers.Real):
        return float(obj)
    if hasattr(obj, "tolist"):
        return _jsonable(obj.tolist())
    return obj


if __name__ == "__main__":
    main()
