from __future__ import annotations

import re
import subprocess
from collections import Counter
from dataclasses import dataclass
from itertools import chain, groupby, permutations, product
from operator import itemgetter
from shutil import which
from typing import TYPE_CHECKING

import orjson

from pymatgen.analysis.prototypes._data import (
    CRYSTAL_LATTICE_PARAMETERS_COUNTS,
    WYCKOFF_MULTIPLICITY_DICT,
    WYCKOFF_POSITION_PARAM_DICT,
    WYCKOFF_POSITION_RELAB_DICT,
)
from pymatgen.analysis.prototypes._optional import HAS_MOYOPY, HAS_PYXTAL, MoyoAdapter, moyopy, pyxtal
from pymatgen.analysis.prototypes.formula import (
    _get_anonymous_formula_dict,
    get_anonymous_formula_from_prototype_formula,
    get_prototype_formula_from_composition,
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
from pymatgen.core import Composition
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.due import Doi, due

if TYPE_CHECKING:
    from typing import Literal

    from pymatgen.core import Structure


@dataclass(frozen=True)
class AflowPrototypeLabel:
    """Parsed AFLOW prototype label."""

    prototype_formula: str
    pearson_symbol: str
    space_group: str
    element_wyckoffs: tuple[str, ...]

    @classmethod
    def from_str(cls, aflow_label: str) -> AflowPrototypeLabel:
        """Parse an AFLOW prototype label without a chemical-system suffix."""
        label = aflow_label.split(":", maxsplit=1)[0]
        prototype_formula, pearson_symbol, spg_num, *element_wyckoffs = label.split("_")
        return cls(
            prototype_formula=prototype_formula,
            pearson_symbol=pearson_symbol,
            space_group=spg_num,
            element_wyckoffs=tuple(element_wyckoffs),
        )

    def __str__(self) -> str:
        return "_".join((self.prototype_formula, self.pearson_symbol, self.space_group, *self.element_wyckoffs))


@dataclass(frozen=True)
class ProtostructureLabel:
    """Parsed protostructure label with chemical system."""

    aflow_label: AflowPrototypeLabel
    chemical_system: str

    @classmethod
    def from_str(cls, protostructure_label: str) -> ProtostructureLabel:
        aflow_label, chemical_system = protostructure_label.split(":", maxsplit=1)
        return cls(AflowPrototypeLabel.from_str(aflow_label), chemical_system)

    def __str__(self) -> str:
        return f"{self.aflow_label}:{self.chemical_system}"


def parse_aflow_prototype_label(aflow_label: str) -> AflowPrototypeLabel:
    """Parse an AFLOW prototype label without a chemical-system suffix."""
    return AflowPrototypeLabel.from_str(aflow_label)


def parse_protostructure_label(protostructure_label: str) -> ProtostructureLabel:
    """Parse a protostructure label with chemical-system suffix."""
    return ProtostructureLabel.from_str(protostructure_label)


@due.dcite(
    Doi("10.1126/sciadv.abn4117"),
    description="Rapid discovery of stable materials by coordinate-free coarse graining.",
)
def get_protostructure_label(
    struct: Structure,
    method: Literal["aflow", "spglib", "moyopy"],
    raise_errors: bool = False,
    **kwargs,
) -> str | None:
    """Get protostructure label for a pymatgen Structure.

    Args:
        struct (Structure): pymatgen Structure
        method (Literal["aflow", "spglib", "moyopy"]): Method to use for symmetry
            detection
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.
        **kwargs: Additional arguments for the specific method

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    if method == "aflow":
        return get_protostructure_label_from_aflow(struct, raise_errors, **kwargs)
    if method == "spglib":
        return get_protostructure_label_from_spglib(struct, raise_errors, **kwargs)
    if method == "moyopy":
        return get_protostructure_label_from_moyopy(struct, raise_errors, **kwargs)
    raise ValueError(f"Invalid method: {method}")


def get_protostructure_label_from_aflow(
    struct: Structure,
    raise_errors: bool = False,
    aflow_executable: str | None = None,
) -> str:
    """Get protostructure label for a pymatgen Structure. Make sure you're running a
    recent version of the aflow CLI as there's been several breaking changes. This code
    was tested under v3.2.12. The protostructure label is constructed as
    `aflow_label:chemsys`.

    Install guide: https://aflow.org/install-aflow/#install_aflow
        http://aflow.org/install-aflow/install-aflow.sh -o install-aflow.sh
        chmod 555 install-aflow.sh
        ./install-aflow.sh --slim

    Args:
        struct (Structure): pymatgen Structure
        aflow_executable (str): path to aflow executable. Defaults to which("aflow").
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    if aflow_executable is None:
        aflow_executable = which("aflow")

    if which(aflow_executable or "") is None:
        raise FileNotFoundError(
            "AFLOW could not be found, please specify path to its binary with aflow_executable='...'"
        )

    cmd = f"{aflow_executable} --prototype --print=json cat".split()

    output = subprocess.run(
        cmd,
        input=struct.to(fmt="poscar"),
        text=True,
        capture_output=True,
        check=True,
    )

    aflow_proto = orjson.loads(output.stdout)

    aflow_label = aflow_proto["aflow_prototype_label"]
    chemsys = struct.chemical_system
    prototype_form, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    element_dict = {}
    for elem, wyk_letters_per_elem in zip(chemsys.split("-"), element_wyckoffs, strict=True):
        wyk_letters_normalized = re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, wyk_letters_per_elem)
        sep_el_wyks = split_alpha_numeric(wyk_letters_normalized)
        element_dict[elem] = count_values_for_wyckoff(
            sep_el_wyks["alpha"],
            sep_el_wyks["numeric"],
            spg_num,
            WYCKOFF_MULTIPLICITY_DICT,
        )

    element_wyckoffs = "_".join(element_wyckoffs)
    element_wyckoffs = canonicalize_element_wyckoffs(element_wyckoffs, spg_num)

    protostructure_label = f"{prototype_form}_{pearson_symbol}_{spg_num}_{element_wyckoffs}:{chemsys}"

    observed_formula = Composition(element_dict).reduced_formula
    expected_formula = struct.composition.reduced_formula
    if observed_formula != expected_formula:
        err_msg = (
            f"Invalid WP multiplicities - {protostructure_label}, expected {observed_formula} to be {expected_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)

        return err_msg

    return protostructure_label


def _get_all_wyckoffs_substring_and_element_dict(
    equivalent_wyckoff_labels: list[tuple[int, str, str]],
    spg_num: int | str,
):
    """Get Wyckoff position substring and element dict from equivalent Wyckoff labels.

    Args:
        equivalent_wyckoff_labels (list[tuple[int, str, str]]): List of tuples containing
            (multiplicity, element symbol, Wyckoff letter).
        spg_num (int | str): Space group number.

    Returns:
        tuple[str, dict]: Tuple containing:
            - str: Wyckoff position substring
            - dict: Dictionary mapping element symbols to their multiplicities
    """
    equivalent_wyckoff_labels = sorted(equivalent_wyckoff_labels, key=lambda x: (x[1], x[2]))

    element_dict = {}
    element_wyckoffs = []
    for el, group in groupby(equivalent_wyckoff_labels, key=lambda x: x[1]):
        list_group = list(group)
        element_dict[el] = sum(WYCKOFF_MULTIPLICITY_DICT[str(spg_num)][e[2]] for e in list_group)
        element_wyckoffs.append(
            "".join(
                f"{len(list(occurrences))}{wyk_letter}"
                for wyk_letter, occurrences in groupby(list_group, key=lambda x: x[2])
            )
        )
    all_wyckoffs = "_".join(element_wyckoffs)
    all_wyckoffs = canonicalize_element_wyckoffs(all_wyckoffs, spg_num)

    return all_wyckoffs, element_dict


def get_protostructure_label_from_spg_analyzer(
    spg_analyzer: SpacegroupAnalyzer,
    raise_errors: bool = False,
) -> str:
    """Get protostructure label for pymatgen SpacegroupAnalyzer.

    Args:
        spg_analyzer (SpacegroupAnalyzer): pymatgen SpacegroupAnalyzer object.
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    sym_struct = spg_analyzer.get_symmetrized_structure()

    spg_num = spg_analyzer.get_space_group_number()
    pearson_symbol = spg_analyzer.get_pearson_symbol()
    prototype_form = get_prototype_formula_from_composition(sym_struct.composition)
    chemsys = sym_struct.chemical_system

    equivalent_wyckoff_labels = [
        (len(s), s[0].species_string, wyk_letter.translate(remove_digits))
        for s, wyk_letter in zip(sym_struct.equivalent_sites, sym_struct.wyckoff_symbols, strict=True)
    ]

    all_wyckoffs, element_dict = _get_all_wyckoffs_substring_and_element_dict(equivalent_wyckoff_labels, spg_num)

    protostructure_label = f"{prototype_form}_{pearson_symbol}_{spg_num}_{all_wyckoffs}:{chemsys}"

    observed_formula = Composition(element_dict).reduced_formula
    expected_formula = sym_struct.composition.reduced_formula
    if observed_formula != expected_formula:
        err_msg = (
            f"Invalid WP multiplicities - {protostructure_label}, expected {observed_formula} to be {expected_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)

        return err_msg

    return protostructure_label


def get_protostructure_label_from_spglib(
    struct: Structure,
    raise_errors: bool = False,
    init_symprec: float = 0.1,
    fallback_symprec: float | None = 1e-5,
) -> str:
    """Get AFLOW prototype label for pymatgen Structure.

    Args:
        struct (Structure): pymatgen Structure object.
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.
        init_symprec (float): Initial symmetry precision for spglib. Defaults to 0.1.
        fallback_symprec (float): Fallback symmetry precision for spglib if first
            symmetry detection failed. Defaults to 1e-5.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    attempt_to_recover = False
    try:
        spg_analyzer = SpacegroupAnalyzer(struct, symprec=init_symprec, angle_tolerance=5)
        try:
            aflow_label_with_chemsys = get_protostructure_label_from_spg_analyzer(spg_analyzer, raise_errors)

            if ("Invalid" in aflow_label_with_chemsys) and fallback_symprec is not None:
                attempt_to_recover = True
        except ValueError:
            if fallback_symprec is None:
                raise
            attempt_to_recover = True

        if attempt_to_recover:
            spg_analyzer = SpacegroupAnalyzer(
                spg_analyzer.get_refined_structure(),
                symprec=fallback_symprec,
                angle_tolerance=-1,
            )
            aflow_label_with_chemsys = get_protostructure_label_from_spg_analyzer(spg_analyzer, raise_errors)
        return aflow_label_with_chemsys  # type: ignore[possibly-undefined]

    except ValueError as exc:
        if not raise_errors:
            return str(exc)
        raise


def get_protostructure_label_from_moyopy(
    struct: Structure,
    raise_errors: bool = False,
    symprec: float = 0.1,
) -> str | None:
    """Get AFLOW prototype label using Moyopy for symmetry detection.

    Args:
        struct (Structure): pymatgen Structure object.
        raise_errors (bool): Whether to raise errors or annotate them. Defaults to
            False.
        symprec (float): Initial symmetry precision for Moyopy. Defaults to 0.1.

    Returns:
        str: protostructure_label which is constructed as `aflow_label:chemsys` or
            explanation of failure if symmetry detection failed and `raise_errors`
            is False.
    """
    if not HAS_MOYOPY:
        raise ImportError("moyopy not found, run pip install moyopy")

    moyo_cell = MoyoAdapter.from_structure(struct)
    moyo_data = moyopy.MoyoDataset(moyo_cell, symprec=symprec)

    spg_num = moyo_data.number
    pearson_symbol = moyo_data.pearson_symbol
    prototype_form = get_prototype_formula_from_composition(struct.composition)
    chemsys = struct.chemical_system

    equivalent_wyckoff_labels = []
    orbit_groups: dict[int, list[int]] = {}

    for idx, orbit_id in enumerate(moyo_data.orbits):
        if orbit_id not in orbit_groups:
            orbit_groups[orbit_id] = []
        orbit_groups[orbit_id].append(idx)

    for orbit in orbit_groups.values():
        wyckoff = moyo_data.wyckoffs[orbit[0]]
        element = struct.species[orbit[0]]
        equivalent_wyckoff_labels += [(len(orbit), element.symbol, wyckoff.translate(remove_digits))]

    all_wyckoffs, element_dict = _get_all_wyckoffs_substring_and_element_dict(equivalent_wyckoff_labels, spg_num)

    protostructure_label = f"{prototype_form}_{pearson_symbol}_{spg_num}_{all_wyckoffs}:{chemsys}"

    observed_formula = Composition(element_dict).reduced_formula
    expected_formula = struct.composition.reduced_formula
    if observed_formula != expected_formula:
        err_msg = (
            f"Invalid WP multiplicities - {protostructure_label}, expected {observed_formula} to be {expected_formula}"
        )
        if raise_errors:
            raise ValueError(err_msg)
        return err_msg

    return protostructure_label


def count_distinct_wyckoff_letters(protostructure_label: str) -> int:
    """Count number of distinct Wyckoff letters in protostructure_label.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: number of distinct Wyckoff letters in protostructure_label
    """
    aflow_label, _ = protostructure_label.split(":")
    _, _, _, element_wyckoffs = aflow_label.split("_", 3)
    element_wyckoffs = element_wyckoffs.translate(remove_digits).replace("_", "")
    return len(set(element_wyckoffs))


def count_wyckoff_positions(protostructure_label: str) -> int:
    """Count number of Wyckoff positions in protostructure_label.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: number of distinct Wyckoff positions in protostructure_label
    """
    aflow_label, _ = protostructure_label.split(":")
    wyk_letters = aflow_label.split("_", maxsplit=3)[-1]
    wyk_letters = wyk_letters.replace("_", "")
    wyk_list = re.split("[A-z]", wyk_letters)[:-1]

    return sum(1 if len(x) == 0 else int(x) for x in wyk_list)


def count_crystal_dof(protostructure_label: str) -> int:
    """Count number of free parameters in coarse-grained protostructure_label
    representation: how many degrees of freedom would remain to optimize during
    a crystal structure relaxation.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: Number of free-parameters in given prototype
    """
    aflow_label, _ = protostructure_label.split(":")
    _, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    return (
        _count_from_dict(element_wyckoffs, WYCKOFF_POSITION_PARAM_DICT, spg_num)
        + CRYSTAL_LATTICE_PARAMETERS_COUNTS[pearson_symbol[0]]
    )


def count_crystal_sites(protostructure_label: str) -> int:
    """Count number of sites from protostructure_label.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        int: Number of free-parameters in given prototype
    """
    aflow_label, _ = protostructure_label.split(":")
    _, _, spg_num, *element_wyckoffs = aflow_label.split("_")

    return _count_from_dict(element_wyckoffs, WYCKOFF_MULTIPLICITY_DICT, spg_num)


def _count_from_dict(element_wyckoffs: list[str], lookup_dict: dict, spg_num: str) -> int:
    """Count number of sites from protostructure_label."""
    n_params = 0

    for wyckoffs in element_wyckoffs:
        sep_el_wyks = split_alpha_numeric(re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, wyckoffs))
        n_params += count_values_for_wyckoff(
            sep_el_wyks["alpha"],
            sep_el_wyks["numeric"],
            spg_num,
            lookup_dict,
        )

    return int(n_params)


def get_prototype_from_protostructure(protostructure_label: str) -> str:
    """Get a canonicalized string for the prototype. This prototype should be
    the same for all isopointal protostructures.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.

    Returns:
        str: Canonicalized AFLOW-style prototype label
    """
    aflow_label, _ = protostructure_label.split(":")
    prototype_formula, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    anonymous_formula = get_anonymous_formula_from_prototype_formula(prototype_formula)
    counts = [
        int(x)
        for x in split_alpha_numeric(re.sub(RE_ELEMENT_NO_SUFFIX, RE_SUBST_ONE_SUFFIX, prototype_formula))["numeric"]
    ]

    counts, element_wyckoffs = map(list, zip(*sorted(zip(counts, element_wyckoffs, strict=True)), strict=True))
    all_wyckoffs = "_".join(element_wyckoffs)
    all_wyckoffs = re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, all_wyckoffs)
    if len(counts) == len(set(counts)):
        all_wyckoffs = canonicalize_element_wyckoffs(all_wyckoffs, int(spg_num))
        return f"{anonymous_formula}_{pearson_symbol}_{spg_num}_{all_wyckoffs}"

    all_wyckoffs_permutations = [
        "_".join(list(map(itemgetter(1), chain.from_iterable(p))))
        for p in product(
            *[
                permutations(g)
                for _, g in groupby(sorted(zip(counts, all_wyckoffs.split("_"), strict=True)), key=lambda x: x[0])
            ]
        )
    ]

    isopointal_all_wyckoffs = list(
        {
            all_wyckoffs.translate(str.maketrans(trans))
            for all_wyckoffs in all_wyckoffs_permutations
            for trans in WYCKOFF_POSITION_RELAB_DICT[spg_num]
        }
    )

    scored_all_wyckoffs = [
        sort_and_score_element_wyckoffs(element_wyckoffs) for element_wyckoffs in isopointal_all_wyckoffs
    ]

    all_wyckoffs = min(scored_all_wyckoffs, key=lambda x: (x[1], x[0]))[0]

    return f"{anonymous_formula}_{pearson_symbol}_{spg_num}_{all_wyckoffs}"


def _find_translations(dict1: dict[str, int], dict2: dict[str, int]) -> list[dict[str, str]]:
    """Find all possible translations between two dictionaries."""
    if Counter(dict1.values()) != Counter(dict2.values()):
        return []

    keys2 = list(dict2.keys())
    used = set()

    def backtrack(translation, index):
        if index == len(dict1):
            return [translation.copy()]

        key1 = list(dict1.keys())[index]
        value1 = dict1[key1]
        valid_translations = []

        for key2 in keys2:
            if key2 not in used and dict2[key2] == value1:
                used.add(key2)
                translation[key1] = key2
                valid_translations.extend(backtrack(translation, index + 1))
                used.remove(key2)
                del translation[key1]

        return valid_translations

    return backtrack({}, 0)


def get_protostructures_from_aflow_label_and_composition(aflow_label: str, composition: Composition) -> list[str]:
    """Get a canonicalized string for the prototype.

    Args:
        aflow_label (str): AFLOW-style prototype label
        composition (Composition): pymatgen Composition object

    Returns:
        list[str]: List of possible protostructure labels that can be generated
            from combinations of the input aflow_label and composition.
    """
    anonymous_formula, pearson_symbol, spg_num, *element_wyckoffs = aflow_label.split("_")

    ele_amt_dict = composition.get_el_amt_dict()
    proto_formula = get_prototype_formula_from_composition(composition)
    anom_amt_dict = _get_anonymous_formula_dict(anonymous_formula)

    translations = _find_translations(ele_amt_dict, anom_amt_dict)
    anom_ele_to_wyk = dict(zip(anom_amt_dict.keys(), element_wyckoffs, strict=True))
    anonymous_formula = RE_ANONYMOUS.sub(RE_SUBST_ONE_PREFIX, anonymous_formula)

    protostructures = set()
    for t in translations:
        wyckoff_part = "_".join(
            RE_WYCKOFF.sub(RE_SUBST_ONE_PREFIX, anom_ele_to_wyk[t[elem]]) for elem in sorted(t.keys())
        )
        canonicalized_wyckoff = canonicalize_element_wyckoffs(wyckoff_part, spg_num)
        chemical_system = "-".join(sorted(t.keys()))

        protostructures.add(f"{proto_formula}_{pearson_symbol}_{spg_num}_{canonicalized_wyckoff}:{chemical_system}")

    return list(protostructures)


def get_random_structure_for_protostructure(protostructure_label: str, **kwargs) -> Structure:
    """Generate a random structure for a given prototype structure.

    NOTE that due to the random nature of the generation, the output structure
    may be higher symmetry than the requested prototype structure.

    Args:
        protostructure_label (str): label constructed as `aflow_label:chemsys` where
            aflow_label is an AFLOW-style prototype label chemsys is the alphabetically
            sorted chemical system.
        **kwargs: Keyword arguments to pass to pyxtal().from_random()
    """
    if not HAS_PYXTAL:
        raise ImportError("pyxtal is required for this function")

    aflow_label, chemsys = protostructure_label.split(":")
    _, _, spg_num, *element_wyckoffs = aflow_label.split("_")

    sep_el_wyks = [split_alpha_numeric(re.sub(RE_WYCKOFF_NO_PREFIX, RE_SUBST_ONE_PREFIX, w)) for w in element_wyckoffs]

    species_sites = [
        [
            site
            for count, wyckoff_letter in zip(d["numeric"], d["alpha"], strict=True)
            for site in [f"{WYCKOFF_MULTIPLICITY_DICT[spg_num][wyckoff_letter]}{wyckoff_letter}"] * int(count)
        ]
        for d in sep_el_wyks
    ]

    species_counts = [
        sum(
            WYCKOFF_MULTIPLICITY_DICT[spg_num][wyckoff_letter] * int(count)
            for count, wyckoff_letter in zip(d["numeric"], d["alpha"], strict=True)
        )
        for d in sep_el_wyks
    ]

    p = pyxtal()
    p.from_random(
        dim=3,
        group=int(spg_num),
        species=chemsys.split("-"),
        numIons=species_counts,
        sites=species_sites,
        **kwargs,
    )
    return p.to_pymatgen()
