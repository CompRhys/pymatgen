from __future__ import annotations

import gzip
import json
from typing import TYPE_CHECKING

import numpy as np
import pytest
from monty.json import MontyDecoder, MontyEncoder
from numpy.testing import assert_allclose

from pymatgen.core import SETTINGS, Composition, Lattice, Species, Structure
from pymatgen.io.aims.inputs import (
    ALLOWED_AIMS_CUBE_TYPES,
    ALLOWED_AIMS_CUBE_TYPES_STATE,
    AimsControlIn,
    AimsCube,
    AimsGeometryIn,
    AimsSpeciesFile,
    SpeciesDefaults,
)
from pymatgen.util.testing import TEST_FILES_DIR

from .conftest import compare_single_files as compare_files

if TYPE_CHECKING:
    from pathlib import Path

TEST_DIR = TEST_FILES_DIR / "io/aims/input_files"


def test_read_write_si_in(tmp_path: Path):
    si = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.si.gz")

    in_lattice = np.array([[0, 2.715, 2.716], [2.717, 0, 2.718], [2.719, 2.720, 0]])
    in_coords = np.array([[0, 0, 0], [0.25, 0.24, 0.26]])

    assert all(sp.symbol == "Si" for sp in si.structure.species)
    assert_allclose(si.structure.lattice.matrix, in_lattice)
    assert_allclose(si.structure.frac_coords.flatten(), in_coords.flatten())

    si_test_from_struct = AimsGeometryIn.from_structure(si.structure)
    assert si.structure == si_test_from_struct.structure

    si_test_from_struct.write_file(directory=tmp_path, overwrite=True)
    with pytest.raises(ValueError, match="geometry.in file exists in "):
        si_test_from_struct.write_file(directory=tmp_path, overwrite=False)

    compare_files(TEST_DIR / "geometry.in.si.ref", f"{tmp_path}/geometry.in")

    si.structure.to(tmp_path / "si.in", fmt="aims")
    compare_files(TEST_DIR / "geometry.in.si.ref", f"{tmp_path}/si.in")

    si_from_file = Structure.from_file(f"{tmp_path}/geometry.in")
    assert all(sp.symbol == "Si" for sp in si_from_file.species)

    with gzip.open(f"{TEST_DIR}/si_ref.json.gz", mode="rt") as si_ref_json:
        si_from_dct = json.load(si_ref_json, cls=MontyDecoder)

    assert si.structure == si_from_dct.structure


def test_read_h2o_in(tmp_path: Path):
    h2o = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.h2o.gz")

    in_coords = [
        [0, 0, 0.119262],
        [0, 0.763239, -0.477047],
        [0, -0.763239, -0.477047],
    ]

    assert all(sp.symbol == symb for sp, symb in zip(h2o.structure.species, ["O", "H", "H"], strict=True))
    assert_allclose(h2o.structure.cart_coords, in_coords)

    h2o_test_from_struct = AimsGeometryIn.from_structure(h2o.structure)
    assert h2o.structure == h2o_test_from_struct.structure

    h2o_test_from_struct.write_file(directory=tmp_path, overwrite=True)

    with pytest.raises(ValueError, match="geometry.in file exists in "):
        h2o_test_from_struct.write_file(directory=tmp_path, overwrite=False)

    compare_files(TEST_DIR / "geometry.in.h2o.ref", f"{tmp_path}/geometry.in")

    with gzip.open(f"{TEST_DIR}/h2o_ref.json.gz", mode="rt") as h2o_ref_json:
        h2o_from_dct = json.load(h2o_ref_json, cls=MontyDecoder)

    assert h2o.structure == h2o_from_dct.structure


def test_write_spins(tmp_path: Path):
    mg2mn4o8 = Structure(
        lattice=Lattice(
            [
                [5.06882343, 0.00012488, -2.66110167],
                [-1.39704234, 4.87249911, -2.66110203],
                [0.00986091, 0.01308528, 6.17649359],
            ]
        ),
        species=[
            "Mg",
            "Mg",
            Species("Mn", spin=5.0),
            Species("Mn", spin=5.0),
            Species("Mn", spin=5.0),
            Species("Mn", spin=5.0),
            "O",
            "O",
            "O",
            "O",
            "O",
            "O",
            "O",
            "O",
        ],
        coords=[
            [0.37489726, 0.62510274, 0.75000002],
            [0.62510274, 0.37489726, 0.24999998],
            [-0.00000000, -0.00000000, 0.50000000],
            [-0.00000000, 0.50000000, 0.00000000],
            [0.50000000, -0.00000000, 0.50000000],
            [-0.00000000, -0.00000000, 0.00000000],
            [0.75402309, 0.77826750, 0.50805882],
            [0.77020285, 0.24594779, 0.99191316],
            [0.22173254, 0.24597689, 0.99194116],
            [0.24597691, 0.22173250, 0.49194118],
            [0.24594765, 0.77020288, 0.49191313],
            [0.22979715, 0.75405221, 0.00808684],
            [0.75405235, 0.22979712, 0.50808687],
            [0.77826746, 0.75402311, 0.00805884],
        ],
    )

    geo_in = AimsGeometryIn.from_structure(mg2mn4o8)

    magmom_lines = [line for line in geo_in.content.split("\n") if "initial_moment" in line]
    assert len(magmom_lines) == 4

    magmoms = np.array([float(line.strip().split()[-1]) for line in magmom_lines])
    assert_allclose(magmoms, 5.0)

    mg2mn4o8 = Structure(
        lattice=mg2mn4o8.lattice,
        species=mg2mn4o8.species,
        coords=mg2mn4o8.frac_coords,
        site_properties={"magmom": np.zeros(mg2mn4o8.num_sites)},
    )
    with pytest.raises(
        ValueError,
        match="species.spin and magnetic moments don't agree. Please only define one",
    ):
        geo_in = AimsGeometryIn.from_structure(mg2mn4o8)


def check_wrong_type_aims_cube(cube_type, exp_err):
    with pytest.raises(ValueError, match=exp_err):
        AimsCube(type=cube_type)


def test_aims_cube():
    check_wrong_type_aims_cube(cube_type="INCORRECT_TYPE", exp_err="Cube type undefined")

    for cube_type in ALLOWED_AIMS_CUBE_TYPES_STATE:
        check_wrong_type_aims_cube(
            cube_type=cube_type,
            exp_err=f"{cube_type=} must have a state associated with it",
        )

    for cube_type in ALLOWED_AIMS_CUBE_TYPES:
        check_wrong_type_aims_cube(
            cube_type=f"{cube_type} 1",
            exp_err=f"{cube_type=} can not have a state associated with it",
        )

    with pytest.raises(
        ValueError,
        match="TEST_ERR is invalid. Cube files must have a format of",
    ):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], format="TEST_ERR")

    with pytest.raises(ValueError, match=r"Spin state must be one of \(1, 2, None\)"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], spin_state=3)

    with pytest.raises(ValueError, match="The cube origin must have 3 components"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], origin=[0])

    with pytest.raises(ValueError, match="Only three cube edges can be passed"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], edges=[[0, 0, 0.1]])

    with pytest.raises(ValueError, match="Each cube edge must have 3 components"):
        AimsCube(
            type=ALLOWED_AIMS_CUBE_TYPES[0],
            edges=[[0, 0, 0.1], [0.1, 0, 0], [0.1, 0]],
        )

    with pytest.raises(
        ValueError,
        match="elf_type is only used when the cube type is elf. Otherwise it must be None",
    ):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], elf_type=1)

    with pytest.raises(ValueError, match="The number of points per edge must have 3 components"):
        AimsCube(type=ALLOWED_AIMS_CUBE_TYPES[0], points=[100, 100, 100, 100])

    test_cube = AimsCube(
        type="elf",
        origin=[0, 0, 0],
        edges=[[0.01, 0, 0], [0, 0.01, 0], [0, 0, 0.01]],
        points=[100, 100, 100],
        spin_state=1,
        kpoint=1,
        filename="test.cube",
        format="cube",
        elf_type=1,
    )

    test_cube_block = [
        "output cube elf",
        "    cube origin  0.000000000000e+00  0.000000000000e+00  0.000000000000e+00",
        "    cube edge 100  1.000000000000e-02  0.000000000000e+00  0.000000000000e+00",
        "    cube edge 100  0.000000000000e+00  1.000000000000e-02  0.000000000000e+00",
        "    cube edge 100  0.000000000000e+00  0.000000000000e+00  1.000000000000e-02",
        "    cube format cube",
        "    cube spinstate 1",
        "    cube kpoint 1",
        "    cube filename test.cube",
        "    cube elf_type 1",
        "",
    ]
    assert test_cube.control_block == "\n".join(test_cube_block)

    test_cube_from_dict = json.loads(json.dumps(test_cube.as_dict(), cls=MontyEncoder), cls=MontyDecoder)
    assert test_cube_from_dict.control_block == test_cube.control_block


def test_aims_control_in(tmp_path: Path):
    parameters = {
        "cubes": [
            AimsCube(type="eigenstate 1", points=[10, 10, 10]),
            AimsCube(type="total_density", points=[10, 10, 10]),
        ],
        "xc": "LDA",
        "smearing": ["fermi-dirac", 0.01],
        "vdw_correction_hirshfeld": True,
        "compute_forces": True,
        "relax_geometry": ["trm", "1e-3"],
        "batch_size_limit": 200,
        "species_dir": "light",
    }

    aims_control = AimsControlIn(parameters.copy())

    for key, val in parameters.items():
        assert aims_control[key] == val

    del aims_control["xc"]
    assert "xc" not in aims_control.parameters
    aims_control.parameters = parameters

    h2o = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.h2o.gz").structure

    si = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.si.gz").structure
    aims_control.write_file(h2o, directory=tmp_path, overwrite=True)

    compare_files(TEST_DIR / "control.in.h2o", f"{tmp_path}/control.in")

    with pytest.raises(ValueError, match="k-grid must be defined for periodic systems"):
        aims_control.write_file(si, directory=tmp_path, overwrite=True)
    aims_control["k_grid"] = [1, 1, 1]
    aims_control["xc"] = "libxc LDA_X+LDA_C_PW"

    with pytest.raises(ValueError, match="control.in file already in "):
        aims_control.write_file(si, directory=tmp_path, overwrite=False)

    aims_control["output"] = "band 0 0 0 0.5 0 0.5 10 G X"
    aims_control["output"] = "band 0 0 0 0.5 0.5 0.5 10 G L"

    aims_control_from_dict = json.loads(json.dumps(aims_control.as_dict(), cls=MontyEncoder), cls=MontyDecoder)
    for key, val in aims_control.parameters.items():
        if key in ["output", "cubes"]:
            np.all(aims_control_from_dict[key] == val)
        assert aims_control_from_dict[key] == val

    aims_control_from_dict.write_file(si, directory=tmp_path, verbose_header=True, overwrite=True)
    compare_files(TEST_DIR / "control.in.si", f"{tmp_path}/control.in")


def test_species_file(monkeypatch: pytest.MonkeyPatch):
    """Tests an AimsSpeciesFile class"""
    monkeypatch.setitem(SETTINGS, "AIMS_SPECIES_DIR", str(TEST_DIR.parent / "species_directory"))
    species_file = AimsSpeciesFile.from_element_and_basis_name("Si", "light", label="Si_surface")
    assert species_file.label == "Si_surface"
    assert species_file.element == "Si"


def test_species_defaults(monkeypatch: pytest.MonkeyPatch):
    """Tests an AimsSpeciesDefaults class"""
    monkeypatch.setitem(SETTINGS, "AIMS_SPECIES_DIR", str(TEST_DIR.parent / "species_directory"))
    si = AimsGeometryIn.from_file(TEST_DIR / "geometry.in.si.gz").structure
    species_defaults = SpeciesDefaults.from_structure(si, "light")
    assert species_defaults.labels == [
        "Si",
    ]
    assert species_defaults.elements == {"Si": "Si"}

    si[0].species = Composition({"Si0+": 1}, strict=True)
    species_defaults = SpeciesDefaults.from_structure(si, "light")
    assert species_defaults.labels == ["Si", "Si0+"]
    assert species_defaults.elements == {"Si": "Si", "Si0+": "Si"}
    assert "Si0+" in str(species_defaults)
    assert "Si" in str(species_defaults)
