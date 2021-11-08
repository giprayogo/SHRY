# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name,missing-function-docstring,wrong-import-order,unused-import

"""Test core operations."""

import functools
import glob
import io
import os
import pickle
import shutil
import subprocess
from pprint import pprint

import numpy as np
import pandas as pd
import pytest
import sympy
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from shry.core import PatternMaker, Polya, Substitutor
from shry.main import LabeledStructure, ScriptHelper
from sympy.tensor.indexed import IndexedBase

from helper import chdir

# Test fixtures


@pytest.fixture
def run():
    structure_file = np.array
    structure_file = "canaug/PbSnTe2.cif"
    scaling_matrix = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    return ScriptHelper(structure_file=structure_file, scaling_matrix=scaling_matrix)


@pytest.fixture
@chdir("canaug")
def pattern_generator(request):
    task = ScriptHelper.from_file("config.ini")
    sites = task.substitutor.disorder_sites["0"]
    permutations = task.substitutor.filter_permutation(sites)

    symbols = np.unique(permutations)
    if request.param:
        dmat = task.substitutor.dmats["0"]
        return PatternMaker(permutations, symbols, dmat)
    return PatternMaker(permutations, symbols)


@pytest.fixture
def generator():
    symbols = np.array([8, 9, "a", 11])
    permutations = np.array(
        [
            [8, 9, "a", 11],
            [9, "a", 8, 11],
            ["a", 8, 9, 11],
        ]
    )
    return PatternMaker(permutations, symbols)


def test_signature():
    """Test "signature" of permutations.

    So that we can reuse the generator for the same symmetry!
    p.s. turn out this is not trivial...
    """
    # Invariance with row swaps.
    perm_a = np.array(
        [
            [8, 9, "a", 11],
            [9, "a", 8, 11],
            ["a", 8, 9, 11],
        ]
    )
    perm_b = np.array(
        [
            [8, 9, "a", 11],
            ["a", 8, 9, 11],
            [9, "a", 8, 11],
        ]
    )

    # Invariance with symbol change
    perm_c = np.array(
        [
            [9, 8, "a", 11],
            [8, "a", 9, 11],
            ["a", 9, 8, 11],
        ]
    )
    perm_d = np.array(
        [
            [21, "x", "a", 0],
            ["x", "a", 21, 0],
            ["a", 21, "x", 0],
        ]
    )

    # Technically the same permutation though at different site.
    # (map needs to be more flexible)
    perm_e = np.array(
        [
            [11, 8, 9, "a"],
            [11, 9, "a", 8],
            [11, "a", 8, 9],
        ]
    )
    # Symbol change + row change
    perm_f = np.array(
        [
            [8, "a", 9, 11],
            [9, 8, "a", 11],
            ["a", 9, 8, 11],
        ]
    )

    pgs = [PatternMaker(x) for x in (perm_a, perm_b, perm_c, perm_d, perm_e, perm_f)]
    assert all(pg.label == pgs[0].label for pg in pgs)
    assert all(
        pg.label == PatternMaker.get_label(x)
        for pg, x in zip(pgs, (perm_a, perm_b, perm_c, perm_d, perm_e, perm_f))
    )


@chdir("canaug")
def test_symm_mat():
    """Check the correctness of the generated symmetry matrix."""
    task = ScriptHelper.from_file("split_orbit.ini")

    # Needs correct split orbit.
    mat = np.load("raw_symm.npy")
    print(mat[np.lexsort(mat.T)[::-1]])
    print(task.substitutor.permutations())
    assert np.array_equal(mat, task.substitutor.permutations())


def test_desymbolized_permutations(generator):
    norm_p = generator.desymbolized_permutations()
    assert (
        norm_p
        == np.array(
            [
                [1, 2, 0, 3],
                [2, 0, 1, 3],
                [0, 1, 2, 3],
            ]
        )
    ).all()
    print(norm_p.min(axis=0))
    assert norm_p.dtype == "int64"


def test_exp_permutations(generator):
    exp_p = generator.exp_permutations()
    print(exp_p)
    assert (
        exp_p
        == np.array(
            [
                [2, 4, 1, 8],
                [4, 1, 2, 8],
                [1, 2, 4, 8],
            ]
        )
    ).all()
    assert exp_p.dtype == "float64"


def test_search_noninvar(generator):
    correct_answer = {
        0: [set()],
        1: [{"11"}, {"8"}],
        2: [{"8", "11"}, {"8", "9"}],
        3: [{"11", "9", "a"}, {"8", "9", "a"}],
        4: [{"8", "9", "a", "11"}],
    }
    # Test: should not intefere
    generator.get_pattern(1)
    print(generator.get_patterns())
    for n, patterns in generator.get_patterns().items():
        for pattern in patterns:
            assert set(pattern) in correct_answer[n]


@pytest.mark.parametrize(
    "target_compositions",
    [
        "Pb0.5Sn0.5",
        "Pb0.5Sn0.125Ge0.375",
        "Pb2+0.375Pb0.125Pb3+0.5",
        "Pb0.25Fl0.25Sn0.25Ge0.25",
        "Pb0.25Fl0.25Sn0.125Ge0.375",
    ],
)
@chdir("canaug")
def test_multi_substitution(target_compositions):
    """Verify multi-substitution tree with Polya enumeration."""
    task = ScriptHelper(
        structure_file="PbSnTe2.cif",
        from_species="Pb",
        to_species=target_compositions,
        scaling_matrix="1 1 2",
    )
    patterns = task.substitutor.patterns()
    n_patterns = [len(x) for x in patterns.values()]
    print(f"Patterns: {patterns}")
    print(f"Pattern counts: {n_patterns}")
    counts = task.substitutor.counts()
    n_counts = list(counts.values())
    print(f"Counts: {counts}")

    assert n_patterns == n_counts


@chdir("canaug")
def test_get_pattern():
    """Temporary name."""
    task = ScriptHelper(
        structure_file="PbSnTe2.cif",
        from_species="Pb",
        to_species="Pb0.5Sn0.125Ge0.375",
        scaling_matrix="1 2 2",
    )
    print(task.substitutor._pattern)
    # print(task.analyzer._pattern("0"))


@chdir("canaug")
def test_need_supercell():
    task = ScriptHelper(
        structure_file="2_2_2_ICSD_CollCode430061.cif",
        from_species="Ti",
        to_species="Ti0.51Bi0.49",
        scaling_matrix="2 0 0 0 1 0,0 0,3",
    )
    assert task.write() is None
    assert task.count() is None


@pytest.mark.parametrize(
    "structure_file, supercell, substitute, to, answer",
    [
        ("01aP_ICSD_CollCode20537.cif", [3, 1, 1], "V", "V4+2Nb4+", 90),
        ("01aP_ICSD_CollCode20537.cif", [3, 1, 1], "V1,V2", "V4+2Nb4+,V4+2Nb4+", 42),
        ("01aP_ICSD_CollCode20537-1-1-1.cif", [3, 1, 1], "", "", 42),
        ("05oA_ICSD_CollCode67631.cif", [1, 1, 4], "Ba", "Sr2+Ba2+3", 130),
        ("05oA_ICSD_CollCode67631-1-1-1.cif", [1, 4, 1], "", "", 130),
        ("12cP_ICSD_CollCode14000.cif", [2, 1, 1], "Ag", "AgAu", 72),
        ("12cP_ICSD_CollCode14000-1-1-1.cif", [2, 1, 1], "", "", 72),
        # ("14cF_ICSD_CollCode40966.cif", "1 1 3", "Na", "NaKRbCs"),
        ("14cF_ICSD_CollCode40966.cif", [1, 1, 3], "Na", "NaKRb2", 425),
        ("PbSnTe2.cif", [2, 2, 2], ["Pb"], ["Pb6CuSn"], 499129),
    ],
)
# @chdir("structures")
@chdir("../examples")
def test_edges(structure_file, supercell, substitute, to, answer):
    # TODO: Separate into individual tests.
    task = ScriptHelper(
        structure_file=structure_file,
        from_species=substitute,
        to_species=to,
        scaling_matrix=supercell,
    )
    # patterns = task.substitutor.patterns()
    task.substitutor.make_patterns()
    weights = task.substitutor.weights()
    letters = task.substitutor.configurations()
    counts = task.substitutor.count()

    # They need to be in correct dimension
    # patterns_shapes = set(x.shape for x in patterns.values())
    # assert len(patterns_shapes) == 1
    # assert patterns_shapes.pop()[0] == answer

    # TODO: I changed the format
    # letters_shapes = set(len(x) for x in letters.values())
    # assert len(letters_shapes) == 1
    # assert letters_shapes.pop() == answer

    assert len(weights) == answer
    assert counts == answer

    for cifwriter in task.substitutor.cifwriters(symprec=1e-5):
        pass


@pytest.mark.parametrize(
    "structure_file, supercell",
    [
        ("01aP_ICSD_CollCode20537-1-1-1.cif", [3, 1, 1]),
        ("02mP_ICSD_CollCode3923-1-1-1.cif", [1, 3, 1]),
        ("03mC_ICSD_CollCode484-1-1-1.cif", [1, 2, 2]),
        ("04oP_ICSD_CollCode2208-1-1-1.cif", [1, 1, 3]),
        ("05oA_ICSD_CollCode67631-1-1-1.cif", [1, 4, 1]),
        ("06ol_ICSD_CollCode8132-1-1-1.cif", [1, 1, 1]),
        ("07oF_ICSD_CollCode2521-1-1-1.cif", [1, 1, 1]),
        ("08tP_ICSD_CollCode430061-1-1-1.cif", [2, 2, 2]),
        ("09tl_ICSD_CollCode603237-1-1-1.cif", [2, 2, 1]),
        ("10hP_ICSD_CollCode251-1-1-1.cif", [2, 1, 1]),
        ("11hR_ICSD_CollCode2085-1-1-1.cif", [1, 1, 2]),
        ("12cP_ICSD_CollCode14000-1-1-1.cif", [2, 1, 1]),
        ("13cI_ICSD_CollCode1286-1-1-1.cif", [1, 1, 1]),
        ("14cF_ICSD_CollCode40966_2-1-1-1.cif", [1, 1, 3]),
    ],
)
@chdir("structures")
def test_bravais(structure_file, supercell):
    """
    Consistency check with supercell for 14 bravais lattices
    """

    def read_log(path):
        """
        Parse inspection log into dataframe
        """
        log = pd.read_csv(path, delim_whitespace=True)

        # capital -> small
        to_small = lambda str: str.lower().replace("_", "")
        log["GroupName"] = log["GroupName"].map(to_small)

        # Sort
        log = log.sort_values(
            ["Weight", "GroupName"],
            ascending=[True, True],
            kind="mergesort",  # Needs to be a stable sort.
        )
        log.reset_index(drop=True, inplace=True)
        return log

    def cleanup():
        for shry_output_dirs in glob.glob("shry*"):
            if os.path.isdir(shry_output_dirs):
                shutil.rmtree(shry_output_dirs)
        for supercell_cif in glob.glob("supercell_*.cif"):
            os.remove(supercell_cif)

    # SHRY
    # Clear previous artefacts
    cleanup()

    try:
        helper = ScriptHelper(
            structure_file=structure_file, scaling_matrix=supercell, write_symm=True
        )
        # helper.write_log()
        helper.write()

        # Preserve log
        shry_logs = glob.glob("shry*/sub.log")
        shutil.copy(shry_logs[0], "shry.log")

        # Supercell
        # Make a log in the same format
        logio = io.StringIO()
        print("N Weight Configuration GroupName", file=logio)

        supercell_string = "x".join(map(str, supercell))
        subprocess.run(
            ["supercell", "-m", "-i", structure_file, "-s", supercell_string],
            check=True,
        )

        # Build the log
        for output_cif in glob.glob("supercell_*.cif"):
            i = str(int(output_cif.split("_")[1].strip("i")))
            configuration = "n/a"
            weight = output_cif.rstrip(".cif").split("_")[-1].lstrip("w")

            # Get symmetry information from pymatgen
            structure = Structure.from_file(output_cif)
            space_group = SpacegroupAnalyzer(structure).get_space_group_symbol()
            print(i, weight, configuration, space_group, file=logio)

        # Write log
        with open("supercell.log", "w") as f:
            logio.seek(0)
            shutil.copyfileobj(logio, f)

        # Compare
        shry_log = read_log("shry.log")
        supercell_log = read_log("supercell.log")

        # Check consistency
        # Check length (number of irreducible structure)
        assert len(shry_log) == len(supercell_log)
        # Check weight
        assert shry_log["Weight"].sum() == supercell_log["Weight"].sum()

        # Check consistency in detail
        assert (shry_log["Weight"] == supercell_log["Weight"]).all()
        nottogether = [
            i
            for i, xy in enumerate(
                zip(shry_log["GroupName"], supercell_log["GroupName"])
            )
            if xy[0] != xy[1]
        ]
        print(">>>>", nottogether)
        for i in nottogether:
            # copy up
            ...

        assert (shry_log["GroupName"] == supercell_log["GroupName"]).all()
    finally:
        # Clear artefacts
        cleanup()


@chdir("../examples")
def test_sequential():
    """
    Test sequential use of Substitutor
    """
    structure = LabeledStructure.from_file("SmFe12.cif")
    structure *= [1, 2, 1]
    structure1 = structure.copy()
    structure2 = structure.copy()
    structure1.replace_species({"Fe1": "Fe7Ti"})
    structure2.replace_species({"Fe2": "Fe6Ti2"})

    substitutor = Substitutor(structure1)
    assert len(substitutor.weights()) == 15
    substitutor.structure = structure2
    assert len(substitutor.weights()) == 147


@pytest.mark.parametrize(
    "structure_file, supercell, substitute, to",
    [
        ("01aP_ICSD_CollCode20537.cif", "3 1 1", "V", "V4+2Nb4+"),
        ("14cF_ICSD_CollCode40966.cif", "1 1 3", "Na", "NaKRb2"),
    ],
)
@chdir("bravais")
def test_order_structure(structure_file, supercell, substitute, to):
    """Test pymatgen.Structure generator."""
    task = ScriptHelper(
        structure_file=structure_file,
        from_species=substitute,
        to_species=to,
        scaling_matrix=supercell,
    )
    task.substitutor.patterns()

    # Generator should generate.
    g = task.substitutor.order_structure()
    result_g = list(x for x in g)
    assert result_g[0] != result_g[1]

    # Generator should restart on a new call
    h = task.substitutor.order_structure()
    assert next(h) == result_g[0]

    # The length should match count.
    assert len(result_g) == task.substitutor.count()


# Tests related to Polya


@pytest.fixture
def polya():
    """Fixture returning a Polya instance."""
    perm_a = np.array(
        [
            [8, 9, 10, 11],
            [9, 10, 8, 11],
            [10, 8, 9, 11],
        ]
    )
    perm_b = np.array(
        [
            [3, 4, 5],
            [5, 4, 3],
            [3, 4, 5],
        ]
    )
    perms_list = [perm_a, perm_b]
    return Polya(perms_list)


def test_ci(polya):
    """Test cycle index calculation."""
    a = IndexedBase("a")
    b = IndexedBase("b")
    assert polya.ci() == {
        0: a[1] ** 4 * b[1] ** 3,
        1: a[1] * a[3] * b[1] * b[2],
        2: a[1] * a[3] * b[1] ** 3,
    }


def test_cgf(polya):
    """Test configuration generation function calculation."""
    cgf = polya.cgf([2, 3])
    with open("cgf.pkl", "rb") as f:
        answer = pickle.load(f)
    assert cgf == answer


def test_count(polya):
    """Test counting of pattern. One should be enough representative."""
    assert polya.count([[3, 1], [2, 1]]) == 5


@chdir("../examples")
def test_no_disorder():
    """Test behaviour when the supplied Structure has no disorder sites."""
    structure = LabeledStructure.from_file("SmFe12.cif")
    substitutor = Substitutor(structure)
    assert substitutor.patterns() == dict()
    assert substitutor.configurations() == dict()
    assert substitutor.weights() == []
    assert substitutor.count() == 0
