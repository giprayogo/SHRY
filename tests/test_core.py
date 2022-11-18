# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name,missing-function-docstring,wrong-import-order,unused-import,invalid-name,protected-access
"""Test core operations."""

#python modules
import os
import filecmp
import glob
import shutil
import numpy as np
import pandas as pd
import pytest
from sympy.tensor.indexed import IndexedBase

#pymatgen
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core import Structure

#shry
from shry.core import NeedSupercellError, PatternMaker, Polya, Substitutor, TooBigError
from shry.main import LabeledStructure, ScriptHelper
from helper import chdir

# Tolerances
SHRY_TOLERANCE=0.01      #angstrom
SHRY_ANGLE_TOLERANCE=5.0 #degree

# PatternMaker basic functions.

def test_perm_label():
    """
    Test rough canonization attempt on the permutations.
    This allows reuse of patterns on multiple color/orbit.
    """
    # Invariance with row swaps.
    perm_a = np.array([[8, 9, "a", 11], [9, "a", 8, 11], ["a", 8, 9, 11],])
    perm_b = np.array([[8, 9, "a", 11], ["a", 8, 9, 11], [9, "a", 8, 11],])

    # Invariance with symbol change
    perm_c = np.array([[9, 8, "a", 11], [8, "a", 9, 11], ["a", 9, 8, 11],])
    perm_d = np.array([[21, "x", "a", 0], ["x", "a", 21, 0], ["a", 21, "x", 0],])

    # Technically the same permutation though at different site.
    # (map needs to be more flexible)
    perm_e = np.array([[11, 8, 9, "a"], [11, 9, "a", 8], [11, "a", 8, 9],])
    # Symbol change + row change
    perm_f = np.array([[8, "a", 9, 11], [9, 8, "a", 11], ["a", 9, 8, 11],])

    pgs = [PatternMaker(x) for x in (perm_a, perm_b, perm_c, perm_d, perm_e, perm_f)]
    assert all(pg.label == pgs[0].label for pg in pgs)
    assert all(
        pg.label == PatternMaker.get_label(x)
        for pg, x in zip(pgs, (perm_a, perm_b, perm_c, perm_d, perm_e, perm_f))
    )


@pytest.fixture
def pm():
    perms = np.array([[8, 9, "a", 11], [9, "a", 8, 11], ["a", 8, 9, 11],])
    return PatternMaker(perms)


def test_perm_rep(pm):
    """In canon form."""
    p = pm._perms
    assert (p == np.array([[1, 2, 0, 3], [2, 0, 1, 3], [0, 1, 2, 3],])).all()
    assert p.dtype == "int64"


def test_bit_perm_rep(pm):
    """In bit form."""
    bp = pm._bit_perm
    assert (bp == np.array([[2, 4, 1, 8], [4, 1, 2, 8], [1, 2, 4, 8],])).all()
    assert bp.dtype == "int64"


def test_large_bit_perm_rep():
    """When large, object type should change to object"""
    # NOTE: not a proper group (incomplete), but this test
    # should not be affected by that.
    large_perm = np.array([list(range(64)), list(range(64))[::-1]])
    bp = PatternMaker(large_perm)._bit_perm
    assert bp.dtype == "object"


def test_search(pm):
    correct_answer = {
        0: [set()],
        1: [{"11"}, {"8"}],
        2: [{"8", "11"}, {"8", "9"}],
        3: [{"11", "9", "a"}, {"8", "9", "a"}],
        4: [{"8", "9", "a", "11"}],
    }
    index_map = pm.get_index_map()
    for n, answers in correct_answer.items():
        for _, p in pm.ap(n):
            assert set(index_map[x] for x in p) in answers


# Substitutor functions


@pytest.mark.parametrize(
    "from_species, to_species",
    [
        (("Fe1",), ("Fe7Ti",)),
        (("Fe1",), ("FeTi",)),
        (("Fe1", "Fe2"), ("Fe7Ti", "Fe3Ti")),
        (("Fe",), ("FeTiSnAu",)),
    ],
)
@chdir("../examples")
def test_all(from_species, to_species):
    """Integrated test with multi-color multi-orbit structure.."""
    if to_species == ("FeTiSnAu",):
        with pytest.raises(TooBigError):
            sh = ScriptHelper(
                structure_file="SmFe12.cif",
                from_species=from_species,
                to_species=to_species,
            )
            assert len(list(sh.substitutor.make_patterns())) == sh.count()
    else:
        sh = ScriptHelper(
            structure_file="SmFe12.cif",
            from_species=from_species,
            to_species=to_species,
        )
        assert len(list(sh.substitutor.make_patterns())) == sh.count()


@chdir("../examples")
def test_need_supercell():
    """
    Test whether program correctly exits if the structure
    needs supercell
    """
    with pytest.raises(NeedSupercellError):
        ScriptHelper(
            structure_file="SmFe12.cif", from_species=("Fe1",), to_species=("Fe12Ti1",),
        )


@chdir("../examples")
def test_non_cif():
    """
    Not-cif should raise ValueError
    """
    with pytest.raises(ValueError) as excinfo:
        ScriptHelper(
            structure_file="example1.py",
            from_species=("Fe1",),
            to_species=("Fe12Ti1",),
        )
        assert "only accept CIF" in str(excinfo.value)


@chdir("../examples")
def test_sequential():
    """
    Test sequential use of Substitutor;
    basically testing the setter of Substitutor.structure
    """
    structure = LabeledStructure.from_file("SmFe12.cif")
    structure1 = structure.copy()
    structure2 = structure.copy()
    structure1.replace_species({"Fe1": "Fe7Ti"})
    structure2.replace_species({"Fe2": "Fe6Ti2"})
    structure1 *= [1, 2, 1]
    structure2 *= [1, 2, 1]

    substitutor = Substitutor(structure1)
    assert substitutor.count() == 11
    assert len(list(substitutor.weights())) == 11
    substitutor.structure = structure2
    assert substitutor.count() == 147
    assert len(list(substitutor.weights())) == 147

@chdir("../examples")
def test_sequential_scaling_diagonal_one_scalar_value():
    """
    Test sequential use of Substitutor;
    basically testing the setter of Substitutor.structure
    """
    structure = LabeledStructure.from_file("SmFe12.cif")
    structure1 = structure.copy()
    structure2 = structure.copy()
    structure1.replace_species({"Fe1": "Fe7Ti"})
    structure2.replace_species({"Fe2": "Fe6Ti2"})
    structure1 *= [2]
    structure2 *= [2]

    substitutor = Substitutor(structure1)
    assert substitutor.count() == 17324048
    substitutor.structure = structure2
    assert substitutor.count() == 1909076572380

@chdir("../examples")
def test_sequential_scaling_diagonal_three_scalar_values():
    """
    Test sequential use of Substitutor;
    basically testing the setter of Substitutor.structure
    """
    structure = LabeledStructure.from_file("SmFe12.cif")
    structure1 = structure.copy()
    structure2 = structure.copy()
    structure1.replace_species({"Fe1": "Fe7Ti"})
    structure2.replace_species({"Fe2": "Fe6Ti2"})
    structure1 *= [1, 2, 1]
    structure2 *= [1, 2, 1]

    substitutor = Substitutor(structure1)
    assert substitutor.count() == 11
    assert len(list(substitutor.weights())) == 11
    substitutor.structure = structure2
    assert substitutor.count() == 147
    assert len(list(substitutor.weights())) == 147

@chdir("../examples")
def test_sequential_scaling_diagonal_matrix():
    """
    Test sequential use of Substitutor;
    basically testing the setter of Substitutor.structure
    """
    structure = LabeledStructure.from_file("SmFe12.cif")
    structure1 = structure.copy()
    structure2 = structure.copy()
    structure1.replace_species({"Fe1": "Fe7Ti"})
    structure2.replace_species({"Fe2": "Fe6Ti2"})
    structure1 *= [[1, 0, 0],[0, 2, 0],[0, 0, 1]]
    structure2 *= [[1, 0, 0],[0, 2, 0],[0, 0, 1]]

    substitutor = Substitutor(structure1)
    assert substitutor.count() == 11
    assert len(list(substitutor.weights())) == 11
    substitutor.structure = structure2
    assert substitutor.count() == 147
    assert len(list(substitutor.weights())) == 147

@chdir("../examples")
def test_sequential_scaling_nondiagonal_matrix():
    """
    Test sequential use of Substitutor;
    basically testing the setter of Substitutor.structure
    """
    structure = LabeledStructure.from_file("SmFe12.cif")
    structure1 = structure.copy()
    structure2 = structure.copy()
    structure1.replace_species({"Fe1": "Fe7Ti"})
    structure2.replace_species({"Fe2": "Fe6Ti2"})
    structure1 *= [[0, 1, 0],[2, 0, 0],[0, 0, 1]]
    structure2 *= [[0, 1, 0],[2, 0, 0],[0, 0, 1]]

    substitutor = Substitutor(structure1)
    assert substitutor.count() == 11
    assert len(list(substitutor.weights())) == 11
    substitutor.structure = structure2
    assert substitutor.count() == 147
    assert len(list(substitutor.weights())) == 147

@chdir("../examples")
def test_no_disorder():
    """No disorder sites should results in original structure."""
    structure = LabeledStructure.from_file("SmFe12.cif")
    substitutor = Substitutor(structure)
    assert substitutor.count() == 1
    assert list(substitutor.letters()) == [""]
    assert list(substitutor.weights()) == [1]
    assert len(list(substitutor.cifwriters())) == 1

@chdir("../examples")
def test_cifwriter():
    #Test cifwriter implementation.
    sh = ScriptHelper("SmFe7Ti.cif")
    sh.write()
    cifs = glob.glob("shry-SmFe*/slice*/*.cif")
    ref_cifs = glob.glob("../tests/test_cifs/smfe7ti/slice*/*.cif")

    def give_arbitrary_charge(filename):
        structure = LabeledStructure.from_file(filename)
        structure.add_oxidation_state_by_element({"Sm": 1, "Fe": 2, "Ti": 3})
        return structure

    try:
        esums = [EwaldSummation(give_arbitrary_charge(x)).total_energy for x in cifs]
        esums_ref = [
            EwaldSummation(give_arbitrary_charge(x)).total_energy for x in ref_cifs
        ]
        assert len(set(esums)) == 16
        assert set(esums) == set(esums_ref)
    finally:
        # Cleanup
        shry_outdirs = glob.glob("shry-SmFe*")
        for outdir in shry_outdirs:
            shutil.rmtree(outdir)

    sh = ScriptHelper("SmFeTi.cif", write_symm=True)
    sh.write()
    structures = [Structure.from_file(x) for x in glob.glob("shry-SmFe*/slice*/*.cif")]
    ref_structures = [Structure.from_file(x) for x in glob.glob("../tests/test_cifs/smfe7ti_sym/slice*/*.cif")]

    try:
        assert any(any(x==structure for x in ref_structures) for structure in structures)
    finally:
        # Cleanup
        shry_outdirs = glob.glob("shry-SmFe*")
        for outdir in shry_outdirs:
            shutil.rmtree(outdir)

@chdir("../examples")
def test_cifwriter2():
    """Test cifwriter edge cases."""
    structure = LabeledStructure.from_file("r3m.cif")
    s = Substitutor(structure)
    assert len(list(s.make_patterns())) == s.count()


@chdir("../examples")
def test_structure():
    """Test Structure generation."""
    parent_structure = LabeledStructure.from_file("r3m.cif")
    s = Substitutor(parent_structure)
    structures = []
    for structure in s.structure_writers():
        structures.append(structure)
    # TODO: Should properly check the content
    assert len(structures) == s.count()


@chdir("../examples")
def test_ewald():
    """Test ewald energy calculation."""

    def give_arbitrary_charge(filename):
        structure = LabeledStructure.from_file(filename)
        structure.add_oxidation_state_by_element({"Sm": 1, "Fe": 2, "Ti": 3})
        return structure

    structure = give_arbitrary_charge("SmFe7Ti.cif")
    s = Substitutor(structure)
    esums = list(s.ewalds())
    assert len(set(esums)) == 16

    structure = LabeledStructure.from_file("SmFe7Ti.cif")
    s = Substitutor(structure)
    with pytest.raises(ValueError) as excinfo:
        list(s.ewalds())
        assert "defined oxidation" in str(excinfo.value)


@pytest.mark.skip(reason="Feature dropped.")
@chdir("../examples")
def test_matheval():
    """
    Test ScriptHelper._math_eval() for various ScriptHelper.sample specification
    """
    sh = ScriptHelper("SmFe12.cif", sample="2/3*10000")
    assert sh.sample == 6666


@pytest.fixture
def polya():
    """Fixture returning a Polya instance."""
    perm_a = np.array([[8, 9, 10, 11], [9, 10, 8, 11], [10, 8, 9, 11],])
    perm_b = np.array([[3, 4, 5], [5, 4, 3], [3, 4, 5],])
    perms_list = [perm_a, perm_b]
    return Polya(perms_list)


def test_ci(polya):
    """Test cycle index calculation."""
    a = IndexedBase("a")
    b = IndexedBase("b")
    assert polya.sym_ci() == {
        0: a[1] ** 4 * b[1] ** 3,
        1: a[1] * a[3] * b[1] * b[2],
        2: a[1] * a[3] * b[1] ** 3,
    }


def test_count(polya):
    """Test counting of pattern. One should be enough representative."""
    assert polya.count(((3, 1), (2, 1))) == 5


#@pytest.mark.skip(reason="Comprehensive but time consuming. It will be activated later.")
@chdir("../benchmarks/03scailing_benchmark")
def test_benchmark():
    """benchmark / the number of symmetry-inequivalent structures."""
    df = pd.read_excel("./benchmark_SG_all.xls")
    for row, zipped in enumerate(
        zip(
            df["Supercell"],
            df["File"],
            df["Substitutions"],
            df["Equivalent Structures"],
            df["Checked"],
        )
    ):
        supercell, filename, substitution, equivalent, checked = zipped

        cif_basename = os.path.basename(filename).replace(".cif", "")
        filename = filename.replace(".cif", "_partial.cif")
        print(f"filename={filename}")
        structure = LabeledStructure.from_file(filename)
        supercell_size = list(map(int,supercell.split("x")))
        print(supercell_size)
        structure *= supercell_size
        s = Substitutor(structure,symprec=SHRY_TOLERANCE,angle_tolerance=SHRY_ANGLE_TOLERANCE)
        count_obtained=s.count()
        count_ref=equivalent
        print(f"count_obtained={count_obtained}")
        print(f"count_ref={count_ref}")
        assert count_obtained == count_ref
