# -*- coding: utf-8 -*-
# pylint: disable=wrong-import-order

"""
Test command line script functions
"""

import subprocess

import pytest
from pymatgen.core.composition import Composition
from shry import const
from shry.main import ScriptHelper

from helper import chdir, pre_cleanup

STRUCTURE_FILE = "01aP_ICSD_CollCode20537.cif"


@chdir("structures")
def test_scripthelper_matheval():
    """
    Test ScriptHelper._math_eval() for various ScriptHelper.sample specification
    """
    script_helper = ScriptHelper(STRUCTURE_FILE, sample="2/3*10000")
    assert script_helper.sample == 6666


@chdir("../examples")
def test_scripthelper_fromfile():
    """
    Test ScriptHelper.from_file() defaults and const.FLEXIBLE_SEPARATOR
    """
    script_helper = ScriptHelper.from_file("example.ini")
    assert script_helper.from_species == ["Pb", "Sn"]
    assert script_helper.to_species == ["Pb2+0.3Pb3+0.5Sn3+0.1", "Sn2+"]
    assert script_helper.scaling_matrix == [2, 2, 2]
    assert script_helper.symmetrize == const.DEFAULT_SYMMETRIZE
    assert script_helper.sample == const.DEFAULT_SAMPLE
    assert script_helper.symprec == const.DEFAULT_SYMPREC
    assert script_helper.dir_size == const.DEFAULT_DIR_SIZE
    assert script_helper.write_symm == const.DEFAULT_WRITE_SYMM

@chdir("../examples")
def test_script_ini():
    """
    Test *.ini loading from script (rough)
    """
    subprocess.run(["shry", "example.ini"], check=True)


@chdir("../examples")
def test_inputfile():
    """Test input file reading."""
    task = ScriptHelper.from_file("complete.ini")

    # config
    # mandatories
    assert task.from_species == ["Fe"]
    assert task.to_species == [Composition("Fe54TiCo")]

    # optionals
    assert task.structure_file == "NdFeB.POSCAR"
    assert task.sample == 10
    assert task.symprec == 1.0e-1


@chdir("../examples")
def test_inputfile_defaults():
    """Test defaults when providing minimum input file."""
    task = ScriptHelper.from_file("defaults.ini")

    # config
    # mandatories
    assert task.from_species == ["Fe"]
    assert task.to_species == [Composition("Fe54TiCo1")]

    # optionals
    assert task.structure_file == "POSCAR"
    assert task.sample == -1
    assert task.symprec == 1.0e-2


@chdir("../examples")
def test_target_site():
    """New target_site format for allowing symmetry-selective substitutions."""
    assert ScriptHelper.from_file("new_target_site.ini")


@chdir("../examples")
def test_no_file():
    """Non-existent input file."""
    with pytest.raises(FileNotFoundError):
        ScriptHelper.from_file("xyz.ini")


@chdir("../examples")
def test_invalid_tasks():
    """Don't accept invalid tasks"""
    with pytest.raises(ValueError, match="task: Invalid task"):
        ScriptHelper.from_file("err_invalid_tasks.ini")


@chdir("../examples")
def test_invalid_sections():
    """Don't accept invalid sections"""
    with pytest.raises(ValueError, match=r"input file: Invalid section\(s\) found"):
        ScriptHelper.from_file("err_invalid_sections.ini")


@pre_cleanup
@chdir("cmdline_parser")
def test_all():
    """Test all arguments."""
    subprocess.run(
        [
            "shry",
            "-t",
            "config",
            "-s",
            "Sb",
            "-c",
            "Sb0.875Bi0.125",
            "-a",
            "1.0e-3",
            "-p",
            "1.0e-3",
            "-f",
            "2_2_2_ICSD_CollCode430061.cif",
            "-r",
            "canaug",
            "-n",
            "0",
            "-m",
            "2",
            "1",
            "1",
        ],
        check=True,
    )


pbsn_site = [False, False, False, False, True, True, True, True]
te_site = [True, True, True, True, False, False, False, False]


@pytest.mark.parametrize(
    "structure_file, substitute, to, target_sites, scaling_matrix",
    [
        ("PbSnTe2.cif", "Pb", "Pb3+", [pbsn_site], "2 2 1"),  # simple, single
        ("PbSnTe2.cif", "Sn2+", "Pb2+Pb3+", [pbsn_site], "2 2 1"),  # same site
        (
            "PbSnTe2.cif",
            "Te1",
            "Pb0.5",
            [te_site],
            "2 2 1",
        ),  # some vacancy -> Should be X
        ("PbSnTe2.cif", "Te1", "X", [te_site], "2 2 1"),  # vacancy -> Should be X
        (
            "PbSnTe2.cif",
            "Sn1,Te2-",
            ["Pb2+0.625Pb3-0.125Sn0.25", "Se"],
            [pbsn_site, te_site],
            "2 1 1",
        ),  # whoa
    ],
)
@chdir("canaug")
@pre_cleanup
def test_substitution(structure_file, substitute, to, target_sites, scaling_matrix):
    task = ScriptHelper(
        structure_file=structure_file,
        from_species=substitute,
        to_species=to,
        scaling_matrix=scaling_matrix,
        extension="POSCAR",
    )
    assert task.target_sites == target_sites
    # print(task.adjusted_structure.site_properties)
    task.write()


@chdir("canaug")
def test_substitution_overlap():
    with pytest.raises(ValueError):
        task = ScriptHelper(
            structure_file="PbSnTe2.cif",
            from_species="Pb2+,Sn",
            to_species=["Pb", "Ge"],
        )


@chdir("canaug")
def test_substitution_fail(exception):
    # Like above but for failing case
    ...


def test_part_sites(run):
    ...
    # Multiple Te orbit, but choose only those with specific
    # target_sites = "Te2+"
    # target_compositions = "Pb2+0.625Pb3+0.125Sn0.125"
