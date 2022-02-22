"""Random select substitution; save substituted structure and JSON info"""
import warnings
warnings.simplefilter('ignore')

import errno
import functools
import glob
import math
import os
import random
import re
import signal
import sys

import numpy as np
import pandas as pd
import pymatgen
import shry
from ase import Atoms
from ase.io import read, write
from ase.spacegroup import Spacegroup
from pymatgen.core.composition import Composition
from pymatgen.core.periodic_table import get_el_sp
from pymatgen.io.cif import CifParser
from pymatgen.util.string import formula_double_format

shry.const.DISABLE_PROGRESSBAR = True
from shry.core import (NeedSupercellError, PatchedSpacegroupAnalyzer,
                       Substitutor, TooBigError)
from shry.main import LabeledStructure

CONFIG_IRREDUCIBLE_MAX = 3e3
CONFIG_IRREDUCIBLE_MIN = 1e1
TOTAL_SUBSTITUTION_MIN = 1e1
TOTAL_SUBSTITUTION_MAX = 1e4
CONFIG_TIMEOUT = 120

CONFIG_CELL_MAX = 500
CONFIG_CELLVEC_MIN = 1
CONFIG_CELLVEC_MAX = 3
CONFIG_RETRY_N = 500

MAX_SUBBED_MIN=3      # default 3
MAX_NSPIECIES_MIN=4   # default 4

SHRY_TOLERANCE=0.01      #angstrom
SHRY_ANGLE_TOLERANCE=5.0 #degree

# "Manual" periodic table for searching element within the same group
# (just so that it looks like make sense)
PT = [
    "H,Li,Na,K,Rb,Cs,Fr".split(","),  # Lazy writing
    "D,Li,Na,K,Rb,Cs,Fr".split(","),  # Workaround for Deuterium
    ",Be,Mg,Ca,Sr,Ba,Ra".split(","),
    ",,,Sc,Y,Lu,Lr".split(","),
    ",,,Ti,Zr,Hf,Rf".split(","),
    ",,,V,Nb,Ta,Db".split(","),
    ",,,Cr,Mo,W,Sg".split(","),
    ",,,Mn,Tc,Re,Bh".split(","),
    ",,,Fe,Ru,Os,Hs".split(","),
    ",,,Co,Rh,Ir,Mt".split(","),
    ",,,Ni,Pd,Pt,Ds".split(","),
    ",,,Cu,Ag,Au,Rg".split(","),
    ",,,Zn,Cd,Hg,Cn".split(","),
    ",B,Al,Ga,In,Tl,Nh".split(","),
    ",C,Si,Ge,Sn,Pb,Fl".split(","),
    ",N,P,As,Sb,Bi,Mc".split(","),
    ",O,S,Se,Te,Po,Lv".split(","),
    ",F,Cl,Br,I,At,Ts".split(","),
    "He,Ne,Ar,Kr,Xe,Rn,Og".split(","),
    # Here I just grouped like this; but anyway combinatorially identical
    "La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb".split(
        ","
    ),  # Lu omitted because it's there
    "Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No".split(","),  # Lr is up there
]
PTDF = pd.DataFrame(PT, dtype=str).T
PTDF.fillna("", inplace=True)
PTDF = PTDF.convert_dtypes()

# Copied from SHRY
# Operates on chemical formula
COMPONENT = re.compile(r"[A-z][a-z]*[0-9.\+\-]*[0-9.]*")
# Operates on single component
AMOUNT = re.compile(r"(?<=[\+\-A-Za-z])[0-9.]+(?![\+\-])")
# WARNING: Does not process string with spaces properly
SPECIES = re.compile(r"[A-Z][a-z]*")
OXSTATE = re.compile(r"[0-9][-+]*")


class TimeoutError(Exception):
    pass


def timeout(seconds=CONFIG_TIMEOUT, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.setitimer(signal.ITIMER_REAL, seconds)  # used timer instead of alarm
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator


# def random_select_from_group(number, howmuch, butnot):
def random_select_from_group(species, howmuch, exclude):
    number = np.where(PTDF == species)[1][0]

    group = PTDF.iloc[:, number]
    # notnone = list(group[(group != "None") & (group != butnot)])
    notnone = group != ""
    notbutnot = functools.reduce(
        lambda x, y: x & y,
        [group != x for x in exclude],
        [True] * len(notnone),
    )
    # notnone = list(group[(group != "None") & (group not in butnot)])
    okgroup = list(group[notnone & notbutnot])
    # return random.sample(notnone, howmuch)
    if howmuch > len(okgroup):
        # Contingency: select anything from anywhere
        group = PTDF.to_numpy().flatten()
        notnone = group != ""
        notbutnot = functools.reduce(
            lambda x, y: x & y,
            [group != x for x in exclude],
            np.ones(notnone.shape, dtype=bool),
        )
        okgroup = list(group[notnone & notbutnot])
    return random.sample(okgroup, howmuch)


def get_oxstate(string):
    ox = OXSTATE.findall(string)
    if ox:
        ox = ox[0]
        if "-" in ox:
            return -1 * int(ox.strip("+-"))
        return int(ox.strip("+-"))
    else:
        return 0


def to_oxstate_string(number):
    if number >= 0:
        return f"{number}+"
    return f"{abs(number)}-"

def remove_label(string):
    return re.sub(OXSTATE, '', string)

def random_scale_and_substitute(cif_filename):
    """
    Random select supercell and substitution
    """

    @timeout()
    def get_count(ct_structure):
        substitutor = Substitutor(
            ct_structure,
            symprec=SHRY_TOLERANCE,
            angle_tolerance=SHRY_ANGLE_TOLERANCE
            )
        return substitutor.count()

    def beautiful_inted_formula(composition):
        inted_element_composition = composition.inted_composition.element_composition
        # Because buggy formula!
        sym_amt = inted_element_composition.get_el_amt_dict()
        syms = sorted(sym_amt.keys(), key=lambda sym: get_el_sp(sym).X)
        formula = [s + formula_double_format(int(sym_amt[s]), False) for s in syms]
        return "".join(formula)
        # return composition.inted_composition.element_composition.formula.replace(
        #     " ", ""
        # )

    parser = CifParser(cif_filename)

    key = list(parser.as_dict().keys())[0]
    cif_dict = parser.as_dict()[key]
    chem_formula = cif_dict["_chemical_formula_sum"].replace(" ", "")

    structure = LabeledStructure.from_file(cif_filename)
    #print(structure)
    sga = PatchedSpacegroupAnalyzer(structure,
                symprec=SHRY_TOLERANCE,
                angle_tolerance=SHRY_ANGLE_TOLERANCE
            )
    if "_space_group_IT_number" in cif_dict:
        cif_sg_number = cif_dict["_space_group_IT_number"]
    else:
        cif_sg_number = cif_dict["_symmetry_Int_Tables_number"]
    assert sga.get_space_group_number() == int(cif_sg_number)
    
    lattice_type = sga.get_lattice_type()
    point_group = sga.get_point_group_symbol()
    space_group_num = sga.get_space_group_number()
    space_group = sga.get_space_group_symbol()
    #structure_id = "icsd-" + cif_dict["_database_code_ICSD"]
    structure_id = "cod-" + cif_dict["_cod_database_code"]
    
    for attempt in range(CONFIG_RETRY_N):
        # Select random supercell, that fits below max
        # Brute force; not taking time anyway
        print(f"---attempt={attempt}---")
        natom = len(structure)
        if natom <= CONFIG_CELL_MAX:
            while True:
                scaling_matrix = [random.randint(CONFIG_CELLVEC_MIN, CONFIG_CELLVEC_MAX) for i in range(3)]
                enlargement = functools.reduce(lambda x, y: x * y, scaling_matrix)
                if (enlargement * natom) <= CONFIG_CELL_MAX:
                    break
        else:
            scaling_matrix = [1, 1, 1]
            enlargement = 1

        print(f"scaling_matrix={scaling_matrix}")

        #print([x.is_ordered for x in structure])
        fully_ordered = all(x.is_ordered for x in structure)
        if not fully_ordered:
            print("Partial occupancies in some of the sites.")
            raise ValueError
            continue

            """
            chk_sup_structure = structure.copy()
            chk_sup_structure *= scaling_matrix
            try:
                substitutor = Substitutor(chk_sup_structure)
                count = substitutor.count()
            except (TimeoutError, TooBigError):
                continue
            print(f"inside fully_ordered.")
            print(f"count={count}")

            if count > CONFIG_IRREDUCIBLE_MAX:
                print(
                    f"Scaling matrix {scaling_matrix} too big! "
                    f"({count} structures without substitution)"
                )
                continue
            
            if count < CONFIG_IRREDUCIBLE_MIN:
                print(
                    f"Scaling matrix {scaling_matrix} is not appropriate (probably too small)! "
                    f"({count} structures without substitution)"
                )
                continue
            """

        # Use SymmetrizedStructure to get some properties
        sym_structure = PatchedSpacegroupAnalyzer(structure, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE).get_symmetrized_structure()
        # Multiplicities of each equivalent sites
        multiplicities = [len(x) for x in sym_structure.equivalent_indices]
        supercell_mul = [enlargement * len(x) for x in sym_structure.equivalent_indices]
        # Index of which sites has _more than one_ multiplicities
        mulsites = [i for i, m in enumerate(supercell_mul) if m > 1]

        # Select random amount of positions to be substituted as index
        # From the index of m > 1 sites, randomly select 1-all sites to be substituted
        # Cap the selected sites to 3 as more tend to break the limit!
        # If no multiplicities: try other supercell
        if not len(mulsites):
            continue
        max_subbed = min((len(mulsites), MAX_SUBBED_MIN))
        #print(f"max_subbed={max_subbed}")
        subbed = random.sample(mulsites, random.randint(1, max_subbed))

        partitions = []
        for i in subbed:
            # Substitute to random N amount of final species
            max_nspecies = min((supercell_mul[i], MAX_NSPIECIES_MIN))
            nspecies = random.randint(2, max_nspecies)  # limit to 4
            #print(f"i={i}, nspecies={nspecies}")

            # Random select each amount of final species
            ranges = []
            left = 0  # Left part of the range
            pad = nspecies - 2  # Padding to ensure valid range is always selected
            for _ in range(nspecies - 1):  # select but the last
                margin = supercell_mul[i] - pad
                right = random.randrange(left + 1, margin)  # at least choose one
                ranges.append(range(left, right))
                pad -= 1
                left = right
            right = supercell_mul[i]  # move right
            ranges.append(range(left, right))
            partitions.append([len(x) for x in ranges])

        # Convert the final ratios to fractions (or pymatgen will protest)
        fraced = [[x / sum(p) for x in p] for p in partitions]

        # Actually substitute
        exclude_species = set()
        exclude_species |= set(SPECIES.findall(chem_formula.replace(" ", "")))
        for e, s in enumerate(subbed):
            lead_i = sym_structure.equivalent_indices[s][0]
            composition = sym_structure[lead_i].species
            composition_str = str(composition).replace(" ", "")

            # Better
            components = COMPONENT.findall(composition_str)
            specieses = [SPECIES.findall(x)[0] for x in components]
            oxstates = [get_oxstate(x) for x in components]
            oxstate_strings = [to_oxstate_string(x) for x in oxstates]

            if len(components) > 1:
                for species in specieses:
                    exclude_species -= {species}

                # Just choose first one
                species = specieses[0]
                oxstate_string = oxstate_strings[0]
            else:
                species = specieses[0]
                oxstate_string = oxstate_strings[0]

                # Don't exclude the initial species from current position
                exclude_species -= {species}

            sub_species = random_select_from_group(
                species, len(fraced[e]), exclude_species
            )

            # Perhaps (hopefully) doesn't matter
            # exclude_species |= set(sub_species)

            target_composition_dict = {
                get_el_sp(s + oxstate_string): f for s, f in zip(sub_species, fraced[e])
            }
            # Use SHRY's patch to Composition
            target_composition = Composition(target_composition_dict)

            for i in sym_structure.equivalent_indices[s]:
                sym_structure.replace(
                    i, target_composition, properties=sym_structure[i].properties
                )

        combinations = []
        for p in partitions:
            n_part = sum(p)
            p_max = max(p)
            p_max_i = p.index(p_max)
            pc = p.copy()
            pc.pop(p_max_i)
            c = 1
            n_r = n_part
            for e in pc:
                c *= math.comb(n_r, e)
                n_r -= e
            combinations.append(c)
        substitutions = functools.reduce(lambda x, y: x * y, combinations)

        if substitutions >= TOTAL_SUBSTITUTION_MAX:
            print(f"Total substitutions = {substitutions:.3e} is too large >= {TOTAL_SUBSTITUTION_MAX:.3e}: rejected!", flush=True)
            print("============================================")
            # Instead return the substitution configuration
            continue
        elif substitutions <= TOTAL_SUBSTITUTION_MIN:
            print(f"Total substitutions = {substitutions:.3e} is too small <= {TOTAL_SUBSTITUTION_MIN:.3e}: rejected!", flush=True)
            print("============================================")
            # Instead return the substitution configuration
            continue
        else:
            print("TOTAL_SUBSTITUTION:")
            print(f"{TOTAL_SUBSTITUTION_MIN:.3e} < {substitutions:.3e} < {TOTAL_SUBSTITUTION_MAX:.3e}", flush=True)

        equivalent_labels = [
            list(x[0].properties["_atom_site_label"])[0]
            for x in sym_structure.equivalent_sites
        ]

        equivalent_formulas = [
            set(beautiful_inted_formula(sym_structure[i].species) for i in g)
            for g in sym_structure.equivalent_indices
        ]
        assert all(len(x) == 1 for x in equivalent_formulas)
        equivalent_formulas = [list(x)[0] for x in equivalent_formulas]

        wyckoffs = [
            f"{remove_label(equivalent_labels[x])} ({sym_structure.wyckoff_symbols[x]}) $\\rightarrow$ {equivalent_formulas[x]}"
            for x in subbed
        ]
        sub_chem_formula = beautiful_inted_formula(sym_structure.composition)

        config = {
            "Compound": chem_formula,
            "ID": structure_id,
            "LatticeType": lattice_type,
            "PointGroup": point_group,
            "SpaceGroup_No": space_group_num,
            "SpaceGroup": space_group,
            "Compositions": sub_chem_formula,
            "Wyckoffs": ", ".join(wyckoffs),
            "Supercell": "x".join(map(str, scaling_matrix)),
            "Substitutions": substitutions,
            "Note": "autosub:success",
            "Checked": None
        }

        try:
            large_structure = sym_structure.copy()
            large_structure *= scaling_matrix
            print(f"get_count starts...", flush=True)
            count = get_count(large_structure)
            #print(f"count = {count}")
            config["Equivalent Structures"] = count
        except (TimeoutError, TooBigError):
            # If too long, likely too big!
            print("TIMEOUT/ sub. was too long => likely too big")
            #print("============================================")
            continue
        except MemoryError:
            print("MemoryError/ memory overflow => likely too big")
            #print("============================================")
            continue
        except NeedSupercellError as e:
            print(config)
            print(enlargement)
            print(multiplicities)
            raise e

        if count <= CONFIG_IRREDUCIBLE_MAX and count >= CONFIG_IRREDUCIBLE_MIN:
            print(f"Expected {count:.3e} structures with this substitution: approved!")
            print("============================================")
            # Instead return the substitution configuration
            return sym_structure, config
        print(
            f"Expected {count:.3e} structures with this substitution: try again! (attempt {attempt})"
        )
        #print("============================================")

    # If all failed then return the default config: No supercell no substitution plain single structure
    config = {
        "Compound": chem_formula,
        "ID": structure_id,
        "LatticeType": lattice_type,
        "PointGroup": point_group,
        "SpaceGroup_No": space_group_num,
        "SpaceGroup": space_group,
        "Compositions": chem_formula,
        "Wyckoffs": "",
        "Supercell": "1x1x1",
        "Substitutions": 1,
        "Equivalent Structures": 1,
        "Note": "autosub:failure",
        "Checked": None
    }
    return sym_structure, config

def remove_glob(pathname, recursive=True):
    for p in glob.glob(pathname, recursive=recursive):
        if os.path.isfile(p):
            os.remove(p)
            
def main():
    paths = sys.argv[1:]
    #paths = [f"SG{i}" for i in range(1,231)]
    #paths = [f"SG{i}" for i in range(10,11)]
    series = []
    try:
        for path in paths:
            #print(os.path.join(os.getcwd(),path))
            if os.path.isdir(path):
                root_dir = path
                remove_glob(os.path.join(root_dir,'*_partial.cif'))
                cifs = glob.glob(os.path.join(root_dir, "*.cif"))
            else:
                continue
                #root_dir = os.path.dirname(path)
                #cifs = [path]

            sg = Spacegroup(int(root_dir.lstrip("SG").rstrip("/")))
            print("============================================")
            print(f"Space group: {sg.no} ({sg.symbol})")
            print("============================================")
            
            for cif in cifs:
                print(f"cif = {cif}")
                cif_dir = os.path.dirname(cif)
                cif_filename = os.path.basename(cif).rstrip(".cif")
                out_cif = os.path.join(cif_dir, cif_filename + "_partial.cif")

                #cif_instance = LabeledStructure.from_file(cif)
                #reduced_str = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(
                #    cif_instance
                #)
                #print("------------------------------------------")
                #print("This is the space group symmbol")
                #print(f"Space group= {reduced_str.get_space_group_symbol()}")
                #sym_str = reduced_str.get_symmetrized_structure()
                #print("This is the symmetrized structure")
                #print(sym_str)
                #print(f"=========================================")
                #print(f" ")
                #print(f" ")

                sub_str, subconfig = random_scale_and_substitute(cif)
                subconfig["File"] = cif
                serie = pd.Series(subconfig)
                series.append(serie)

                # Critical: do *NOT* refine_struct
                sub_str.to(filename=out_cif, symprec=0.01, refine_struct=False)
                out_parser = CifParser(out_cif)
                key = list(out_parser.as_dict().keys())[0]
                # Space group should _not_ change
                assert int(sg.no) == int(
                    out_parser.as_dict()[key]["_symmetry_Int_Tables_number"]
                )
            print("")

    finally:
        # Even if crash, the dataframe should be saved.
        df = pd.DataFrame(series)
        print(df)
        df.to_excel(f"autosub_SG{sg.no}.xls")


if __name__ == "__main__":
    main()
