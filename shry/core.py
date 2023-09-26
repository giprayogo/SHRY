# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Core operations, pattern generation, etc.
"""

# python modules
import collections
import functools
import itertools
import logging
import math
import random
import sys
from typing import OrderedDict, Tuple

# python modules
import numpy as np
import spglib
import sympy
import tqdm
from monty.fractions import gcd_float
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core.composition import Composition, reduce_formula
from pymatgen.io.cif import CifBlock, CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list_pbc
from pymatgen.util.string import transformation_to_string
from scipy.special import comb
from sympy.utilities.iterables import multiset_permutations
from tabulate import tabulate
from pymatgen.core.periodic_table import get_el_sp

# shry modules
from . import const

# shry version control
try:
    from ._version import version as shry_version
except (ModuleNotFoundError, ImportError):
    shry_version = "unknown"

np.seterr(all="raise")
np.set_printoptions(linewidth=const.LINEWIDTH, threshold=16)
np.set_printoptions(linewidth=1000, threshold=sys.maxsize)


def get_integer_formula_and_factor(
    self,
    max_denominator: int = int(1 / const.DEFAULT_ATOL),
    iupac_ordering: bool = False,
) -> Tuple[str, float]:
    """
    The default Composition groups together different ox states which is not ideal...
    """
    el_amt = self.as_dict()
    g = gcd_float(list(el_amt.values()), 1 / max_denominator)

    d = {k: round(v / g) for k, v in el_amt.items()}
    (formula, factor) = reduce_formula(d, iupac_ordering=iupac_ordering)
    if formula in Composition.special_formulas:
        formula = Composition.special_formulas[formula]
        factor /= 2
    return formula, factor * g


# Patched extra functionalities and bug fixes on top of Pymatgen's classes.
def to_int_dict(self):
    """
    Returns:
        Dict with element symbol and integer amount
    """
    _, factor = self.get_integer_formula_and_factor(
        max_denominator=int(1 / const.DEFAULT_ATOL)
    )
    int_dict = {e: int(a) for e, a in (self / factor).as_dict().items()}

    # be safe: Composition groups together different ox states which is not ideal...
    for x, y in zip(int_dict.values(), self.as_dict().values()):
        if not np.isclose(x * factor, y, atol=const.DEFAULT_ATOL):
            raise ValueError(
                "Composition (Occupancy) is not rational! Please try to increase significant digits "
                "e.g., 1/3 = 0.3333 -> 1/3 = 0.3333333333333."
            )

    return int_dict


@property
def inted_composition(self):
    """
    Return Composition instance with integer formula
    """
    _, factor = self.get_integer_formula_and_factor(
        max_denominator=int(1 / const.DEFAULT_ATOL)
    )
    int_comp = self / factor

    # be safe
    int_dict = {e: int(a) for e, a in int_comp.as_dict().items()}
    if not all(
        np.isclose(x * factor, y, atol=const.DEFAULT_ATOL)
        for x, y in zip(int_dict.values(), self.as_dict().values())
    ):
        raise ValueError(
            "Composition (Occupancy) is not rational! Please try to increase significant digits "
            "e.g., 1/3 = 0.3333 -> 1/3 = 0.3333333333333."
        )

    return int_comp


def formula_double_format_tol(
    afloat, ignore_ones=True, tol: float = const.DEFAULT_ATOL * 10
):
    """
    This function is used to make pretty formulas by formatting the amounts.
    Instead of Li1.0 Fe1.0 P1.0 O4.0, you get LiFePO4.

    Args:
        afloat (float): a float
        ignore_ones (bool): if true, floats of 1 are ignored.
        tol (float): Tolerance to round to nearest int. i.e. 2.0000000001 -> 2

    Returns:
        A string representation of the float for formulas.
    """
    if ignore_ones and afloat == 1:
        return ""
    if abs(afloat - round(afloat)) < tol:
        return round(afloat)
    return round(afloat, 8)


@property
def formula(self) -> str:
    """
    Returns a formula string, with elements sorted by electronegativity,
    e.g., Li4 Fe4 P4 O16.
    """
    sym_amt = self.get_el_amt_dict()
    syms = sorted(sym_amt, key=lambda sym: get_el_sp(sym).X)
    formula = [
        f"{s}{formula_double_format_tol(sym_amt[s], False)}" for s in syms
    ]
    return " ".join(formula)


Composition.to_int_dict = to_int_dict
Composition.get_integer_formula_and_factor = get_integer_formula_and_factor
Composition.inted_composition = inted_composition
Composition.formula = formula


class PatchedSymmetrizedStructure(SymmetrizedStructure):
    """
    Fixed site_properties display
    """

    def __str__(self):
        outs = [
            "SymmetrizedStructure",
            f"Full Formula ({self.composition.formula})",
            f"Reduced Formula: {self.composition.reduced_formula}",
        ]

        def to_s(x):
            return f"{x:0.6f}"

        outs.append(
            "abc   : "
            + " ".join([to_s(i).rjust(10) for i in self.lattice.abc])
        )
        outs.append(
            "angles: "
            + " ".join([to_s(i).rjust(10) for i in self.lattice.angles])
        )
        if self._charge:
            if self._charge >= 0:
                outs.append(f"Overall Charge: +{self.charge}")
            else:
                outs.append(f"Overall Charge: -{self._charge}")
        outs.append(f"Sites ({len(self)})")
        data = []
        props = self.site_properties  # This should be updated!
        keys = sorted(props.keys())
        for i, sites in enumerate(self.equivalent_sites):
            site = sites[0]
            row = [str(i), site.species_string]
            row.extend([to_s(j) for j in site.frac_coords])
            row.append(self.wyckoff_symbols[i])
            for k in keys:
                row.append(site.properties[k])  # This line
            data.append(row)
        outs.append(
            tabulate(
                data,
                headers=["#", "SP", "a", "b", "c", "Wyckoff"] + keys,
            )
        )
        return "\n".join(outs)


class PatchedSpacegroupAnalyzer(SpacegroupAnalyzer):
    """
    Patched to use PatchedSymmetrizedStructure in get_symmetrized_structure()
    """

    def get_symmetrized_structure(self):
        """
        Use PatchedSymmetrizedStructure
        """
        ds = self.get_symmetry_dataset()
        sg = SpacegroupOperations(
            self.get_space_group_symbol(),
            self.get_space_group_number(),
            self.get_symmetry_operations(),
        )
        return PatchedSymmetrizedStructure(
            self._structure, sg, ds["equivalent_atoms"], ds["wyckoffs"]
        )


class AltCifBlock(CifBlock):
    """
    CifBlock but optimized for when there are
    many, slightly different blocks:

    - Pre-format fields on init, value update, and from_string
    - Cache unchanged strings
    - Loop to list instead of to string
    """

    def _loop_to_list(self, loop):
        s = ["loop_"] + [f" {l}" for l in loop]

        # Join the fields
        lines = "\n".join(
            "  ".join(x) for x in zip(*[self.data[k] for k in loop])
        ).replace(";", "\n")

        # Check for maximum lengths
        for line in lines.split("\n"):
            if len(line) < self.maxlen:
                s.append(line)
            else:
                sublines = [
                    line[i : i + self.maxlen]
                    for i in range(0, len(line), self.maxlen)
                ]
                s.extend(sublines)
        return s

    def __str__(self):
        """
        Returns the cif string for the data block
        """
        s = [f"data_{self.header}"]
        keys = self.data.keys()
        written = []
        for k in keys:
            if k in written:
                continue
            # can this be optimized (won't affect much though)
            for l in self.loops:
                # search for a corresponding loop
                if k in l:
                    loop_id = "loop_\n " + "\n ".join(l)
                    if self.string_cache[loop_id] is not None:
                        s.extend(self.string_cache[loop_id])
                    else:
                        looplist = self._loop_to_list(l)
                        s.extend(looplist)
                        self.string_cache[loop_id] = looplist
                    written.extend(l)
                    break
            if k not in written:
                # k didn't belong to a loop
                if self.string_cache[k] is not None:
                    s.extend(self.string_cache[k])
                else:
                    v = self.data[k]
                    if len(k) + len(v) + 3 < self.maxlen:
                        s.append(f"{k}   {v}")
                    else:
                        s.extend([k, v])
        return "\n".join(s)

    def __setitem__(self, key, value):
        for l in self.loops:
            if key in l:
                self.data[key] = [self._format_field(x) for x in value]
                loop_id = "loop_\n " + "\n ".join(l)
                self.string_cache[loop_id] = None
                return
        self.data[key] = self._format_field(value)
        self.string_cache[key] = None

    def __init__(self, data, loops, header):
        super().__init__(data, loops, header)
        # Pre-format everything.
        formatted = []
        self.string_cache = OrderedDict()
        for k in self.data.keys():
            if k in formatted:
                continue
            for l in self.loops:
                if k in l:
                    for _k in l:
                        self.data[_k] = list(
                            map(self._format_field, self.data[_k])
                        )
                    loop_id = "loop_\n " + "\n ".join(l)
                    self.string_cache[loop_id] = self._loop_to_list(l)
                    formatted.extend(l)
                    break
            if k not in formatted:
                self.data[k] = self._format_field(self.data[k])
                # Cache preparation
                v = self.data[k]
                s = []
                if len(k) + len(v) + 3 < self.maxlen:
                    s.append(f"{k}   {v}")
                else:
                    s.extend([k, v])
                self.string_cache[k] = s


# I just want distinguishable exceptions


class NeedSupercellError(Exception):
    """
    Just so that this error will raise a distinguishable exception
    """

    ...


class TooBigError(Exception):
    """
    When the number of substitution is too large
    """

    ...


# Some numerical methods


def rec_asc(a, n, m, k, length):
    """
    Recursive part of the thing below
    """
    x = m
    while 2 * x <= n and k < length:
        a[k - 1] = x
        yield from rec_asc(a, n - x, x, k + 1, length)
        x += 1
    a[k - 1] = n
    # Pad
    b = a[:k]
    b += [0] * (length - k)
    # Permute
    for z in multiset_permutations(b):
        yield z


@functools.lru_cache(None)
def aR_array(n, length=None):
    """
    aR() cache.
    """
    a = [1] * n
    if length is None:
        length = n
    return np.array(list(rec_asc(a, n, 1, 1, length)))


@functools.lru_cache(None)
def multinomial_coeff(a):
    """
    Get multinomial coefficient of (sum(a), (a_1, a_2, ...))
    """
    return functools.reduce(
        lambda a, b: a * b,
        [comb(cx, x, exact=True) for cx, x in zip(np.cumsum(a), a)],
    )


class Substitutor:
    """
    Makes unique ordered Structure(s) for the given disorder Structure.

    Since the PatternMaker instances are saved,
    it is advantageous to reuse the same Substitutor
    for symmetrically similar structures.

    This includes different subsitute concentrations
    from the same basic structure, among others.

    Args:
        structure (pymatgen.Structure): Input structure.
        symprec (float): Symmetry precision.
        angle_tolerance (float): Angle tolerance for symmetry search.
        groupby (function): Function to group disordered sites.
            Defaults to lambda x: x.properties["_atom_site_label"], with the x loops
            over all PeriodicSites within the Structure.

            The proper function depends on what you want; do you want all sites of
            a certain species or crystallographic orbit to be subsituted together?

            Fallback to crystallographic orbit when failed.
        sample (int): Randomly sample the generated structures.
        no_dmat (bool): Whether or not to use distance matrix as
            permutation invariant (faster, default=True)
        t_kind (): (to be written)
        cache (bool or None): By default(=None), cache patterns when involving either
            multiple orbits or multiple species, but otherwise don't.
            Caching allows "reuse" of previously generated patterns,
            at the cost of memory.
            Set to False if memory is limited,
            but note that pattern generation will be much slower.
    """

    __slots__ = (
        "_symprec",
        "_angle_tolerance",
        "_groupby",
        "_shuffle",
        "_seed",
        "_atol",
        "_symmops",
        "_pms",
        "_enumerator_collection",
        "disorder_groups",
        "_group_dmat",
        "_group_perms",
        "_group_indices",
        "_group_bits",
        "_group_bit_perm",
        "_structure",
        "_sample",
        "cache",
        "_no_dmat",
        "_t_kind",
        "_segmenter",
        "_charset",
        "_template_cifwriter",
        "_template_structure",
        "_rg",
    )

    def __init__(
        self,
        structure,
        symprec=const.DEFAULT_SYMPREC,
        angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        atol=const.DEFAULT_ATOL,
        groupby=None,
        shuffle=False,
        seed=const.DEFAULT_SEED,
        no_dmat=const.DEFAULT_NO_DMAT,
        t_kind=const.DEFAULT_T_KIND,
        cache=None,  # "True", "False", "None" (default)
    ):
        self._symprec = symprec
        self._angle_tolerance = angle_tolerance
        self._no_dmat = no_dmat
        self._t_kind = t_kind
        # TODO: These two are consequential to self._structure, so should be property
        self._atol = atol
        self._shuffle = shuffle
        self._seed = seed
        if groupby is None:
            self._groupby = lambda x: x.properties["_atom_site_label"]
        else:
            self._groupby = groupby

        # random.Random instance.
        # The default RNG is used in many other modules,
        # so setting seed does not make sense there.
        rg = random.Random()
        rg.seed(self._seed)
        self._rg = rg

        self.cache = cache
        self._symmops = None
        self._pms = dict()
        self._enumerator_collection = PolyaCollection()

        self.disorder_groups = dict()
        self._group_dmat = dict()
        self._group_perms = dict()
        self._group_indices = dict()
        self._group_bits = dict()
        self._group_bit_perm = dict()
        # TODO: Decouple output-format specific attributes.
        self._template_cifwriter = None
        self._template_structure = None

        self.structure = structure

    @property
    def structure(self):
        """
        Input Structure containing disorder sites
        """
        return self._structure

    @property
    def groupby(self):
        """
        Function used to divide disorder sites
        """
        return self._groupby

    @groupby.setter
    def groupby(self, groupby):
        self._groupby = groupby
        # Reload division part
        self.structure = self.structure

    @structure.setter
    def structure(self, structure):
        logging.info("\nSetting Substitutor with Structure")
        logging.info(f"{structure}")

        # TODO: Avoid these kinds of manual invocations.
        self.disorder_groups.clear()
        self._group_dmat.clear()
        self._group_perms.clear()
        self._group_indices.clear()
        self._group_bits.clear()
        self._group_bit_perm.clear()
        self._template_cifwriter = None
        self._template_structure = None

        # Read.
        self._structure = structure.copy()

        sga = PatchedSpacegroupAnalyzer(
            self._structure,
            symprec=self._symprec,
            angle_tolerance=self._angle_tolerance,
        )
        try:
            self._symmops = sga.get_symmetry_operations()
        except TypeError as exc:
            raise RuntimeError("Couldn't find symmetry.") from exc

        logging.info(
            f"Space group: {sga.get_hall()} ({sga.get_space_group_number()})"
        )
        logging.info(f"Total {len(self._symmops)} symmetry operations")
        logging.info(sga.get_symmetrized_structure())
        equivalent_atoms = sga.get_symmetry_dataset()["equivalent_atoms"]

        # Identify and label disorder sites.
        # Index is used for later recoloring the sites.
        logging.info("\nProcessing disorder sites.")
        disorder_sites = []
        for i, site in enumerate(self._structure):
            site.properties["index"] = i
            site.properties["equivalent_atoms"] = equivalent_atoms[i]
            if not site.is_ordered:
                disorder_sites.append(site)
                # Ad hoc fix: if occupancy is less than 1, stop.
                # TODO: Automatic vacancy handling
                if not np.isclose(
                    site.species.num_atoms, 1.0, atol=self._atol
                ):
                    logging.warning(
                        f"The occupancy of the site {site.species} is {site.species.num_atoms}."
                    )
                    logging.warning(
                        f"This should be 1 within the torelance, atol={self._atol}."
                    )
                    logging.warning(
                        "If you want to consider vacancy sites, please add pseudo atoms."
                    )
                    raise RuntimeError(
                        "The sum of number of occupancies is not 1."
                    )
        if not disorder_sites:
            logging.warning("No disorder sites found within the Structure.")

        try:
            disorder_sites.sort(key=self._groupby)
        except KeyError:
            # Fallback to crystallographic orbit.
            self._groupby = lambda x: x.properties["equivalent_atoms"]
            disorder_sites.sort(key=self._groupby)

        for orbit, sites in itertools.groupby(
            disorder_sites, key=self._groupby
        ):
            # Can it fit?
            sites = tuple(sites)
            composition = sites[0].species.to_int_dict()
            integer_formula = "".join(
                e + str(a) for e, a in composition.items()
            )
            formula_unit_sum = sum(composition.values())
            if len(sites) % formula_unit_sum:
                raise NeedSupercellError(
                    (
                        f"Can't fit {integer_formula} "
                        f"within {len(sites)} sites "
                        f"(enlarge by {formula_unit_sum/len(sites):.4f}x). "
                        f"If the integer composition ({integer_formula}) is not what you expected, "
                        "please try to increase the precision of the occupancy "
                        "e.g., 1/3 = 0.3333 -> 1/3 = 0.3333333333333."
                    )
                )
            self.disorder_groups[orbit] = sites

            # DMAT
            indices = [x.properties["index"] for x in sites]
            self._group_indices[orbit] = indices
            group_dmat = self._structure.distance_matrix[
                np.ix_(indices, indices)
            ]
            self._group_dmat[orbit] = self.ordinalize(
                group_dmat, atol=self._atol
            )

            # PERM
            coords = [x.frac_coords for x in sites]
            group_perms = []

            for symmop in tqdm.tqdm(
                self._symmops,
                desc="Building symmetry list",
                **const.TQDM_CONF,
                disable=const.DISABLE_PROGRESSBAR,
            ):
                o_coords = symmop.operate_multi(coords)
                atol = const.FIND_COORDS_LIST_ATOL_INIT
                step = const.FIND_COORDS_LIST_ATOL_MUL_STEP_INIT
                is_up = True

                for _ in range(const.FIND_COORDS_LIST_ATTEMPTS):
                    matches = [
                        find_in_coord_list_pbc(coords, coord, atol=atol)
                        for coord in o_coords
                    ]
                    # Sometimes the operated coordinates fell
                    # outside of the symmetry tolerance.
                    # On the other hand, if the tolerance is too low,
                    # one site can be matched with multiple (wrong) sites.
                    #
                    # Robust: Attempt to guess proper tolerance
                    # when one of the above case was detected.
                    # Increment/decrement tolerance with fixed steps,
                    # then binary search if the tolerance becomes too low
                    # after initially too high, and vice versa.
                    # Raise exception if can't find proper tolerance after 10 attempts
                    match_lengths = [len(x) for x in matches]
                    if all(x == 1 for x in match_lengths):
                        break
                    elif any(not x for x in match_lengths):
                        # Too strict
                        if not is_up:
                            step /= 2
                            is_up = True
                        atol *= step
                    elif any(x > 1 for x in match_lengths):
                        # Too loose
                        if is_up:
                            step /= 2
                            is_up = False
                        atol /= step
                else:
                    raise RuntimeError("Failed to build symmetry list.")

                indices = [x[0] for x in matches]
                group_perms.append(indices)
            group_perms = np.array(group_perms)
            self._group_perms[orbit] = group_perms
            row_members = np.sort(group_perms, axis=1)
            if not (row_members == row_members[0]).all():
                raise RuntimeError("Invalid permutation(s) generated.")

            bits = [2**i for i in range(group_perms.shape[1])]
            self._group_bits[orbit] = bits

            try:
                bit_perm = np.array(
                    [[2 ** int(i) for i in row] for row in group_perms],
                    dtype=int,
                )
            except OverflowError:
                bit_perm = np.array(
                    [[2 ** int(i) for i in row] for row in group_perms],
                    dtype=object,
                )
            self._group_bit_perm[orbit] = bit_perm

        # If not explicitly set, turn on caching if there are multiple orbits/colors,
        # and turn off otherwise
        if self.cache is None:
            da = list(self._disorder_amounts().values())
            if not da:  # no disorder
                pass
            else:
                if len(da) > 1 or len(da[0]) > 2:
                    self.cache = True
                else:
                    self.cache = False

        # Letter-related functions
        self._segmenter = self._disorder_amounts().values()
        n_segments = sum(len(x) for x in self._segmenter)
        self._charset = [chr(97 + i) for i in range(n_segments)]

    @staticmethod
    def ordinalize(array, atol=const.DEFAULT_ATOL):
        """
        Ordinalize array elements to the specified absolute tolerance.

        Args:
            array (np.array): The array to be ordinalized.
            atol (float): Absolute tolerance. Defaults to 1e-5.

        Returns:
            np.array: The ordinalized array.
        """
        flat = array.flatten()
        dflat = []
        points = [flat[0]]
        for element in flat:
            equal = 0
            for i, point in enumerate(points):
                dist = np.abs(point - element)
                if dist < atol:
                    equal = i
                    break
            else:
                equal = len(points)
                points.append(element)
                dflat.append(equal)
                continue
            dflat.append(equal)

        return np.array(dflat).reshape(array.shape)

    @property
    def pattern_makers(self):
        """
        PatternMaker instances in key-value pairs with its label as key
        """
        return self._pms

    @pattern_makers.setter
    def pattern_makers(self, pms):
        self._pms = pms

    def _sorted_compositions(self):
        """
        Integer amount of sites substituted to meet desired specification.
        Ordered to minimize the size of the substituted part.
        """
        orbit_compositions = dict()
        for orbit, sites in self.disorder_groups.items():
            ref = sites[0]
            composition = dict(
                sorted(
                    ref.species.inted_composition.items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
            )
            fu = int(sum(composition.values()))
            # Note: __init__ checks makes sure this has no remainder.
            fus = len(sites) // fu

            orbit_compositions[orbit] = {
                e: fus * int(a) for e, a in composition.items()
            }
        return orbit_compositions

    def _disorder_elements(self):
        return {
            orbit: tuple(x.keys())
            for orbit, x in self._sorted_compositions().items()
        }

    def _disorder_amounts(self):
        return {
            orbit: tuple(x.values())
            for orbit, x in self._sorted_compositions().items()
        }

    def make_patterns(self):
        """
        "Raw" pattern of the all disorder sites.

        Returns:
            list: list of list of sites with a distinct species.
        """

        def rscum(iterable):
            """Cumulative sum of an iterable, but skip the last element, and reverse."""
            cum = 0
            cums = []
            for e in iterable[:-1]:
                cum += e
                cums.append(cum)
            for cum in cums[::-1]:
                yield cum

        def nocache_get_pm(subperm, dmat):
            pm = PatternMaker(
                subperm,
                invar=dmat,
                enumerator_collection=self._enumerator_collection,
                t_kind=self._t_kind,
                shuffle=self._shuffle,
                rg=self._rg,
            )
            row_map = pm.get_row_map()
            index_map = pm.get_index_map()
            return pm, row_map, index_map

        def cached_get_pm(subperm, dmat):
            label = PatternMaker.get_label(subperm)
            if label in self._pms:
                pm = self._pms[label]
                row_map, index_map = pm.update_index(subperm)
            else:
                pm = PatternMaker(
                    subperm,
                    invar=dmat,
                    enumerator_collection=self._enumerator_collection,
                    t_kind=self._t_kind,
                    cache=self.cache,
                    shuffle=self._shuffle,
                    rg=self._rg,
                )
                row_map = pm.get_row_map()
                index_map = pm.get_index_map()
                self._pms[label] = pm
            return pm, row_map, index_map

        def maker_recurse_unit(aut, pattern, orbit, amount):
            """PatternMaker aut/pattern generation recursion unit."""
            group_perms = self._group_perms[orbit]
            group_dmat = self._group_dmat[orbit]

            # TODO: Do pattern "plus" at the end.
            subpattern = pattern[-1]
            subperm = group_perms[np.ix_(aut, subpattern)]
            if not self._no_dmat:
                dmat = group_dmat[np.ix_(subpattern, subpattern)]
            else:
                dmat = None

            pm, row_map, index_map = get_pm(subperm, dmat)
            for _aut, _subpattern in pm.ap(amount):
                _aut = np.sort(row_map[_aut])
                _subpattern = index_map[_subpattern]
                yield [aut[_aut], pattern + [_subpattern]]

        def maker_recurse_c(aut, pattern, orbit, chain):
            if len(chain) > 0:
                amount = chain.pop()
                for aut, pattern in maker_recurse_unit(
                    aut, pattern, orbit, amount
                ):
                    _chain = chain.copy()
                    yield from maker_recurse_c(aut, pattern, orbit, _chain)
            else:
                yield aut, pattern

        def maker_recurse_o(aut, pattern, ochain):
            if len(ochain) > 0:
                orbit, sites = ochain.pop()

                chain = list(rscum(self._disorder_amounts()[orbit][::-1]))[
                    ::-1
                ]
                indices = np.arange(len(sites))

                for aut, pattern in maker_recurse_c(
                    aut, pattern + [indices], orbit, chain
                ):
                    # TODO: something cheaper?
                    _ochain = ochain.copy()
                    yield from maker_recurse_o(aut, pattern, _ochain)
            else:
                yield aut, pattern

        logging.info("Making patterns.")
        # Safety kill
        count = self.count()
        if count > const.MAX_IRREDUCIBLE:
            raise TooBigError(f"({count} irreducible expected)")

        # To make the branches shallow.
        if self.cache:
            get_pm = cached_get_pm
        else:
            get_pm = nocache_get_pm

        for aut, pattern in maker_recurse_o(
            np.arange(len(self._symmops)),
            [],
            # TODO: Temporary fix.
            list(self.disorder_groups.items())[::-1],
        ):
            yield (aut, pattern)

    def total_count(self):
        """
        Total number of combinations.
        """
        ocount = (
            multinomial_coeff(x) for x in self._disorder_amounts().values()
        )
        return functools.reduce(lambda x, y: x * y, ocount, 1)

    def count(self):
        """
        Final number of patterns.
        """
        logging.info(
            f"\nCounting unique patterns for {self.structure.formula}"
        )

        if len(self._symmops):
            enumerator = self._enumerator_collection.get(
                # orbit_symmetries,
                list(self._group_perms.values()),
                len(self._symmops),
            )
        else:  # No-permutations-supplied
            enumerator = self._enumerator_collection.get([])
        count = enumerator.count(tuple(self._disorder_amounts().values()))
        return count

    def quantities(self, q, symprec=None):
        """
        Mixed quantities generator.
        Yield tuple of selected quantities.

        Args:
            q: (list) valid options: ("cifwriter", "weight", "letter", "ewald", "structure")
                will return in this order.
        """
        is_c = "cifwriter" in q
        is_s = "structure" in q
        is_w = "weight" in q
        is_l = "letter" in q
        is_e = "ewald" in q

        packet = collections.defaultdict(lambda: None)
        for a, p in self.make_patterns():
            packet.clear()
            if is_c:
                packet["cifwriter"] = self._get_cifwriter(p, symprec)
            if is_s:
                packet["structure"] = self._get_structure(p)
            if is_e:
                packet["ewald"] = self._get_ewald(p, symprec)
            if is_w:
                packet["weight"] = self._get_weight(a)
            if is_l:
                packet["letter"] = self._get_letters(p)
            # filter on packet!
            # naive: if filter(packet["ewald"]): yield else skip
            yield packet

    def letters(self):
        """
        Substitution alphabet representation generator.

        Return:
            dict: (tuple -> np.array) pairs of all pattern letters,
                separated by (specified) orbit.

                Each row of the array is a pattern for the key orbit.
        """
        for _, p in self.make_patterns():
            yield self._get_letters(p)

    def weights(self):
        """
        Pattern weights generator.

        Return:
            list: List of weights of all patterns.
        """
        for a, _ in self.make_patterns():
            yield self._get_weight(a)

    def cifwriters(self, symprec=None):
        """
        Cifwriters generator.
        """
        for _, p in self.make_patterns():
            yield self._get_cifwriter(p, symprec)

    def structure_writers(self, symprec=None):
        """
        Pymatgen Structures generator.
        """
        # This one does not need symprec.
        # Just to keep the signature the same.
        del symprec
        for _, p in self.make_patterns():
            yield self._get_structure(p)

    def ewalds(self, symprec=None):
        """
        Ewald energy generator.
        """
        for _, p in self.make_patterns():
            yield self._get_ewald(p, symprec)

    def _get_letters(self, p):
        """
        (sub. data structure) turn list of list
        into list of some generic symbols from symbols_list.
        Used by `configurations()` and `cifwriter()`.
        Probably not the best idea. Graphical explanation:

        [[0, 1, 2], [0, 2], [0, 1], [1]]  # pattern
        +
        ((2, 1), (1, 1))  # segmenter
        +
        [a, b, c, d]  # symbols_list
        =
        [b, a, b, c, d]  # il
        """
        si = iter(self._charset)
        pi = iter(p)
        il = []
        for se in self._segmenter:
            sei = iter(se)
            word = [next(si)] * sum(se)
            next(pi)
            next(sei)
            # Iterate disoder sites
            for _ in sei:
                letter = next(si)
                # Iterate subpatterns
                for j in next(pi):
                    word[j] = letter
            il.extend(word)
        return "".join(il)

    def _get_weight(self, a):
        return len(self._symmops) // a.size

    def _get_cifwriter(self, p, symprec=None):
        """
        Return cifwriter for the given pattern.

        Args:
            p: Substitution pattern.
        """
        # TODO: Simplify or clarify.
        des = self._disorder_elements()
        orbits = des.keys()
        gis = self._group_indices

        # Build template CifWriter. First run only.
        if self._template_cifwriter is None:
            # In this template, all sites are ordered / single-occupied.
            template_structure = self._structure.copy()
            pi = iter(p)
            for orbit in orbits:
                indices = gis[orbit]
                de = des[orbit]
                for e in de:
                    subpattern = next(pi)
                    for i in subpattern:
                        template_structure.sites[indices[i]].species = e
            cifwriter = CifWriter(template_structure)

            # Use faster CifBlock implementation
            cfkey = cifwriter.ciffile.data.keys()
            cfkey = list(cfkey)[0]
            block = AltCifBlock.from_string(str(cifwriter.ciffile.data[cfkey]))
            cifwriter.ciffile.data[cfkey] = block

            self._template_cifwriter = cifwriter
            self._template_structure = template_structure
        else:
            cifwriter = self._template_cifwriter
            template_structure = self._template_structure
            cfkey = cifwriter.ciffile.data.keys()
            cfkey = list(cfkey)[0]
            block = cifwriter.ciffile.data[cfkey]

        if symprec is None:
            type_symbol = block["_atom_site_type_symbol"].copy()
            label = block["_atom_site_label"].copy()

            # These two blocks are always "1" in the final structure
            block["_atom_site_symmetry_multiplicity"] = ["1"] * len(label)
            block["_atom_site_occupancy"] = ["1.0"] * len(label)

            pi = iter(p)
            for orbit in orbits:
                indices = gis[orbit]
                de = des[orbit]
                for e in de:
                    subpattern = next(pi)
                    for s in subpattern:
                        gi = indices[s]
                        type_symbol[gi] = str(e)
                        label[gi] = f"{e.symbol}{gi}"

            block["_atom_site_type_symbol"] = type_symbol
            block["_atom_site_label"] = label
        else:
            format_str = "{:.8f}"
            latt = template_structure.lattice.matrix
            positions = template_structure.frac_coords

            # this only actually
            cell_specie = list(set(x.species for x in template_structure))
            # Flattened list of species @ disorder sites
            specie = [y for x in des.values() for y in x]
            z_map = [
                cell_specie.index(Composition({specie[j]: 1}))
                for j in range(len(p))
            ]
            zs = [cell_specie.index(x.species) for x in template_structure]

            pi = iter(p)
            zi = iter(z_map)
            for orbit in orbits:
                indices = gis[orbit]
                de = des[orbit]
                for e in de:
                    subpattern = next(pi)
                    z = next(zi)
                    for s in subpattern:
                        gi = indices[s]
                        zs[gi] = z

            space_group_data = spglib.get_symmetry_dataset(
                (latt, positions, zs),
                symprec=self._symprec,
                angle_tolerance=self._angle_tolerance,
            )

            ops = [
                transformation_to_string(rot, trans, delim=", ")
                for rot, trans in zip(
                    space_group_data["rotations"],
                    space_group_data["translations"],
                )
            ]
            u, inv = np.unique(
                space_group_data["equivalent_atoms"], return_inverse=True
            )
            equivalent_indices = [[] for _ in range(len(u))]
            for j, inv in enumerate(inv):
                equivalent_indices[inv].append(j)
            unique_indices = [
                (
                    sorted(
                        j,
                        key=lambda s: tuple(
                            abs(x)
                            for x in template_structure.sites[s].frac_coords
                        ),
                    )[0],
                    len(j),
                )
                for j in equivalent_indices
            ]
            unique_indices = sorted(
                unique_indices,
                key=lambda t: (
                    cell_specie[zs[t[0]]].average_electroneg,  # careful here
                    -t[1],
                    template_structure.sites[t[0]].a,
                    template_structure.sites[t[0]].b,
                    template_structure.sites[t[0]].c,
                ),
            )

            block["_symmetry_space_group_name_H-M"] = space_group_data[
                "international"
            ]
            block["_symmetry_Int_Tables_number"] = space_group_data["number"]
            block["_symmetry_equiv_pos_site_id"] = [
                str(i) for i in range(1, len(ops) + 1)
            ]
            block["_symmetry_equiv_pos_as_xyz"] = ops

            block["_atom_site_type_symbol"] = []
            block["_atom_site_label"] = []
            block["_atom_site_symmetry_multiplicity"] = []
            block["_atom_site_fract_x"] = []
            block["_atom_site_fract_y"] = []
            block["_atom_site_fract_z"] = []
            block["_atom_site_occupancy"] = []
            count = 0
            for j, mult in unique_indices:
                # Careful: The structure species itself is not updated
                sp = cell_specie[zs[j]].elements[0]  # careful here
                site = template_structure.sites[j]
                block["_atom_site_type_symbol"].append(sp.__str__())
                block["_atom_site_label"].append(f"{sp.symbol}{count}")
                block["_atom_site_symmetry_multiplicity"].append(str(mult))
                block["_atom_site_fract_x"].append(format_str.format(site.a))
                block["_atom_site_fract_y"].append(format_str.format(site.b))
                block["_atom_site_fract_z"].append(format_str.format(site.c))
                block["_atom_site_occupancy"].append("1.0")
                count += 1

        return cifwriter

    def _get_structure(self, p):
        """
        Get Pymatgen structure for the given substitution pattern.
        """
        des = self._disorder_elements()
        orbits = des.keys()
        gis = self._group_indices

        structure = self._structure.copy()
        pi = iter(p)
        for orbit in orbits:
            indices = gis[orbit]
            de = des[orbit]
            for e in de:
                subpattern = next(pi)
                for i in subpattern:
                    structure.sites[indices[i]].species = e
        return structure

    def _get_ewald(self, p, symprec):
        """
        Get ewald sums of a structure.
        """
        # TODO: less ad hoc implementation.
        cifwriter = self._get_cifwriter(p, symprec)
        cifparser = CifParser.from_string(str(cifwriter))
        structure = cifparser.get_structures(primitive=False)[0]
        try:
            if not np.isclose(structure.charge, 0.0):
                logging.warn(
                    f"Unit cell is charged: (total charge = {structure.charge})."
                )
            return EwaldSummation(structure).total_energy
        except TypeError as exc:
            raise ValueError(
                "Ewald summation required CIFs with defined oxidation states."
            ) from exc


class PatternMaker:
    """
    Generate permutationally-unique patterns
    """

    __slots__ = (
        "_search",
        "ap",
        "_enumerator_collection",
        "_patterns",
        "_auts",
        "_subobj_ts",
        "_get_subobj_ts",
        "_nix",
        "invar",
        "_perms",
        "_nperm",
        "_column_index",
        "_relabel_index",
        "_row_index",
        "_bits",
        "_bit_perm",
        "label",
        "_cache",
        "_sieve",
        "_get_not_reject_mask",
        "_get_accept_mask",
        "_get_mins",
        "_iterate",
        "_rg",
    )

    def __init__(
        self,
        perm_list,
        invar=None,
        enumerator_collection=None,
        t_kind=const.DEFAULT_T_KIND,
        cache=False,
        shuffle=False,
        rg=None,
    ):
        # Algo select bits
        if invar is None:
            self._search = self._invarless_search
        else:
            self._search = self._invar_search

        if t_kind == "sum":
            self._get_subobj_ts = self._get_sum_subobj_ts
            self._get_not_reject_mask = self._get_not_reject_mask_int
            self._get_accept_mask = self._get_accept_mask_int
            self._get_mins = self._get_mins_int
        elif t_kind == "plsum":
            self._get_subobj_ts = self._get_sum_subobj_ts
            self._get_not_reject_mask = self._get_not_reject_mask_float
            self._get_accept_mask = self._get_accept_mask_float
            self._get_mins = self._get_mins_float

            # Modify invar
            self._sieve = [[3, 3]]
            self._fill_sieve(invar.max() + 1)
            invar = self._logprime(invar)
        elif t_kind == "det":
            self._get_subobj_ts = self._get_det_subobj_ts
            self._get_not_reject_mask = self._get_not_reject_mask_float
            self._get_accept_mask = self._get_accept_mask_float
            self._get_mins = self._get_mins_float

            # Modify invar
            # self._sieve = [[3, 3]]
            # self._fill_sieve(invar.max() + 1)
            # invar = self._logprime(invar)
        else:
            raise RuntimeError(f'Unrecognized T kind "{t_kind}"')

        self._cache = cache
        if cache:
            self.ap = self.cached_ap
        else:
            self.ap = self.nocache_ap

        # Optional random.Random instance.
        if rg is None:
            self._rg = random
        else:
            self._rg = rg

        # This is to avoid an additional branch in the generator itself.
        if shuffle:

            def _iterate(x):
                nonzeros = np.flatnonzero(x)
                self._rg.shuffle(nonzeros)
                for i in nonzeros:
                    yield i

            self._iterate = _iterate
        else:
            self._iterate = np.flatnonzero

        if enumerator_collection is None:
            self._enumerator_collection = PolyaCollection()
        else:
            self._enumerator_collection = enumerator_collection

        (
            column_index,
            relabel_index,
            row_index,
            indexed_perm_list,
        ) = PatternMaker.reindex(perm_list)

        if invar is not None:
            self.invar = invar[np.ix_(column_index, column_index)]
        else:
            self.invar = invar

        self._nperm, self._nix = indexed_perm_list.shape

        # The confusing rotations;
        self._column_index = column_index
        self._relabel_index = relabel_index
        self._row_index = row_index
        self._perms = indexed_perm_list

        # Bits below is for set-like comparison upon arrays.
        # Use pure python's unbounded integer for unlimited array length.
        self._bits = [2**i for i in range(self._nix)]
        # Convert to list to avoid overflow
        # Will change type to object for very large integer
        # Careful: each are still numpy object type...
        try:
            self._bit_perm = np.array(
                [[2 ** int(i) for i in row] for row in indexed_perm_list],
                dtype=int,
            )
        except OverflowError:
            self._bit_perm = np.array(
                [[2 ** int(i) for i in row] for row in indexed_perm_list],
                dtype=object,
            )

        # Containers for search results.
        self._patterns = collections.defaultdict(list)
        self._patterns[0] = [{}]
        self._auts = collections.defaultdict(list)
        self._auts[0] = np.ones((self._nperm,), dtype="bool")
        self._subobj_ts = collections.defaultdict(list)
        self._subobj_ts[0] = [np.array([0])]

        self.label = self._perms.tobytes()

    @staticmethod
    def reindex(perm_list):
        """Standardize symmetry ordering for reuse (rough)"""
        perm_list = perm_list.copy()

        # Column sort
        stab_map = perm_list == perm_list[0]
        # TODO: column_index is not much used...
        column_index = np.lexsort(stab_map)
        perm_list = perm_list[:, column_index]

        # TODO: More intuitive if we do this first, then the previous one.
        # Relabel to match column position
        relabel_index = perm_list[0]
        relabel_element = np.vectorize(
            {s: i for i, s in enumerate(relabel_index)}.get
        )
        try:
            perm_list = relabel_element(perm_list)
        except TypeError as exc:
            raise ValueError(
                f"\n{perm_list}\n" "Rows must have same elements."
            ) from exc

        # Row sort
        row_index = np.lexsort(perm_list.T)
        relabeled_perm_list = perm_list[row_index]

        return column_index, relabel_index, row_index, relabeled_perm_list

    def update_index(self, perm_list):
        """
        Update map for the output pattern.

        This is used to update internal/external map in case
        of similar permutations.

        Raises:
            ValueError: If the internal permutation is unrelated to the
                supplied one.
        """
        (
            column_index,
            relabel_index,
            row_index,
            relabeled_perm_list,
        ) = PatternMaker.reindex(perm_list)

        if relabeled_perm_list.tobytes() != self._perms.tobytes():
            raise ValueError(
                "The supplied permutation does not match the internal one.",
            )

        # Update values.
        self._column_index = column_index
        self._relabel_index = relabel_index
        self._row_index = row_index
        return row_index, relabel_index

    @staticmethod
    def get_label(perm_list):
        """
        Rough equivalence. Byte string from a standardized
        representation of the given permutation
        """
        _, _, _, relabeled_perm_list = PatternMaker.reindex(perm_list)
        return relabeled_perm_list.tobytes()

    def get_row_map(self):
        """
        Map between internal order of permutation vs. input permutation
        """
        return self._row_index

    def get_index_map(self):
        """
        Map between internal representation of elements vs. input
        """
        return self._relabel_index

    # TODO: cleanup
    def get_column_map(self):
        """
        Internal column shuffle.
        """
        return self._column_index

    def cached_ap(self, n):
        """
        Get patterns and automorphisms for the specified replacement amount (cached)
        """
        # Patterns are symmetrical
        _n = min(n, self._nix - n)
        # if not _n in self._patterns:
        if n not in self._patterns:
            if _n in self._patterns:
                inverter = np.arange(self._nix)
                self._auts[n] = self._auts[_n]
                self._patterns[n] = [
                    np.setdiff1d(inverter, p) for p in self._patterns[_n]
                ]
            else:
                starts = [i for i in self._patterns.keys() if i < _n]
                if not starts:
                    start = 0
                else:
                    start = max(starts)

                # "Mirror" patterns
                if _n != n:
                    inverter = np.arange(self._nix)
                    ap = [
                        (a, np.setdiff1d(inverter, p))
                        for a, p in self._search(start=start, stop=_n)
                    ]
                else:
                    ap = [
                        (a, p) for a, p in self._search(start=start, stop=_n)
                    ]
                self._auts[n], self._patterns[n] = zip(*ap)

        for a, p in zip(self._auts[n], self._patterns[n]):
            yield a, p

    def nocache_ap(self, n):
        """
        Get patterns and automorphisms for the specified replacement amount
        """
        # Patterns are symmetrical
        _n = min(n, self._nix - n)
        if _n != n:
            inverter = np.arange(self._nix)
            for a, p in self._search(stop=_n):
                rp = np.setdiff1d(inverter, p)
                yield a, rp
        else:
            for a, p in self._search(stop=_n):
                yield a, p

    def _fill_sieve(self, n):
        """
        Incremental sieve of eratosthenes to make first n prime numbers
        """
        i = self._sieve[-1][0]
        while len(self._sieve) < (n - 1):
            i += 2  # no even primes except 2
            stop = int(math.sqrt(i))
            for tup in self._sieve:
                prime, mult = tup
                if prime > stop:
                    self._sieve.append([i, i])
                    break
                while mult < i:
                    mult += prime
                if mult == i:
                    break
            else:
                self._sieve.append([i, i])

    @functools.lru_cache()
    def _primes_list(self):
        return [2] + [x[0] for x in self._sieve]

    def _logprime(self, array):
        primed = [self._primes_list()[x] for x in array.flatten()]
        return np.log(np.array(primed).reshape(array.shape))

    def _invarless_search(self, start=0, stop=None):
        """
        Alternative definition of canonical parent without invariant (much slower).

        Args:
            stop (int): Maximum substitution size before stop
        """
        if stop is None:
            stop = self._nix // 2
        stop = min(stop, self._nix - stop)

        enumerator = self._enumerator_collection.get([self._perms])
        n_pred = enumerator.count(((stop, self._nix - stop),))

        # Progress bar.
        pbar = tqdm.tqdm(
            total=n_pred,
            desc=f"Making patterns ({stop}/{self._nix})",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        )

        stack = []
        if not start:
            # Fill first stack
            noniso_orbits = np.unique(self._perms.min(axis=0))
            for _i in noniso_orbits:
                pattern = np.array([_i])
                bitsum = self._bits[_i]
                patch = np.zeros(self._nix, dtype=int)
                patch[_i] = 1
                bs = self._bit_perm.dot(patch)
                o_bitsums = self._bit_perm.dot(patch)
                aut = np.flatnonzero(o_bitsums == bitsum)
                stack.append((pattern, aut, bs))
        else:
            patterns = self._patterns[start]
            auts = self._auts[start]
            for pattern, aut in zip(patterns, auts):
                bs = self._bit_perm[:, pattern].sum(axis=1)
                stack.append((pattern, aut, bs))

        while stack:
            pattern, aut, pbs = stack.pop()
            if pattern.size == stop:
                yield aut, pattern
                pbar.update()
                continue
            # Tree is expanded by adding un-added sites
            leaf_mask = np.ones(self._nix, dtype="bool")
            leaf_mask[pattern] = False
            leaf_array = np.flatnonzero(leaf_mask)

            # Discard symmetry duplicates from the remaining leaves
            if aut.size > 1 and leaf_mask.sum() > 1:
                leaf_reps = self._perms[np.ix_(aut, leaf_array)].min(axis=0)
                leaf_indices = np.searchsorted(leaf_array, leaf_reps)
                uniq_mask = np.zeros(leaf_array.shape, dtype="bool")
                uniq_mask[leaf_indices] = True
            else:
                uniq_mask = np.ones(leaf_array.shape, dtype="bool")

            # Insertion location
            loci = pattern.searchsorted(leaf_array)

            for i in self._iterate(uniq_mask):
                x = leaf_array[i]
                j = loci[i]
                _pbs = self._bit_perm[:, x] + pbs

                _i = np.concatenate((pattern[:j], [x], pattern[j:]))
                _aut = np.flatnonzero(_pbs == _pbs[-1])

                # Compute canonical parent.
                m = np.argmin(_pbs)
                can_pattern = self._perms[m, _i]
                discard_i = np.where(can_pattern == can_pattern.max())[0]

                # If the new site is discarded then tree parent == canonical parent.
                if j == discard_i:
                    stack.append((_i, _aut, _pbs))
                # Check if tree parent is related to canonical parent.
                elif _i[discard_i] in self._perms[_aut, x]:
                    stack.append((_i, _aut, _pbs))
        pbar.close()
        tqdm.tqdm(disable=const.DISABLE_PROGRESSBAR).write("Done.")

    def _get_sum_subobj_ts(self, pattern, leaf_array, subobj_ts):
        new_rows = self.invar[np.ix_(leaf_array, pattern)]
        new_rows_sums = new_rows.sum(axis=1, keepdims=True)
        part_new_row_sums = subobj_ts + new_rows

        return np.concatenate((part_new_row_sums, new_rows_sums), axis=1)

    @functools.lru_cache(None)
    def _get_det(self, pattern):
        pattern_invar = self.invar[np.ix_(pattern, pattern)]
        return np.abs(np.linalg.det(np.atleast_2d(pattern_invar)))
        # return np.linalg.det(np.atleast_2d(pattern_invar))

    def _get_det_subobj_ts(self, pattern, leaf_array, subobj_ts):
        tiled_pattern = np.tile(pattern, (leaf_array.size, 1))
        leaves = np.column_stack((tiled_pattern, leaf_array))

        subobj_ts = []
        for leaf in leaves:
            # NOTE: Be careful with the order of these combinations...
            leaf_subobjs = list(itertools.combinations(leaf, leaf.size - 1))
            leaf_subobjs_ts = [self._get_det(x) for x in leaf_subobjs[::-1]]
            subobj_ts.append(leaf_subobjs_ts)

        return np.array(subobj_ts)

    @staticmethod
    def _get_not_reject_mask_int(delta_t):
        return ~(delta_t < 0).any(axis=1)

    def _get_not_reject_mask_float(self, delta_t):
        # return np.ones(delta_t.shape[0], dtype=bool)
        return ~((delta_t < 0.0) & ~np.isclose(delta_t, 0.0)).any(axis=1)

    @staticmethod
    def _get_accept_mask_int(delta_t):
        return (delta_t > 0).all(axis=1)

    def _get_accept_mask_float(self, delta_t):
        # return np.zeros(delta_t.shape[0], dtype=bool)
        return ((delta_t > 0.0) & ~np.isclose(delta_t, 0.0)).all(axis=1)

    @staticmethod
    def _get_mins_int(subobj_ts):
        return np.flatnonzero(subobj_ts == subobj_ts.min())

    @staticmethod
    def _get_mins_float(subobj_ts):
        return np.flatnonzero(np.isclose(subobj_ts, subobj_ts.min()))

    def _invar_search(self, start=0, stop=None):
        """
        Combines invar and lexmax to determine canonical parent.

        Args:
            start: Start seach at this depth.
            stop: Stop search at this depth.
        """
        if stop is None:
            stop = self._nix // 2
        stop = min(stop, self._nix - stop)

        enumerator = self._enumerator_collection.get([self._perms])
        n_pred = enumerator.count(((stop, self._nix - stop),))

        # Progress bar.
        pbar = tqdm.tqdm(
            total=n_pred,
            desc=f"Making patterns ({stop}/{self._nix})",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        )

        stack = []
        if not start:
            # Populate the first layer.
            noniso_orbits = np.unique(self._perms.min(axis=0))
            for _i in noniso_orbits:
                pattern = np.array([_i])
                # TODO: refine
                bitsum = self._bits[_i]
                patch = np.zeros(self._nix, dtype=int)
                patch[_i] = 1
                bs = self._bit_perm.dot(patch)
                aut = np.flatnonzero(bs == bitsum)
                subobj_ts = np.array([0])
                stack.append((subobj_ts, pattern, aut, bs))
        else:
            patterns = self._patterns[start]
            auts = self._auts[start]
            for pattern, aut in zip(patterns, auts):
                bs = self._bit_perm[:, pattern].sum(axis=1)
                subobj_ts = self.invar[np.ix_(pattern, pattern)].sum(axis=1)
                stack.append((subobj_ts, pattern, aut, bs))

        while stack:
            subobj_ts, pattern, aut, pbs = stack.pop()
            if pattern.size == stop:
                yield aut, pattern
                pbar.update()
                continue
            # Invert index
            leaf_mask = np.ones(self._nix, dtype="bool")
            leaf_mask[pattern] = False
            leaf_array = np.flatnonzero(leaf_mask)

            # Calculate subobject Ts for all leaves
            leaf_subobj_ts = self._get_subobj_ts(
                pattern, leaf_array, subobj_ts
            )

            # Reject all leaves where any T is smaller than the new row's T
            delta_t = leaf_subobj_ts[:, :-1] - leaf_subobj_ts[:, -1:]
            not_reject_mask = self._get_not_reject_mask(delta_t)
            if not not_reject_mask.any():
                continue

            # Discard symmetry duplicates from the remaining leaves
            if aut.size > 1 and not_reject_mask.sum() > 1:
                not_reject_leaf = leaf_array[not_reject_mask]
                leaf_reps = self._perms[np.ix_(aut, not_reject_leaf)].min(
                    axis=0
                )
                leaf_indices = leaf_array.searchsorted(leaf_reps)
                uniq_mask = np.zeros(leaf_array.shape, dtype="bool")
                uniq_mask[leaf_indices] = True
            else:
                uniq_mask = not_reject_mask

            # Insertion location TODO: without this possible?
            loci = pattern.searchsorted(leaf_array)

            accept_mask = self._get_accept_mask(delta_t)
            accept_mask &= uniq_mask
            accepts = leaf_array[accept_mask]

            for i in self._iterate(uniq_mask):
                x = leaf_array[i]
                j = loci[i]
                _pbs = self._bit_perm[:, x] + pbs

                _subobj_ts = leaf_subobj_ts[i]
                _subobj_ts[j:] = np.concatenate(
                    (_subobj_ts[-1:], _subobj_ts[j:-1])
                )

                _i = np.concatenate((pattern[:j], [x], pattern[j:]))
                # NOTE: just in case I fail to consistently sort perm
                # bitsum = sum([self._bits[y] for y in _i])
                # _aut = np.flatnonzero(_pbs == bitsum)
                _aut = np.flatnonzero(_pbs == _pbs[-1])

                if x in accepts:
                    stack.append((_subobj_ts, _i, _aut, _pbs))
                    continue

                # Compute canonical parent.
                ts_min_i = self._get_mins(_subobj_ts)
                _sub = _i[ts_min_i]
                m = np.argmin(_pbs)
                can_pattern = self._perms[m, _sub]
                discard_i = ts_min_i[can_pattern == can_pattern.max()]

                # If the new site is discarded then tree parent == canonical parent.
                if j == discard_i:
                    stack.append((_subobj_ts, _i, _aut, _pbs))
                # Check if tree parent is related to canonical parent.
                elif _i[discard_i] in self._perms[_aut, x]:
                    stack.append((_subobj_ts, _i, _aut, _pbs))
        pbar.close()
        tqdm.tqdm(disable=const.DISABLE_PROGRESSBAR).write("Done.")
        # self._gen_flag[stop] = True


class Polya:
    """
    Perform operations related (weighed) Polya Enumeration Theorem.

    Includes enumeration (count) and calculation of cycle index
    and configuration generation function.

    Args:
        perms_list (list, tuple): List of permutations as (np.array).
            When specifying multiple substitution involving separated orbits,
            each of the separated group of sites should be assigned their own
            letter in the cycle index to avoid mixing of the substituted species.

            Rows of the arrays corresponds to a single permutation,
            where each element e permutes column c to e (c -> e).
            All arrays should contain the same number of rows.
    """

    def __init__(
        self,
        perm_list,
        group_size=None,
    ):
        # Used only for dividing the final weight
        if group_size is not None:
            self.group_size = group_size
        else:
            perm_sizes = set(x.shape[0] for x in perm_list)
            if not perm_sizes:
                self.group_size = 1
                logging.warning("No permutations supplied to enumerator.")
            else:
                self.group_size = max(perm_sizes)
        self._perm_list = perm_list
        self._expand_cache = dict()

    def label(self):
        """Use ci as label"""
        return sum(self.sym_ci().values()) / self.group_size

    @functools.lru_cache(None)
    def ci(self):
        """
        Returns the cycle index (as Counters of cycle lengths per permutation).

        Returns:
            dict: Dictionary of cycle index with integer index
                and Counters as a key-value pair.
        """
        # Use indexed base.
        cycle_index = collections.defaultdict(list)
        for permutations in self._perm_list:
            try:
                index_map = {s: i for i, s in enumerate(permutations[0])}
            except IndexError as e:
                logging.error("========================================")
                logging.error(self._perm_list)
                logging.error(permutations)
                logging.error("========================================")
                raise e
            # Iterate by row.
            for i, permutation in enumerate(permutations):
                # Build cycles.
                cycles = []
                visited = []
                for head in permutation:
                    if head in visited:
                        continue
                    cycle = []
                    while head not in cycle:
                        cycle.append(head)
                        visited.append(head)
                        try:
                            head = permutation[index_map[head]]
                        except (KeyError, IndexError) as exc:
                            logging.error("=============================")
                            logging.error(f"HEAD: {head}")
                            logging.error(f"IMAP: {index_map}")
                            logging.error(f"P: {permutation}")
                            logging.error(f"BP: {permutations}")
                            raise RuntimeError(
                                "Check permutation list."
                            ) from exc
                    cycles.append(cycle)
                counter = collections.Counter(len(cycle) for cycle in cycles)
                cycle_index[i].append(counter)
        return cycle_index

    def sym_ci(self):
        """
        Returns the symbolic cycle index.
        Uses sympy.IndexedBase to represent subscript.

        Conventionally, each term is in the form of (x_n**m),
        with (n) denoting cycle length and (m) number of
        cycles with such length for a permutation.

        Returns:
            dict: Dictionary of cycle index with integer index
                and sympy equations as a key-value pair.
        """
        cycle_index = self.ci()
        sym_cycle_index = dict()
        for j, permutations in enumerate(self._perm_list):
            symbol = sympy.IndexedBase(chr(97 + j))
            for i, _ in enumerate(permutations):
                sym_cycle_index.setdefault(i, sympy.Integer(1))
                for length, n in cycle_index[i][j].items():
                    sym_cycle_index[i] *= symbol[length] ** n
        return sym_cycle_index

    @functools.lru_cache(None)
    def count(self, amt_tuple):
        """
        Count unique pattern given.

        Args:
            amt_list (tuple): List of list, with the
                numbers in the inner list are coefficients
                of a chemical formula A_xB_y.

                In order to find results, the formula needs to be
                in the integer amount of the number of actual sites
                in the structure.

                Each elements addresses separate orbits.

        Returns:
            int: Number of unique patterns for the given amount.
        """

        def combine(a, b):
            assert len(a.shape) == 2 and len(b.shape) == 2
            ra, ca = a.shape
            rb, cb = b.shape
            assert ca == cb
            return (a[:, None] + b[None, :]).reshape(ra * rb, ca)

        # Should be cached too somehow
        # @functools.lru_cache(None)
        def exmul(arrays):
            """
            Calculate exponent coefficients "multiplied".
            """
            return functools.reduce(combine, arrays)

        logging.info(f"Counting unique patterns for {amt_tuple}.")

        # Padding
        joint_coeffs = np.array(sum(amt_tuple, ()))
        coeff_sums = [len(x) for x in amt_tuple]
        pads = [
            (sum(coeff_sums[:i]), sum(coeff_sums[i + 1 :]))
            for i in range(len(amt_tuple))
        ]

        o_counts = []
        for o_cycles in self.ci().values():
            o_parts = [
                [aR_array(cnum, len(color)) for cnum in cycles.values()]
                for cycles, color in zip(o_cycles, amt_tuple)
            ]

            # Exponent values of each variables
            exps = [
                np.pad(clen * partition, [(0, 0), pad])
                for cycles, pad, partitions in zip(o_cycles, pads, o_parts)
                for clen, partition in zip(cycles, partitions)
            ]
            config_shape = [x.shape[0] for x in exps]

            mul_exps = exmul(exps)
            # Find matching coefficient; if none then 0
            match = np.flatnonzero((mul_exps == joint_coeffs).all(axis=1))
            match_i = np.array(np.unravel_index(match, config_shape)).T

            # TODO: Basically flatten; can be better written
            f_parts = [part for parts in o_parts for part in parts]
            counts = [
                functools.reduce(
                    lambda x, y: x * y,
                    [
                        multinomial_coeff(tuple(p[j]))
                        for p, j in zip(f_parts, i)
                    ],
                )
                for i in match_i
            ]
            if not counts:
                counts = [0]
            o_counts.append(sum(counts))
        # Special case: no cycles means 1 unique object
        if not o_counts:
            return 1
        return int(sum(o_counts) / self.group_size)


class PolyaCollection:
    """
    Collection of Polya objects. This is useful because identical
    instances are often recreated, with heavy sympy.expand() operations.
    TODO: No longer important so eventually remove.
    """

    def __init__(self):
        self._collection = dict()

    def __getitem__(self, key):
        return self._collection[key]

    def __setitem__(self, key, value):
        self._collection[key] = value

    def get(self, permutation_list, group_size=None):
        """
        Get appropriate Polya instance for the given parameters
        """
        polya = Polya(permutation_list, group_size)
        label = polya.label()
        if label in self._collection:
            return self._collection[label]
        else:
            self._collection[label] = polya
            return polya
