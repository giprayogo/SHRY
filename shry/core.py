# -*- coding: utf-8 -*-
# pylint: disable=unused-import, logging-fstring-interpolation, invalid-name, too-many-lines

# information
__author__ = "Genki Prayogo, and Kosuke Nakano"
__copyright__ = "Copyright (c) 2021-, The SHRY Project"
__credits__ = ["Genki Prayogo", "Kosuke Nakano"]

__license__ = "MIT"
__version__ = "1.0.3"
__maintainer__ = "Genki Prayogo"
__email__ = "g.prayogo@icloud.com"
__date__ = "15. Nov. 2021"
__status__ = "Production"

"""
Core operations, pattern generation, etc.
"""

import collections
import functools
import itertools
import logging
import math
import sys
from pprint import pprint
from typing import OrderedDict, Tuple

import numpy as np
import spglib
import sympy
import tqdm
from memory_profiler import profile
from monty.fractions import gcd_float
from pymatgen.core.composition import Composition, reduce_formula
from pymatgen.core.operations import SymmOp
from pymatgen.io.cif import CifBlock, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.symmetry.structure import SymmetrizedStructure
from pymatgen.util.coord import find_in_coord_list_pbc
from pymatgen.util.string import transformation_to_string
from sympy.tensor.indexed import Indexed
from tabulate import tabulate

from . import const

np.seterr(all="raise")
np.set_printoptions(linewidth=const.LINEWIDTH, threshold=16)
np.set_printoptions(linewidth=1000, threshold=sys.maxsize)


def get_integer_formula_and_factor(
    self, max_denominator: int = 10000, iupac_ordering: bool = False
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
    _, factor = self.get_integer_formula_and_factor()
    int_dict = {e: int(a) for e, a in (self / factor).as_dict().items()}
    # be safe: Composition groups together different ox states which is not ideal...
    if not all(
        np.isclose(x * factor, y)
        for x, y in zip(int_dict.values(), self.as_dict().values())
    ):
        raise ValueError("Composition is not rational!")
    return int_dict


@property
def inted_composition(self):
    """
    Return Composition instance with integer formula
    """
    _, factor = self.get_integer_formula_and_factor()
    int_comp = self / factor

    # be safe
    int_dict = {e: int(a) for e, a in int_comp.as_dict().items()}
    if not all(
        np.isclose(x * factor, y)
        for x, y in zip(int_dict.values(), self.as_dict().values())
    ):
        raise ValueError("Composition is not rational!")

    return int_comp


Composition.to_int_dict = to_int_dict
Composition.get_integer_formula_and_factor = get_integer_formula_and_factor
Composition.inted_composition = inted_composition


class PatchedSymmetrizedStructure(SymmetrizedStructure):
    """
    Fixed site_properties display
    """

    def __str__(self):
        outs = [
            "SymmetrizedStructure",
            "Full Formula ({s})".format(s=self.composition.formula),
            "Reduced Formula: {}".format(self.composition.reduced_formula),
        ]

        def to_s(x):
            return "%0.6f" % x

        outs.append(
            "abc   : " + " ".join([to_s(i).rjust(10) for i in self.lattice.abc])
        )
        outs.append(
            "angles: " + " ".join([to_s(i).rjust(10) for i in self.lattice.angles])
        )
        if self._charge:
            if self._charge >= 0:
                outs.append("Overall Charge: +{}".format(self._charge))
            else:
                outs.append("Overall Charge: -{}".format(self._charge))
        outs.append("Sites ({i})".format(i=len(self)))
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
                    line[i : i + self.maxlen] for i in range(0, len(line), self.maxlen)
                ]
                s.extend(sublines)
        return s

    def __str__(self):
        """
        Returns the cif string for the data block
        """
        s = ["data_{}".format(self.header)]
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
                        s.append("{}   {}".format(k, v))
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
                        self.data[_k] = list(map(self._format_field, self.data[_k]))
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
                    s.append("{}   {}".format(k, v))
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


class Substitutor:
    """
    Makes unique ordered Structure(s) for the given disorder Structure.

    Since the PatternMaker instances are saved,
    it is advantageous to reuse the same Substitutor
    for symmetrically similar structures.

    This includes different subsitute concentrations
    from the same basic structure, among others.

    Args:
        structure (pymatgen.Structure): Input structure containing disorder sites.
            If all sites are ordered then the output structure is exactly the same.
        symprec (float): Precision used for symmetry analysis.
        groupby (function): Function to group disordered sites.
            Defaults to lambda x: x.properties["_atom_site_label"], with the x loops
            over all PeriodicSites within the Structure.

            The proper function depends on what you want; do you want all sites of
            a certain species or crystallographic orbit to be subsituted together?

            Fallback to crystallographic orbit when failed.
        sample (int): Randomly choose some of the generated structures.
    """

    __slots__ = (
        "_symprec",
        "_angle_tolerance",
        "_groupby",
        "_symmops",
        "_patterns",
        "_pattern_automorphisms",
        "_pattern_makers",
        "_enumerator_collection",
        "disorder_groups",
        "_group_dmat",
        "_group_perms",
        "_group_indices",
        "_group_bits",
        "_group_bit_perm",
        "_structure",
        "_sample",
        "_no_dmat",
        "_t_kind",
        "_made_patterns",
    )

    def __init__(
        self,
        structure,
        symprec=const.DEFAULT_SYMPREC,
        angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        groupby=None,
        sample=None,
        no_dmat=const.DEFAULT_NO_DMAT,
        t_kind=const.DEFAULT_T_KIND,
    ):
        self._symprec = symprec
        self._angle_tolerance = angle_tolerance
        self._no_dmat = no_dmat
        self._t_kind = t_kind
        if groupby is None:
            self._groupby = lambda x: x.properties["_atom_site_label"]

        self._symmops = None
        self._patterns = []
        self._pattern_automorphisms = []
        self._pattern_makers = dict()
        self._enumerator_collection = PolyaCollection()

        self.disorder_groups = dict()
        self._group_dmat = dict()
        self._group_perms = dict()
        self._group_indices = dict()
        self._group_bits = dict()
        self._group_bit_perm = dict()

        self.structure = structure
        self.sampled_indices = sample
        self._made_patterns = False

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

    @property
    def sampled_indices(
        self,
    ):
        """
        Used when sampling the patterns
        """
        return self._sample

    @sampled_indices.setter
    def sampled_indices(self, sample):
        try:
            indices = np.arange(self.count())
        except ValueError as e:
            raise TooBigError(f"({self.count()} irreducible structures)") from e
        if sample is None:
            self._sample = indices
        else:
            logging.info(f"Sampling {sample} of {self.count()} structures.")
            try:
                self._sample = np.random.choice(indices, size=sample, replace=False)
            except ValueError:
                logging.warning(
                    "WARNING: Selected sample size exceeded the number of patterns. Using all."
                )
                self._sample = indices

        self.configurations.cache_clear()
        self.weights.cache_clear()

    @structure.setter
    def structure(self, structure, sample=None):
        logging.info("\nSetting Substitutor with Structure")
        logging.info(f"{structure}")
        self._structure = structure.copy()

        sga = PatchedSpacegroupAnalyzer(
            self._structure,
            symprec=self._symprec,
            angle_tolerance=self._angle_tolerance,
        )
        self._symmops = sga.get_symmetry_operations()

        logging.info(f"Space group: {sga.get_hall()} ({sga.get_space_group_number()})")
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
                # Ad hoc fix: if occupancy is less than 1,
                # Stop the program
                # TODO: Automatic handling
                if not np.isclose(site.species.num_atoms, 1):
                    raise RuntimeError("Please fill vacancy sites with pseudo atoms")
        if not disorder_sites:
            logging.warning("No disorder sites found within the Structure.")
            return

        try:
            disorder_sites.sort(key=self._groupby)
        except KeyError:
            # Fallback to crystallographic orbit.
            self._groupby = lambda x: x.properties["equivalent_atoms"]
            disorder_sites.sort(key=self._groupby)

        for orbit, sites in itertools.groupby(disorder_sites, key=self._groupby):
            # Can it fit?
            sites = tuple(sites)
            composition = sites[0].species.to_int_dict()
            integer_formula = "".join(e + str(a) for e, a in composition.items())
            formula_unit_sum = sum(composition.values())
            if len(sites) % formula_unit_sum:
                raise NeedSupercellError(
                    (
                        f"Can't fit {integer_formula} "
                        f"within {len(sites)} sites "
                        f"(enlarge by {formula_unit_sum/len(sites):.4f}x)."
                    )
                )
            self.disorder_groups[orbit] = sites

            # DMAT
            indices = [x.properties["index"] for x in sites]
            self._group_indices[orbit] = indices
            group_dmat = self._structure.distance_matrix[np.ix_(indices, indices)]
            self._group_dmat[orbit] = self.ordinalize(group_dmat)

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

            bits = [2 ** i for i in range(group_perms.shape[1])]
            self._group_bits[orbit] = bits

            try:
                bit_perm = np.array(
                    [[2 ** int(i) for i in row] for row in group_perms], dtype=int
                )
            except OverflowError:
                bit_perm = np.array(
                    [[2 ** int(i) for i in row] for row in group_perms], dtype=object
                )
            self._group_bit_perm[orbit] = bit_perm

        # Clear previous caches.
        self.configurations.cache_clear()
        self.weights.cache_clear()
        self.count.cache_clear()

        self._made_patterns = False
        self.sampled_indices = sample

    @staticmethod
    def ordinalize(array, atol=1e-8):
        """
        Ordinalize array elements to the specified absolute tolerance.

        Args:
            array (np.array): The array to be ordinalized.
            atol (float): Absolute tolerance. Defaults to 1e-8.

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
        return self._pattern_makers

    @pattern_makers.setter
    def pattern_makers(self, pattern_makers):
        self._pattern_makers = pattern_makers

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

            # TODO: stricter ordering.
            orbit_compositions[orbit] = {
                e: fus * int(a) for e, a in composition.items()
            }
        return orbit_compositions

    def _disorder_elements(self):
        return {
            orbit: tuple(x.keys()) for orbit, x in self._sorted_compositions().items()
        }

    def _disorder_amounts(self):
        return {
            orbit: tuple(x.values()) for orbit, x in self._sorted_compositions().items()
        }

    def make_patterns(self):
        """
        "Raw" pattern of the all disorder sites.

        Returns:
            list: list of list of sites with a distinct species.
        """
        if self._made_patterns is True:
            return
        self._made_patterns = True
        logging.info("Making patterns.")

        def rscum(iterable):
            """Cumulative sum of an iterable, but skip the last element, and reverse."""
            cum = 0
            cums = []
            for e in iterable[:-1]:
                cum += e
                cums.append(cum)
            for cum in cums[::-1]:
                yield cum

        in_stack = []
        out_stack = []

        # Initial values
        aut = np.arange(len(self._symmops))

        # Structure: auts, pattern (growing list)
        in_stack.append([aut, []])
        ran = False
        for orbit, sites in self.disorder_groups.items():
            ran = True
            logging.info(f"Making pattern for {orbit}")

            # Reverse to minimize sub. amount.
            chain = list(rscum(self._disorder_amounts()[orbit][::-1]))

            # Feed the new orbit into the stack.
            # (But read the aut from patterns from the previous orbit)
            # x[1] is the pattern
            indices = np.arange(len(sites))
            for x in in_stack:
                x[1].append(indices)

            group_perms = self._group_perms[orbit]
            group_dmat = self._group_dmat[orbit]

            # Intra-orbit chaining.
            for amount in chain:
                logging.info(f"Making pattern for {amount}/{chain}")
                # Progress bar.
                pbar = tqdm.tqdm(
                    total=len(in_stack),
                    desc="Progress",
                    **const.TQDM_CONF,
                    disable=const.DISABLE_PROGRESSBAR,
                )

                for aut, pattern in in_stack:
                    # Operate on the _last_ subpattern, except for the first one
                    subpattern = pattern[-1]
                    subperm = group_perms[np.ix_(aut, subpattern)]
                    # Do we need?
                    # subperm, index = np.unique(subperm, axis=0, return_index=True)
                    if not self._no_dmat:
                        dmat = group_dmat[np.ix_(subpattern, subpattern)]
                    else:
                        dmat = None

                    label = PatternMaker.get_label(subperm)
                    if label in self._pattern_makers:
                        maker = self._pattern_makers[label]
                        maker.update_index(subperm)
                    else:
                        maker = PatternMaker(
                            subperm,
                            invar=dmat,
                            enumerator_collection=self._enumerator_collection,
                            t_kind=self._t_kind,
                        )
                        self._pattern_makers[label] = maker
                    patterns = maker.patterns(amount)
                    auts = maker.auts(amount)
                    for _aut, _subpattern in zip(auts, patterns):

                        _pattern = pattern + [_subpattern]
                        # out_stack.append([index[_aut], _pattern])
                        out_stack.append([aut[_aut], _pattern])
                    pbar.update()
                pbar.close()
                in_stack = out_stack
                out_stack = []
        if ran:
            self._pattern_automorphisms = [x[0] for x in in_stack]
            self._patterns = [x[1] for x in in_stack]
            n_generated = len(in_stack)
            n_expected = self.count()
            if n_generated != n_expected:
                raise RuntimeError(
                    "Mismatch between generated and predicted "
                    f"number of structures ({n_generated}/{n_expected})"
                )

    @functools.lru_cache()
    def configurations(self):
        """
        Return generic alphabetical letters for denoting the substitutions.

        Return:
            dict: (tuple -> np.array) pairs of all pattern letters,
                separated by (specified) orbit.

                Each row of the array is a pattern for the key orbit.
        """
        logging.info("\nBuilding configuration letters.")

        # Use **amount** so that I can initialize the letters
        # for each orbit
        das = self._disorder_amounts().values()
        n_segments = sum(len(x) for x in das)
        letters = [chr(97 + i) for i in range(n_segments)]

        configs = []
        for i in self.sampled_indices:
            sentence = self._ll2il(das, self._patterns[i], letters)
            configs.append("".join(sentence))
        return configs

    @staticmethod
    def _ll2il(segmenter, pattern, symbols_list):
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
        si = iter(symbols_list)
        pi = iter(pattern)
        il = []
        for se in segmenter:
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
        return il

    @functools.lru_cache()
    def weights(self):
        """
        Pattern weights.

        Return:
            list: List of weights of all patterns.
        """
        logging.info("\nObtaining pattern weights.")
        space_group_size = len(self._symmops)
        weights = []
        for i in self.sampled_indices:
            weights.append(space_group_size // self._pattern_automorphisms[i].size)
        return weights

    @functools.lru_cache()
    def count(self):
        """
        Final number of expected patterns
        """
        logging.info("\nCounting total unique patterns for Structure.")

        if len(self._symmops):
            enumerator = self._enumerator_collection.get(
                # orbit_symmetries,
                list(self._group_perms.values()),
                len(self._symmops),
            )
        else:  # No-permutations-supplied
            enumerator = self._enumerator_collection.get([])
        count = enumerator.count(tuple(self._disorder_amounts().values()))

        # Arbitrary safe limit
        if count > const.MAX_IRREDUCIBLE:
            raise TooBigError(f"({count} irreducible expected)")
        return count

    def cifwriters(self, symprec=None):
        """
        Generator for pymatgen.CifWriter of the final output.

        Returns:
            generator: A generator of pymatgen.CifWriter. Each iteration
                returned a symmetrically unique CifWriter for the given
                input disorder.
        """
        logging.info("\nCreating CifWriter instances.")
        template_structure = self._structure.copy()
        try:
            template_pattern = self._patterns[self.sampled_indices[0]]
        except IndexError:
            raise RuntimeError("Patterns have not been generated.")

        # Build template CifWriter
        des = self._disorder_elements()
        orbits = des.keys()
        gis = self._group_indices

        # TODO: This exact pattern has been copy+pasted twice...
        pi = iter(template_pattern)
        for orbit in orbits:
            indices = gis[orbit]
            de = des[orbit]
            for e in de:
                subpattern = next(pi)
                for i in subpattern:
                    template_structure.sites[indices[i]].species = e

        cifwriter = CifWriter(template_structure)
        cfkey = cifwriter.ciffile.data.keys()
        cfkey = list(cfkey)[0]

        # Use faster CifBlock implementation
        block = AltCifBlock.from_string(str(cifwriter.ciffile.data[cfkey]))
        cifwriter.ciffile.data[cfkey] = block

        if symprec is None:
            template_type_symbol = block["_atom_site_type_symbol"]
            template_label = block["_atom_site_label"]

            # These two blocks are always "1" in the final structure
            block["_atom_site_symmetry_multiplicity"] = ["1"] * len(template_label)
            block["_atom_site_occupancy"] = ["1.0"] * len(template_label)

            for i in self.sampled_indices:
                type_symbol = template_type_symbol.copy()
                label = template_label.copy()

                pi = iter(self._patterns[i])
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
                yield cifwriter
        else:
            format_str = "{:.%df}" % 8
            latt = template_structure.lattice.matrix
            positions = template_structure.frac_coords

            cell_specie = list(set(x.species for x in template_structure))
            # Flattened list of species @ disorder sites
            specie = [y for x in des.values() for y in x]
            z_map = [
                cell_specie.index(Composition({specie[j]: 1}))
                for j in range(len(template_pattern))
            ]
            template_zs = [cell_specie.index(x.species) for x in template_structure]

            for i in self.sampled_indices:
                zs = template_zs.copy()

                pi = iter(self._patterns[i])
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
                        space_group_data["rotations"], space_group_data["translations"]
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
                                abs(x) for x in template_structure.sites[s].frac_coords
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
                    block["_atom_site_label"].append("{}{}".format(sp.symbol, count))
                    block["_atom_site_symmetry_multiplicity"].append(str(mult))
                    block["_atom_site_fract_x"].append(format_str.format(site.a))
                    block["_atom_site_fract_y"].append(format_str.format(site.b))
                    block["_atom_site_fract_z"].append(format_str.format(site.c))
                    block["_atom_site_occupancy"].append("1.0")
                    count += 1
                yield cifwriter


class PatternMaker:
    """
    Generate permutationally-unique patterns
    """

    __slots__ = (
        "search",
        "_enumerator_collection",
        "_patterns",
        "_auts",
        "_subobj_ts",
        "_bs",
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
        "_sieve",
        "_get_not_reject_mask",
        "_get_accept_mask",
        "_get_mins",
        "_gen_flag",
    )

    def __init__(
        self,
        perm_list,
        invar=None,
        enumerator_collection=None,
        t_kind=const.DEFAULT_T_KIND,
    ):
        # Algo select bits
        if invar is None:
            self.search = self._invarless_search
        else:
            self.search = self._invar_search

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
        self._bits = [2 ** i for i in range(self._nix)]
        # Convert to list to avoid overflow
        # Will change type to object for very large integer
        # Careful: each are still numpy object type...
        try:
            self._bit_perm = np.array(
                [[2 ** int(i) for i in row] for row in indexed_perm_list], dtype=int
            )
        except OverflowError:
            self._bit_perm = np.array(
                [[2 ** int(i) for i in row] for row in indexed_perm_list], dtype=object
            )

        self._gen_flag = collections.defaultdict(bool)

        # Containers for search results.
        self._patterns = collections.defaultdict(list)
        self._patterns[0] = [{}]
        self._auts = collections.defaultdict(list)
        self._auts[0] = np.ones((self._nperm,), dtype="bool")
        self._subobj_ts = collections.defaultdict(list)
        self._subobj_ts[0] = [np.array([0])]
        self._bs = collections.defaultdict(list) 
        self._bs[0] = [np.zeros(self._nix)]

        self.label = self._perms.tobytes()

    @staticmethod
    def reindex(perm_list):
        """Standardize symmetry ordering for reuse (rough)"""
        perm_list = perm_list.copy()

        # Column sort
        stab_map = perm_list == perm_list[0]
        column_index = np.lexsort(stab_map)
        perm_list = perm_list[:, column_index]

        # Relabel to match column position
        relabel_index = perm_list[0]
        relabel_element = np.vectorize({s: i for i, s in enumerate(relabel_index)}.get)
        perm_list = relabel_element(perm_list)

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

        # Clear caches.
        self.auts.cache_clear()
        self.patterns.cache_clear()

    @staticmethod
    def get_label(perm_list):
        """
        Rough equivalence. Byte string from a standardized
        representation of the given permutation
        """
        _, _, _, relabeled_perm_list = PatternMaker.reindex(perm_list)
        return relabeled_perm_list.tobytes()

    @functools.lru_cache(None)
    def patterns(self, n):
        """
        Get patterns for the specified replacement amount
        """
        # Patterns are symmetrical
        _n = min(n, self._nix - n)
        if not self._gen_flag[_n]:
            lessthan = [i for i in self._gen_flag.keys() if i < _n]
            if not lessthan:
                start = 0
            else:
                start = max(lessthan)
            self.search(start=start, stop=_n)
        # "Mirror" patterns
        if _n != n:
            inverter = np.arange(self._nix)
            return [
                self._relabel_index[np.setdiff1d(inverter, pattern)]
                for pattern in self._patterns[_n]
            ]
        return [self._relabel_index[pattern] for pattern in self._patterns[_n]]

    @functools.lru_cache(None)
    def auts(self, n):
        """
        Get the automorphisms in terms of input permutation
        """
        # Patterns are symmetrical
        _n = min(n, self._nix - n)
        if not self._gen_flag[_n]:
            lessthan = [i for i in self._gen_flag.keys() if i < _n]
            if not lessthan:
                start = 0
            else:
                start = max(lessthan)
            self.search(start=start, stop=_n)
        # Map to original permutation; put identity on 0
        auts = [np.sort(self._row_index[aut]) for aut in self._auts[_n]]
        return auts

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
        lessthan = [i for i in self._gen_flag.keys() if i < stop]
        if not lessthan:
            start = 0
        else:
            start = max(lessthan, start)

        # Cross-check with exact enumeration.
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
            patterns = self._patterns[start].copy()
            auts = self._auts[start].copy()
            bs = self._bs[start].copy()
            for pattern, aut in zip(patterns, auts):
                stack.append((pattern, aut, bs))

        while stack:
            pattern, aut, pbs = stack.pop()
            if pattern.size == stop:
                self._patterns[stop].append(pattern)
                self._auts[stop].append(pattern)
                self._bs[pattern.size].append(pbs)
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

            for i in np.flatnonzero(uniq_mask):
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

        n_gen = len(self._patterns[stop])
        if n_pred != n_gen:
            raise RuntimeError(
                "Mismatch between predicted and generated number of structures.\n"
                f"(at {stop}, {n_gen}/{n_pred} were generated)"
            )
        tqdm.tqdm(disable=const.DISABLE_PROGRESSBAR).write("Done.")
        self._gen_flag[stop] = True

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

    # @profile
    def _invar_search(self, start=0, stop=None):
        """
        Combines invar and lexmax to determine canonical parent.

        Args:
            start: Start seach at this depth.
            stop: Stop search at this depth.
        """
        # TODO: stop and there
        if stop is None:
            stop = self._nix // 2
        stop = min(stop, self._nix - stop)
        lessthan = [i for i in self._gen_flag.keys() if i < stop]
        if not lessthan:
            start = 0
        else:
            start = max(lessthan, start)

        # Cross-check with exact enumeration.
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
            patterns = self._patterns[start].copy()
            auts = self._auts[start].copy()
            sums = self._subobj_ts[start].copy()
            # TODO: create from scratch, instead of filling the memory.
            bs = self._bs[start].copy()
            for subobj_ts, pattern, aut in zip(sums, patterns, auts):
                stack.append((subobj_ts, pattern, aut, bs))

        while stack:
            subobj_ts, pattern, aut, pbs = stack.pop()
            if pattern.size == stop:
                self._patterns[pattern.size].append(pattern)
                self._auts[pattern.size].append(aut)
                self._subobj_ts[pattern.size].append(subobj_ts)
                self._bs[pattern.size].append(pbs)
                pbar.update()
                continue
            # Tree is expanded by adding un-added sites
            leaf_mask = np.ones(self._nix, dtype="bool")
            leaf_mask[pattern] = False
            leaf_array = np.flatnonzero(leaf_mask)

            # Calculate subobject Ts for all leaves
            leaf_subobj_ts = self._get_subobj_ts(pattern, leaf_array, subobj_ts)

            # Reject all leaves where any T is smaller than the new row's T
            delta_t = leaf_subobj_ts[:, :-1] - leaf_subobj_ts[:, -1:]
            not_reject_mask = self._get_not_reject_mask(delta_t)
            if not not_reject_mask.any():
                continue

            # Discard symmetry duplicates from the remaining leaves
            if aut.size > 1 and not_reject_mask.sum() > 1:
                not_reject_leaf = leaf_array[not_reject_mask]
                leaf_reps = self._perms[np.ix_(aut, not_reject_leaf)].min(axis=0)
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
            for i in np.flatnonzero(uniq_mask):
                x = leaf_array[i]
                j = loci[i]
                _pbs = self._bit_perm[:, x] + pbs

                _subobj_ts = leaf_subobj_ts[i]
                _subobj_ts[j:] = np.concatenate((_subobj_ts[-1:], _subobj_ts[j:-1]))

                _i = np.concatenate((pattern[:j], [x], pattern[j:]))
                # NOTE: just in case I fail to consistently sort perm
                # bitsum = sum([self._bits[y] for y in _pattern])
                # _aut = np.flatnonzero(_pbs[t] == bitsum)
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

        n_gen = len(self._patterns[stop])
        if n_pred != n_gen:
            raise RuntimeError(
                "Mismatch between predicted and generated number of structures.\n"
                f"(at {stop}, {n_gen}/{n_pred} were generated)"
            )
        tqdm.tqdm(disable=const.DISABLE_PROGRESSBAR).write("Done.")
        self._gen_flag[stop] = True


class Polya:
    """
    Perform operations related (weighed) Polya Enumeration Theorem.

    Includes enumeration (count) and calcuation of cycle index
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
        return sum(self.ci().values()) / self.group_size

    @functools.lru_cache(None)
    def ci(self):
        """
        Returns the cycle index.

        The terms are separated by each permutation (rows)
        to ease further processing.
        Uses sympy.IndexedBase to represent subscript.

        Conventionally, each term is in the form of (x_n**m),
        with (n) denoting cycle length and (m) number of
        cycles with such length for a permutation.

        Returns:
            dict: Dictionary of cycle index with integer index
                and sympy equations as a key-value pair.
        """
        # Use indexed base.
        cycle_index = dict()
        for j, permutations in enumerate(self._perm_list):
            # TODO: What is the correct behaviour, if one permutation is empty?
            symbol = sympy.IndexedBase(chr(97 + j))
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
                            raise RuntimeError("Check permutation list.") from exc
                    cycles.append(cycle)

                # Cast to (sympy.Expr)
                counter = collections.Counter(len(cycle) for cycle in cycles)
                cycle_index.setdefault(i, sympy.Integer(1))
                for length, n in counter.items():
                    cycle_index[i] *= symbol[length] ** n
        return cycle_index

    @functools.lru_cache(None)
    def cgf(self, len_tuple):
        """
        Returns the configuration generation function.

        Basically, each x_n**m terms of the cycle index
        will be substituted by (a**m + b**m +...), with
        the number of terms corresponds to the final
        number of color / species.

        Terms coming from different orbit will be given
        a distinct letter, but grouped together within
        the same dictionary key.

        Args:
            len_list (tuple): Number of distinct species in the final
            structure, for each separated orbit.

        Returns:
            tuple (dict, list): Combined results of:
                dict: CGF in similar format to the cycle index.
                list: List of letters used in the CGF.
                    Useful for example in the counting where
                    we want to find a specific term.
        """
        len_tuple_string = ", ".join(map(str, len_tuple))
        logging.info(f"Building CGF for {len_tuple_string}.")

        def divide(iterable, sizes):
            i = 0
            for j in sizes:
                yield iterable[i : i + j]
                i += j

        # Sanity check.
        if len(len_tuple) != len(self._perm_list):
            raise ValueError("Invalid sequences length.")

        # Avoid variable name collisions. Start from "a".
        ci_letters = [chr(97 + i) for i in range(len(self._perm_list))]
        sub_letters = [chr(ord(ci_letters[-1]) + 1 + i) for i in range(sum(len_tuple))]
        sub_letters = list(divide(sub_letters, len_tuple))
        replacement_map = {
            sympy.Symbol(l): sum(sympy.Symbol(e) for e in r)
            for l, r in zip(ci_letters, sub_letters)
        }

        cgf = dict()
        for i, expr in tqdm.tqdm(
            self.ci().items(),
            desc="Substituting terms",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        ):
            indexed_symbols = [x for x in expr.free_symbols if isinstance(x, Indexed)]

            replacements = []
            for indexed_symbol in indexed_symbols:
                label = indexed_symbol.base.label
                index = indexed_symbol.indices[0]
                replacement = replacement_map[label]

                raised = []
                for symbol in replacement.free_symbols:
                    raised.append(symbol ** index)

                raised_rep = replacement.xreplace(
                    dict(zip(replacement.free_symbols, raised))
                )
                replacements.append(raised_rep)
            expr = expr.xreplace(dict(zip(indexed_symbols, replacements)))

            # TODO: consider symbol symmetry? and shared cache
            if expr in self._expand_cache:
                expanded = self._expand_cache[expr]
            else:
                expanded = expr.expand() / self.group_size
                self._expand_cache[expr] = expanded
            cgf[i] = expanded
        return cgf, sub_letters

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
        logging.info(f"Counting unique patterns for {amt_tuple}.")
        _sequences = [len(x) for x in amt_tuple]
        cgf, sub_letters = self.cgf(tuple(_sequences))

        term = sympy.Integer(1)
        for letters, nums in zip(sub_letters, amt_tuple):
            for letter, num in zip(letters, nums):
                term *= sympy.Symbol(letter) ** num

        count = int(sum(c.coeff(term) for c in cgf.values()))
        logging.info(count)
        return count


# TODO: shared expand cache
class PolyaCollection:
    """
    Collection of Polya objects. This is useful because identical
    instances are often recreated, with heavy sympy.expand() operations.
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
