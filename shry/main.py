# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Main task abstraction
"""

# python modules
import ast
import collections
import configparser
import datetime
import io
import itertools
import logging
import math
import operator
import os
import re
import shutil
import signal
import sys
from fnmatch import fnmatch

# python modules
import numpy as np
import tqdm
from pymatgen.core import Composition, PeriodicSite, Species, Structure
from pymatgen.core.composition import CompositionError
from pymatgen.core.lattice import Lattice
from pymatgen.io.cif import CifParser, str2float
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.util.coord import in_coord_list_pbc, lattice_points_in_supercell

# shry modules
from . import const
from .core import Substitutor

# shry version control
try:
    from ._version import version as shry_version
except (ModuleNotFoundError, ImportError):
    shry_version = "unknown (not installed)"

# Runtime patches to pymatgen's Composition and Poscar.
# We want some extra functions like separating different oxidation states, etc.


def from_string(composition_string) -> Composition:
    """
    Workaround when working with strings including oxication states
    """

    def composition_builder():
        components = re.findall(
            r"[A-z][a-z]*[0-9.\+\-]*[0-9.]*", composition_string
        )
        amount = re.compile(r"(?<=[\+\-A-Za-z])[0-9.]+(?![\+\-])")
        amount_part = [amount.findall(x) for x in components]
        amount_part = [x[0] if x else "1" for x in amount_part]
        species_part = [
            x.strip(amount) for x, amount in zip(components, amount_part)
        ]
        species_part = [
            x + "0" if not re.search(r"[0-9\+\-]+", x) else x
            for x in species_part
        ]
        amount_part = [float(x) for x in amount_part]

        return Composition(
            {
                Species.from_string(species): amount
                for species, amount in zip(species_part, amount_part)
            }
        )

    try:
        return Composition(composition_string)
    except (CompositionError, ValueError, IndexError):
        # - CompositionError: get_sym_dict() error
        #   when "+" oxidation is used.
        # - ValueError: Wrong parse by get_sym_dict()
        #   if the oxidation negative.
        # - IndexError: Sometimes appear when the string gets more complex
        return composition_builder()


@property
def site_symbols(self):
    """
    On Poscar: sometimes we would like to use a separate pseudopotential
    for each oxidation states.

    Write oxidation states, if, and only if, for the given element,
    there are more than 1 state.
    """
    return [a[0] for a in itertools.groupby(self.syms)]


@property
def syms(self):
    """
    Replaced self.syms in some functions of Poscar
    """
    e_oxstates = collections.defaultdict(set)
    s_symbols = collections.defaultdict(list)

    for site in self.structure:
        e_oxstates[site.specie.symbol].add(str(site.specie))
    for e, oes in e_oxstates.items():
        if len(oes) == 1:
            s_symbols[list(oes)[0]] = e
        else:
            for oe in oes:
                s_symbols[oe] = oe

    return [s_symbols[str(site.specie)] for site in self.structure]


@property
def natoms(self):
    """
    Sequence of number of sites of each type associated with the Poscar.
    Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
    """
    return [len(tuple(a[1])) for a in itertools.groupby(self.syms)]


@staticmethod
def parse_oxi_states(data):
    """
    Parse oxidation states from data dictionary
    """
    ox_state_regex = re.compile(r"\d?[\+,\-]?$")

    try:
        oxi_states = {
            data["_atom_type_symbol"][i]: str2float(
                data["_atom_type_oxidation_number"][i]
            )
            for i in range(len(data["_atom_type_symbol"]))
        }
        # attempt to strip oxidation state from _atom_type_symbol
        # in case the label does not contain an oxidation state
        for i, symbol in enumerate(data["_atom_type_symbol"]):
            oxi_states[ox_state_regex.sub("", symbol)] = str2float(
                data["_atom_type_oxidation_number"][i]
            )
    except (ValueError, KeyError):
        # Some CIF (including pymatgen's output) are like this.
        try:
            oxi_states = dict()
            for i, symbol in enumerate(data["_atom_site_type_symbol"]):
                _symbol = ox_state_regex.sub("", symbol)

                parsed_oxi_state = ox_state_regex.search(symbol).group(0)
                if not parsed_oxi_state:
                    oxi_states[_symbol] = None
                    continue

                sign = re.search(r"[-+]", parsed_oxi_state)
                if sign is None:
                    sign = ""
                else:
                    sign = sign.group(0)
                parsed_oxi_state = parsed_oxi_state.replace("+", "").replace(
                    "-", ""
                )
                parsed_oxi_state = str2float(sign + parsed_oxi_state)
                oxi_states[ox_state_regex.sub("", symbol)] = parsed_oxi_state
        except (ValueError, KeyError):
            oxi_states = None
    return oxi_states


Composition.from_string = from_string
Poscar.site_symbols = site_symbols
Poscar.natoms = natoms
Poscar.syms = syms
CifParser.parse_oxi_states = parse_oxi_states


class ScriptHelper:
    """
    Combine configurations into typical workflow in single methods
    """

    def __init__(
        self,
        structure_file,
        from_species=const.DEFAULT_FROM_SPECIES,
        to_species=const.DEFAULT_TO_SPECIES,
        scaling_matrix=const.DEFAULT_SCALING_MATRIX,
        symmetrize=const.DEFAULT_SYMMETRIZE,
        sample=const.DEFAULT_SAMPLE,
        seed=const.DEFAULT_SEED,
        symprec=const.DEFAULT_SYMPREC,
        atol=const.DEFAULT_ATOL,
        angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        dir_size=const.DEFAULT_DIR_SIZE,
        write_symm=const.DEFAULT_WRITE_SYMM,
        write_ewald=const.DEFAULT_WRITE_EWALD,
        max_ewald=const.DEFAULT_MAX_EWALD,
        no_write=const.DEFAULT_NO_WRITE,
        no_dmat=const.DEFAULT_NO_DMAT,
        no_cache=False,
        t_kind=const.DEFAULT_T_KIND,
    ):
        # TODO: Refactor with descriptors.
        self._timestamp = datetime.datetime.now().timestamp()
        self.no_write = no_write
        self.no_dmat = no_dmat
        self.t_kind = t_kind
        self.structure_file = structure_file
        self._seed = seed

        if len(from_species) != len(to_species):
            raise RuntimeError(
                "from_species and to_species must have the same length."
            )
        self.from_species = from_species
        self.to_species = to_species
        self.scaling_matrix = np.array(scaling_matrix)
        self.symmetrize = symmetrize

        if isinstance(sample, str):
            if sample == "all":
                sample = None
            else:
                sample = int(self._math_eval(sample))
        self.sample = sample
        self.symprec = symprec
        self.atol = atol
        self.angle_tolerance = angle_tolerance
        self.dir_size = dir_size
        self.write_symm = write_symm
        self.write_ewald = write_ewald
        self.max_ewald = max_ewald
        if self.max_ewald is not None:
            self.write_ewald = True

        logging.info("\nRun configurations:")
        logging.info(const.HLINE)
        logging.info(self)
        logging.info(const.HLINE)

        # TODO: these "derivative attributes" should have a proper setter()
        # In fact any that are not simple assignment should.
        if no_cache:
            cache = False
        else:
            cache = None
        self.structure = LabeledStructure.from_file(
            self.structure_file,
            symmetrize=symmetrize,
        )
        self.modified_structure = self.structure.copy()
        # Note: since we don't limit the allowable scaling_matrix,
        # changing the order of enlargement vs. replace can have different meaning.
        self.modified_structure.replace_species(
            dict(zip(self.from_species, self.to_species))
        )
        self.modified_structure *= self.scaling_matrix
        self.substitutor = Substitutor(
            self.modified_structure,
            symprec=self.symprec,
            atol=self.atol,
            angle_tolerance=self.angle_tolerance,
            shuffle=self.sample is not None,
            seed=seed,
            no_dmat=self.no_dmat,
            t_kind=self.t_kind,
            cache=cache,
        )

    def __str__(self):
        string = ""
        print_format = "  * {} = {}\n"
        string += print_format.format("SHRY version", shry_version)
        string += print_format.format("structure_file", self.structure_file)
        string += print_format.format(
            "from_species", ", ".join(map(str, self.from_species))
        )
        string += print_format.format(
            "to_species", ", ".join(map(str, self.to_species))
        )
        string += print_format.format(
            "scaling_matrix", np.array(self.scaling_matrix).flatten()
        )
        string += print_format.format("symmetrize", self.symmetrize)
        string += print_format.format("sample", self.sample)
        string += print_format.format("symprec", self.symprec)
        string += print_format.format("angle_tolerance", self.angle_tolerance)
        string += print_format.format("dir_size", self.dir_size)
        string += print_format.format("write_symm", self.write_symm)
        string += print_format.format("t-kind", self.t_kind)
        return string

    @staticmethod
    def _math_eval(expr):
        """
        Like eval() but limit to +-/*^** expressions for security
        """
        operators = {
            ast.Add: operator.add,
            ast.Sub: operator.sub,
            ast.Mult: operator.mul,
            ast.Div: operator.truediv,
            ast.Pow: operator.pow,
            ast.BitXor: operator.xor,
            ast.USub: operator.neg,
        }

        def _tracer(node):
            if isinstance(node, ast.Num):
                return node.n
            elif isinstance(node, ast.BinOp):
                return operators[type(node.op)](
                    _tracer(node.left), _tracer(node.right)
                )
            elif isinstance(node, ast.UnaryOp):
                return operators[type(node.op)](_tracer(node.operand))
            else:
                raise TypeError("Invalid characters on math expression")

        return _tracer(ast.parse(expr, mode="eval").body)

    @classmethod
    def from_file(cls, input_file):
        """
        Reads *.ini file containing command line arguments
        """
        # TODO: simplify
        parser = configparser.ConfigParser(
            empty_lines_in_values=False, allow_no_value=False
        )
        encoding = getattr(io, "LOCALE_ENCODING", "utf8")
        with open(
            input_file, "r", encoding=encoding, errors="surrogateescape"
        ) as f:
            parser.read_file(f)

        structure_file = parser.get("DEFAULT", "structure_file")
        from_species = parser.get(
            "DEFAULT", "from_species", fallback=const.DEFAULT_FROM_SPECIES
        )
        to_species = parser.get(
            "DEFAULT", "to_species", fallback=const.DEFAULT_TO_SPECIES
        )
        scaling_matrix = parser.get(
            "DEFAULT",
            "scaling_matrix",
            fallback=const.DEFAULT_SCALING_MATRIX_STR,
        )
        # Allow ",", ";", and whiteline as separator (';' not allowed by *.ini)
        from_species = const.FLEXIBLE_SEPARATOR.split(from_species)
        to_species = const.FLEXIBLE_SEPARATOR.split(to_species)
        scaling_matrix = [
            int(x) for x in const.FLEXIBLE_SEPARATOR.split(scaling_matrix)
        ]
        symmetrize = parser.getboolean(
            "DEFAULT", "symmetrize", fallback=const.DEFAULT_SYMMETRIZE
        )
        sample = parser.getint(
            "DEFAULT", "sample", fallback=const.DEFAULT_SAMPLE
        )
        symprec = parser.getfloat(
            "DEFAULT", "symprec", fallback=const.DEFAULT_SYMPREC
        )
        angle_tolerance = parser.getfloat(
            "DEFAULT",
            "angle_tolerance",
            fallback=const.DEFAULT_ANGLE_TOLERANCE,
        )
        dir_size = parser.getint(
            "DEFAULT", "dir_size", fallback=const.DEFAULT_DIR_SIZE
        )
        write_symm = parser.getboolean(
            "DEFAULT", "write_symm", fallback=const.DEFAULT_WRITE_SYMM
        )

        no_write = parser.getboolean(
            "DEFAULT", "no_write", fallback=const.DEFAULT_NO_WRITE
        )
        no_dmat = parser.getboolean(
            "DEFAULT", "no_dmat", fallback=const.DEFAULT_NO_DMAT
        )
        t_kind = parser.getboolean(
            "DEFAULT", "t_kind", fallback=const.DEFAULT_T_KIND
        )

        return cls(
            structure_file=structure_file,
            from_species=from_species,
            to_species=to_species,
            scaling_matrix=scaling_matrix,
            symmetrize=symmetrize,
            sample=sample,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            dir_size=dir_size,
            write_symm=write_symm,
            no_write=no_write,
            no_dmat=no_dmat,
            t_kind=t_kind,
        )

    @property
    def _outdir(self):
        """
        Output directory name based on input structure filename and timestamp.

        Returns:
            str: Output directory name.
        """
        structure_file_basename = os.path.basename(self.structure_file).split(
            "."
        )[0]
        return f"shry-{structure_file_basename}-{self._timestamp}"

    def save_modified_structure(self):
        """
        Save adjusted structure
        """
        if self.no_write:
            return
        os.makedirs(self._outdir, exist_ok=True)
        structure_file_basename = os.path.basename(self.structure_file).split(
            "."
        )[0]
        scaling_matrix_string = "-".join(
            self.scaling_matrix.flatten().astype(str)
        )
        filename = os.path.join(
            self._outdir,
            structure_file_basename + "-" + scaling_matrix_string + ".cif",
        )
        self.modified_structure.to(filename=filename, symprec=self.symprec)

    def write(self):
        """
        Save the irreducible structures
        """
        npatterns = self.substitutor.count()
        if not npatterns:
            logging.warning("No expected patterns.")
            return
        if self.sample is not None:
            npatterns = self.sample

        if self.no_write:
            # (for io-less benchmark)
            pbar = tqdm.tqdm(
                total=npatterns,
                desc=f"Generating {npatterns} order structures",
                **const.TQDM_CONF,
                disable=const.DISABLE_PROGRESSBAR,
            )
            for _ in self.substitutor.make_patterns():
                pbar.update()
            pbar.close()
            return

        # Log file stream
        logio = io.StringIO()

        def dump_log():
            """
            Dump log
            """
            logfile = os.path.join(self._outdir, "sub.log")
            encoding = getattr(io, "LOCALE_ENCODING", "utf8")
            with open(
                logfile, "w", encoding=encoding, errors="surrogateescape"
            ) as f:
                logio.seek(0)
                shutil.copyfileobj(logio, f)

        def signal_handler(sig, frame):
            """
            Write log when interrupted
            """
            del sig, frame
            dump_log()
            now = datetime.datetime.now()
            tz = now.astimezone().tzname()
            time_string = now.strftime("%c ") + tz
            logging.info(const.HLINE)
            logging.info("Aborted %s", time_string)
            logging.info(const.HLINE)
            sys.exit(0)

        signal.signal(signal.SIGINT, signal_handler)

        # Filenames
        # Formatting.
        def outdir(i):
            return os.path.join(self._outdir, f"slice{i // self.dir_size}")

        formula = "".join(self.modified_structure.formula.split())
        ndigits = int(math.log10(npatterns)) + 1
        index_f = "_{:0" + str(ndigits) + "d}"
        filenames = [
            os.path.join(outdir(i), formula + index_f.format(i))
            for i in range(npatterns)
        ]

        # Make directories
        ndirs = npatterns // self.dir_size + 1
        for i in range(ndirs):
            os.makedirs(os.path.join(self._outdir, f"slice{i}"), exist_ok=True)

        # Save the structures
        pbar = tqdm.tqdm(
            total=npatterns,
            desc=f"Writing {npatterns} order structures",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        )
        os.makedirs(os.path.join(self._outdir, "slice0"), exist_ok=True)

        # TODO: Refactor
        header = "N Weight Configuration"
        if self.write_ewald:
            header = header + " EwaldEnergy"
            quantities = ("cifwriter", "weight", "letter", "ewald")
        else:
            quantities = ("cifwriter", "weight", "letter")

        if self.write_symm:
            header = header + " GroupName"
            symprec = self.symprec
        else:
            symprec = None

        if self.sample is not None:
            header += f" (seed={self._seed})"

        quantities_generator = enumerate(
            self.substitutor.quantities(quantities, symprec)
        )
        # Limit amount if sampled
        if self.sample is not None:
            quantities_generator = itertools.takewhile(
                lambda x: x[0] < npatterns,
                quantities_generator,
            )

        print(header, file=logio)
        for i, packet in quantities_generator:
            cifwriter = packet["cifwriter"]
            ewald = packet["ewald"]
            weight = packet["weight"]
            letter = packet["letter"]

            if self.max_ewald is not None and ewald > self.max_ewald:
                continue

            line = f"{i} {weight} {letter}"
            if ewald is not None:
                line = line + f" {ewald}"
            if self.write_symm:
                space_group = list(cifwriter.ciffile.data.values())[0][
                    "_symmetry_space_group_name_H-M"
                ]
                line += line + f" {space_group}"
            print(line, file=logio)

            try:
                cifwriter.write_file(filename=filenames[i] + f"_{weight}.cif")
                pbar.update()
            # Warn if too many structures are generated
            except IndexError:
                logging.warning(
                    f"Too many structures generated (>{npatterns}). "
                    "Check `atol` value."
                )
                break
        pbar.close()
        dump_log()

        # Warn if too little structures are generated
        if i < npatterns - 1:
            logging.warning(
                f"Too little structures generated ({i+1}/{npatterns}). "
                "Check `atol` value."
            )

    def count(self) -> None:
        """
        Count the number of unique substituted structures
        """
        count = self.substitutor.count()
        total_count = (
            self.substitutor.total_count()
        )  # pylint: disable=assignment-from-no-return
        logging.info(const.HLINE)
        logging.info(f"Total number of combinations is {total_count}")
        logging.info(f"Expected unique patterns is {count}")
        logging.info(const.HLINE)
        return count


class LabeledStructure(Structure):
    """
    Structure + CIF's _atom_site_label
    """

    def __str__(self):
        return "LabeledStructure\n" + super().__str__()

    def __mul__(self, scaling_matrix):
        """
        The parent method returns Structure instance!
        Overwrite the offending line.
        """

        scale_matrix = np.array(scaling_matrix, np.int16)

        # check the shape of the scaling matrix
        if scale_matrix.shape not in {(1,), (3,), (3, 3)}:
            logging.warning(
                "The scale_matrix.shape should be (1,), (3,), or (3, 3)"
            )
            raise ValueError
        else:
            if scale_matrix.shape != (3, 3):
                scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
            else:
                pass

        new_lattice = Lattice(np.dot(scale_matrix, self._lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        for site in self:
            for v in c_lat:
                s = PeriodicSite(
                    site.species,
                    site.coords + v,
                    new_lattice,
                    properties=site.properties,
                    coords_are_cartesian=True,
                    to_unit_cell=False,
                    skip_checks=True,
                )
                new_sites.append(s)

        new_charge = (
            self._charge * np.linalg.det(scale_matrix)
            if self._charge
            else None
        )
        # This line.
        return self.__class__.from_sites(new_sites, charge=new_charge)

    @classmethod
    def from_file(  # pylint: disable=arguments-differ
        cls,
        filename,
        primitive=False,
        sort=False,
        merge_tol=0.0,
        symmetrize=False,
    ):
        fname = os.path.basename(filename)
        if not fnmatch(fname.lower(), "*.cif*") and not fnmatch(
            fname.lower(), "*.mcif*"
        ):
            raise ValueError("LabeledStructure only accepts CIFs.")
        instance = super().from_file(
            filename, primitive=primitive, sort=sort, merge_tol=merge_tol
        )

        instance.read_label(filename, symmetrize=symmetrize)
        return instance

    def replace_species(self, species_mapping):
        """
        Replace sites from species to another species,
        or labels from _atom_site_label.

        Args:
            species_mapping (dict): from-to map of the species to be replaced.
                (string) to (string).

        Raises:
            RuntimeError: If no sites matches at least one of the from_species
                specification.
        """
        for map_num, mapping in enumerate(species_mapping.items()):
            from_species, to_species = mapping
            replace = False
            to_composition = Composition.from_string(to_species)
            for site in self:
                # Find matching site label
                if from_species in site.properties["_atom_site_label"]:
                    replace = True
                    # If match by label, don't change the label!
                    # site.properties["_atom_site_label"] = tuple(sorted({to_species}))
                    try:
                        site.species = to_composition
                    except ValueError:
                        site.species = to_composition.fractional_composition
                # If failed, try to find matching Element
                elif any(
                    e.symbol == from_species for e in site.species.elements
                ):
                    replace = True
                    # Since all sites are replaced, merge them under one label.
                    # But give a distinct ID in case two or more sites are
                    # replaced into the same species.
                    new_label = tuple(sorted({to_species}) + [map_num])
                    site.properties["_atom_site_label"] = new_label
                    try:
                        site.species = to_composition
                    except ValueError:
                        site.species = to_composition.fractional_composition

            if not replace:
                raise RuntimeError(
                    f"Can't find the specified site ({from_species})."
                )

    def read_label(
        self,
        cif_filename,
        symprec=const.DEFAULT_SYMPREC,
        angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        symmetrize=False,
    ):
        """
        Add _atom_site_label as site_properties.
        This is useful for enforcing a certain concentration over
        a group of sites, which may not necessarily consists
        of a single orbit after enlargement to a supercell.

        Args:
            cif_filename (str): Source CIF file.
            symprec (float): Precision for the symmetry operations.

        Raises:
            RuntimeError: If any sites couldn't be matched
                to one any sites defined within the CIF.
        """
        logging.info(f"Reading _atom_site_label from {cif_filename}")

        def ufloat(string):
            """Remove uncertainties notion from floats (if any)."""
            try:
                return float(string)
            except ValueError:
                return float(string.split("(")[0])

        encoding = getattr(io, "LOCALE_ENCODING", "utf8")
        with open(
            cif_filename, "r", encoding=encoding, errors="surrogateescape"
        ) as f:
            parser = CifParser.from_string(f.read())

        # Since Structure only takes the first structure inside a CIF, do the same.
        cif_dict = list(parser.as_dict().values())[0]
        labels = cif_dict["_atom_site_label"]
        x_list = map(ufloat, cif_dict["_atom_site_fract_x"])
        y_list = map(ufloat, cif_dict["_atom_site_fract_y"])
        z_list = map(ufloat, cif_dict["_atom_site_fract_z"])
        coords = [(x, y, z) for x, y, z in zip(x_list, y_list, z_list)]

        # Merge labels to allow multiple references.
        cif_sites = []
        for coord, zipgroup in itertools.groupby(
            zip(coords, labels), key=lambda x: x[0]
        ):
            labels = tuple(sorted({x[1] for x in zipgroup}))
            cif_sites.append(
                PeriodicSite(
                    "X",
                    coord,
                    self.lattice,
                    properties={"_atom_site_label": labels},
                )
            )

        # Find equivalent sites.
        if symmetrize:
            symm_ops = SpacegroupAnalyzer(
                self, symprec, angle_tolerance
            ).get_space_group_operations()
        else:
            # A bit of trick.
            parser.data = cif_dict
            # Spacegroup symbol and number are not important here.
            symm_ops = SpacegroupOperations(0, 0, parser.get_symops(parser))

        coords = [x.frac_coords for x in self.sites]
        cif_coords = [x.frac_coords for x in cif_sites]

        # List of coordinates that are equivalent to this site
        o_cif_coords = [
            symmop.operate_multi(cif_coords) for symmop in symm_ops
        ]
        o_cif_coords = np.swapaxes(np.stack(o_cif_coords), 0, 1)
        o_cif_coords = [np.unique(np.mod(x, 1), axis=0) for x in o_cif_coords]

        for site in tqdm.tqdm(
            self.sites,
            desc="Matching CIF labels",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        ):
            equivalent = [
                in_coord_list_pbc(o, site.frac_coords) for o in o_cif_coords
            ]

            try:
                equivalent_site = cif_sites[equivalent.index(True)]
                site.properties[
                    "_atom_site_label"
                ] = equivalent_site.properties["_atom_site_label"]
            except ValueError as exc:
                raise RuntimeError("CIF-Structure mismatch.") from exc
