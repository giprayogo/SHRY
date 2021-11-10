#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=logging-fstring-interpolation, logging-not-lazy

# information
__author__ = "Genki Prayogo, and Kosuke Nakano"
__copyright__ = "Copyright (c) 2021-, The SHRY Project"
__credits__ = ["Genki Prayogo", "Kosuke Nakano"]

__license__ = "MIT"
__version__ = "1.0.0"
__maintainer__ = "Genki Prayogo"
__email__ = "g.prayogo@icloud.com"
__date__ = "2. Nov. 2021"
__status__ = "Production"

"""Command line interface."""
import argparse
import datetime
import fnmatch
import logging
# import sys

import tqdm

from . import const

# Disable detailed stack trace.
# sys.excepthook = lambda t, e, _: print(f"{t.__name__}: {e}")


class TqdmLoggingHandler(logging.Handler):
    """
    Prevent logging from overwriting tqdm
    """

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:  # pylint: disable=bare-except
            self.handleError(record)


def print_header():
    """
    Print header text
    """
    now = datetime.datetime.now()
    tz = now.astimezone().tzname()
    time_string = now.strftime("%c ") + tz
    bold = "\033[1m"
    end = "\033[0m"

    # logger.setLevel(logging.INFO)
    handler = TqdmLoggingHandler()
    handler.setLevel(logging.INFO)
    logging.basicConfig(format="%(message)s", level=logging.INFO, handlers=[handler])
    logging.info("********************************\n")
    logging.info(
        f"SHRY: Suite for High-throughput generation of models"
        f"with atomic substitutions implemented by Python"
    )
    logging.info("\n********************************")
    logging.info("Begin " + time_string)


def print_footer():
    """
    Print footer text
    """
    now = datetime.datetime.now()
    tz = now.astimezone().tzname()
    time_string = now.strftime("%c ") + tz
    logging.info(const.HLINE)
    logging.info("Ends " + time_string)
    logging.info(const.HLINE)


def main():  # pylint: disable=missing-function-docstring
    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("Input")
    group.add_argument(
        "input",
        type=str,
        help="Input or structure file containing run configuration"
        " (*.ini, see ./examples/.)",
    )

    group = parser.add_argument_group("Modify structure")
    group.add_argument(
        "--from-species",
        "-f",
        nargs="*",
        type=str,
        help=(
            "Replace species/label, if found within the input"
            "structure, into species defined by `--to`."
            "Matches either `_atom_site_label` or `_atom_site_type_symbol` "
            "defined within the CIF."
        ),
        default=const.DEFAULT_FROM_SPECIES,
    )
    group.add_argument(
        "--to-species",
        "-t",
        nargs="*",
        type=str,
        help=(
            "Final chemical formula of the `--from` sites"
            "after substitution. Accepts various formats such as "
            "SmFe12, Sm0.5Fe0.5, Sm3+Sm2+2Fe3, etc. "
            "Different oxidation states are treated as distinct sites."
        ),
        default=const.DEFAULT_TO_SPECIES,
    )
    group.add_argument(
        "--scaling-matrix",
        "-s",
        nargs="*",
        type=str,  # To allow flexible separator
        help="Scales the unit cell. Comma-or-space-separated 3-or-9-digit values. "
        "Nine digit case is for the more general scaling matrix, "
        "whilst 3 is diagonal-only.",
        default=const.DEFAULT_SCALING_MATRIX_STR,
    )

    group = parser.add_argument_group("Run option")
    group.add_argument(
        "--symmetrize", action="store_true", help="Symmetrize input CIF."
    )
    group.add_argument(
        "--sample",
        default=const.DEFAULT_SAMPLE,
        help="From all generated order structure, sample N of them (no replacement).",
    )
    group.add_argument(
        "--symprec",
        type=float,
        default=const.DEFAULT_SYMPREC,
        help="Precision used in symmetry search.",
    )
    group.add_argument(
        "--angle-tolerance",
        type=float,
        default=const.DEFAULT_ANGLE_TOLERANCE,
        help="Angle tolerance for symmetry search.",
    )
    group.add_argument(
        "--dir-size",
        type=int,
        default=const.DEFAULT_DIR_SIZE,
        help="Divides output order structures in directories of DIR_SIZE files each.",
    )
    group.add_argument(
        "--count-only",
        action="store_true",
        help="Skip order structures generation, instead only"
        " prints the final number of structure.",
    )
    group.add_argument(
        "--write-symm",
        action="store_true",
        help="Write symmetries for the order structures (slower).",
    )
    group.add_argument(
        "--disable-progressbar",
        action="store_true",
        help="Disable progressbar",
    )
    group.add_argument(
        "--no-write",
        action="store_true",
        help="Do not write final structure",
    )
    group.add_argument(
        "--no-dmat",
        action="store_true",
        help="(devel) Alternative algorithm without distance matrix",
    )
    group.add_argument(
        "--t-kind",
        default="sum",
        choices=("sum", "sumfl", "det"),
        help="Type of T function (sum, determinant)"
    )
    args = parser.parse_args()
    const.DISABLE_PROGRESSBAR = args.disable_progressbar

    # Print header first for faster perceived response
    print_header()

    from .main import ScriptHelper  # pylint:disable=import-outside-toplevel

    if fnmatch.fnmatch(args.input, "*.ini"):
        # Note: when using *.ini, other arguments will be ignored.
        helper = ScriptHelper.from_file(args.input)
    else:
        # Trick to allow ",", ";", and whiteline as separator
        from_species = const.FLEXIBLE_SEPARATOR.split(",".join(args.from_species))
        to_species = const.FLEXIBLE_SEPARATOR.split(",".join(args.to_species))
        scaling_matrix = [
            int(x)
            for x in const.FLEXIBLE_SEPARATOR.split(",".join(args.scaling_matrix))
        ]
        from_species = list(filter(None, from_species))
        to_species = list(filter(None, to_species))
        scaling_matrix = list(filter(None, scaling_matrix))
        helper = ScriptHelper(
            structure_file=args.input,
            from_species=from_species,
            to_species=to_species,
            scaling_matrix=scaling_matrix,
            symmetrize=args.symmetrize,
            sample=args.sample,
            symprec=args.symprec,
            angle_tolerance=args.angle_tolerance,
            dir_size=args.dir_size,
            write_symm=args.write_symm,
            no_write=args.no_write,
            no_dmat=args.no_dmat,
            t_kind=args.t_kind,
        )
    helper.count()
    helper.save_modified_structure()
    if not args.count_only:
        helper.write()
    print_footer()
