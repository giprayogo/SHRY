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

"""
Shared constants
"""
import re

# Default input arguments
DEFAULT_FROM_SPECIES = ()
DEFAULT_TO_SPECIES = ()
DEFAULT_SCALING_MATRIX = (1, 1, 1)
DEFAULT_SCALING_MATRIX_STR = "1 1 1"
DEFAULT_SYMMETRIZE = False
DEFAULT_SAMPLE = None
DEFAULT_SYMPREC = 1.0e-2
DEFAULT_ANGLE_TOLERANCE = 5.0
DEFAULT_DIR_SIZE = 10000
DEFAULT_WRITE_SYMM = False
DEFAULT_NO_WRITE = False
DEFAULT_NO_DMAT = False
DEFAULT_T_KIND = "sum"

# Symmetry list building parameters
FIND_COORDS_LIST_ATOL_INIT = 1e-8
FIND_COORDS_LIST_ATOL_MUL_STEP_INIT = 10
FIND_COORDS_LIST_ATTEMPTS = 10

# Global state configuration
DISABLE_PROGRESSBAR = False

# Flexible input format
FLEXIBLE_SEPARATOR = re.compile(r"[;,\s]+")

# For matching symmetric sites
MAX_ATOL = 1.0e-5
MAX_IRREDUCIBLE = 1e9

# Display formatting
LINEWIDTH = 88
HLINE = "-" * LINEWIDTH
TQDM_CONF = {
	"ncols": LINEWIDTH
}
