# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Shared constants
"""

# python modules
import re
import datetime

# Default input arguments
# TODO: Refactor together with ScriptHelper
DEFAULT_FROM_SPECIES = ()
DEFAULT_TO_SPECIES = ()
DEFAULT_SCALING_MATRIX = (1, 1, 1)
DEFAULT_SCALING_MATRIX_STR = "1 1 1"
DEFAULT_SYMMETRIZE = False
DEFAULT_SAMPLE = None
DEFAULT_SYMPREC = 1.0e-2
DEFAULT_ATOL = 1.0e-5
DEFAULT_ANGLE_TOLERANCE = 5.0
DEFAULT_DIR_SIZE = 10000
DEFAULT_WRITE_SYMM = False
DEFAULT_WRITE_EWALD = False
DEFAULT_MAX_EWALD = None
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
TQDM_CONF = {"ncols": LINEWIDTH}

# Date and time, also used for (default) random seed
NOW = datetime.datetime.now()
DEFAULT_SEED = int(NOW.timestamp())
