# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
Equivalent enumlib operations in SHRY
"""

from pymatgen.core import Structure
from pymatgen.command_line.enumlib_caller import EnumlibAdaptor
import shry
from shry import Substitutor

shry.const.DISABLE_PROGRESSBAR = True

# SmFe7Ti structure
# cif_file = 'SmFe7Ti.cif'
# structure = Structure.from_file(cif_file)

# PbSnTe structure
cif_file = 'PbSnTe.cif'
structure = Structure.from_file(cif_file)
structure *= (2, 2, 2)

# A) Enumerate and generate structure with enumlib
adaptor = EnumlibAdaptor(structure)
adaptor.run()
enumlib_structures = adaptor.structures
enumlib_num_structs = len(adaptor.structures)
print(f"Enumlib resulted in {enumlib_num_structs} structures")

# B) Enumerate and generate structure with shry
# 1) supercell-like: group equivalent sites
s = Substitutor(structure)
# Shry uses generator; below is to put the Structures into a list
shry_structures = [x for x in s.structure_writers()]
shry_num_structs = s.count()
print(f"SHRY (group equivalent sites) resulted in {shry_num_structs} structures")

# 2) enumlib-like: group same species
s = Substitutor(structure, groupby=lambda x: x.species)
shry_enumlib_structures = [x for x in s.structure_writers()]
shry_enumlib_num_structs = s.count()
print(f"SHRY (group same species) resulted in {shry_enumlib_num_structs} structures")
