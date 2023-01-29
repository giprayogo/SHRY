# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Basic use: read a CIF file and
write the substituted CIFs with weight.
"""
import os
from pymatgen.core import Structure
from shry import Substitutor

# Read a CIF file
cif_file = "SmFe7Ti.cif"
structure = Structure.from_file(cif_file)

substitutor = Substitutor(structure)

# Save CIF files
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)
for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight"))):
    cifwriter = packet["cifwriter"]
    weight = packet["weight"]

    filename = f"cif_i{i}w{weight}.cif"
    cifwriter.write_file(filename=os.path.join(output_dir, filename))
