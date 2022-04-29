# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
LabeledStructure examples with replace_species
"""
import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)

# Replace Fe1 in the CIF to Fe3Ti
structure.replace_species({'Fe1': 'Fe3Ti'})
# Replace Fe2 in the CIF to Fe7Ti,
# and Fe3 to full Ti
structure.replace_species({'Fe2': 'Fe7Ti', 'Fe3': 'Ti'})

# Change cell shape for a little variation
structure *= ((2, 0, 1), (0, 1, 0), (1, 0, 1))

substitutor = Substitutor(structure)

output_dir = "output"
os.makedirs(output_dir, exist_ok=True)
for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight", "letter"))):
    cifwriter = packet["cifwriter"]
    weight = packet["weight"]

    # Print the configuration letter
    print(packet["letter"])

    filename = f"cif_i{i}w{weight}.cif"
    cifwriter.write_file(filename=os.path.join(output_dir, filename))
