# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
Save Substitutor instance for later use.
"""
import os
import pickle

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)

# Replace Fe1 to Fe7Ti
structure.replace_species({'Fe1': 'Fe7Ti'})

# Enable caching
substitutor = Substitutor(structure, cache=True)

output_dir = "output"
os.makedirs(output_dir, exist_ok=True)
for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight"))):
    cifwriter = packet["cifwriter"]
    weight = packet["weight"]

    filename = f"cif_run0_i{i}w{weight}.cif"
    cifwriter.write_file(filename=os.path.join(output_dir, filename))

# Save pattern makers for later use.
with open('pg.pkl', 'wb') as f:
    pickle.dump(substitutor.pattern_makers, f)
