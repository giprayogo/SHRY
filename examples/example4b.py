# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
Load saved Substitutor instance for further substitution.
"""
import pickle
import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure2 = structure.copy()

structure.replace_species({'Fe': 'Fe3Ti'})
structure2.replace_species({'Fe': 'FeTi3'})

# Enable caching here too
substitutor = Substitutor(structure, cache=True)

# Load pattern from previous Substitutor
with open('pg.pkl', 'rb') as f:
    pattern_makers = pickle.load(f)
substitutor.pattern_makers = pattern_makers

output_dir = "output"
os.makedirs(output_dir, exist_ok=True)
for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight"))):
    cifwriter = packet["cifwriter"]
    weight = packet["weight"]

    filename = f"cif_run1_i{i}w{weight}.cif"
    cifwriter.write_file(filename=os.path.join(output_dir, filename))

# Reuse the Substitutor instance
substitutor.structure = structure2
for i, packet in enumerate(substitutor.quantities(("cifwriter", "weight"))):
    cifwriter = packet["cifwriter"]
    weight = packet["weight"]

    filename = f"cif_run2_i{i}w{weight}.cif"
    cifwriter.write_file(filename=os.path.join(output_dir, filename))

# Update the previous pattern makers.
with open('pg.pkl', 'wb') as f:
    pickle.dump(substitutor.pattern_makers, f)
