# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
Cache SHRY instances for future run.
"""
import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure *= [1, 2, 1]
structure1 = structure.copy()
structure2 = structure.copy()
structure3 = structure.copy()
structure1.replace_species({'Fe1': 'Fe7Ti'})
# Higher Ti concentraction
structure2.replace_species({'Fe1': 'Fe3Ti'})
# Reverse of above
structure3.replace_species({'Fe1': 'FeTi3'})

os.makedirs("output1", exist_ok=True)
os.makedirs("output2", exist_ok=True)
os.makedirs("output3", exist_ok=True)

substitutor = Substitutor(structure1)
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output1/{i}.cif"
    cifwriter.write_file(filename=output_filename)

substitutor.structure = structure2
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output2/{i}.cif"
    cifwriter.write_file(filename=output_filename)

substitutor.structure = structure3
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output3/{i}.cif"
    cifwriter.write_file(filename=output_filename)
