import os
from re import sub

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe7Ti'})

substitutor = Substitutor(structure)
substitutor.sampled_indices = 50
os.makedirs("output", exist_ok=True)
# A generator of Structure instances.
for i, _structure in enumerate(substitutor.order_structure()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"

    _structure.to(filename=output_filename)
