import pickle
import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe9Ti3'})

os.makedirs("output", exist_ok=True)

substitutor = Substitutor(structure)

# Load previous patterns
with open('pg.pkl', 'rb') as f:
    pattern_makers = pickle.load(f)
substitutor.pattern_makers = pattern_makers

# Refresh Substitutor state
substitutor.get_ap()

# Substitute other structures
for i, cifwriter in enumerate(substitutor.cifwriters()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"
    cifwriter.write_file(filename=output_filename)
