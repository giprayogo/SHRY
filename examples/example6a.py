import pickle

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe9Ti3'})

substitutor = Substitutor(structure)
substitutor.make_patterns()

# Save pattern makers for later use.
with open('pg.pkl', 'wb') as f:
    pickle.dump(substitutor.pattern_makers, f)
