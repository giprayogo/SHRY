import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe7Ti'})

# structure *= [1, 2, 1]
# 9-digits case
structure *= [[1, 0, 0], [0, 2, 0], [0, 0, 1]]
# structure *= np.array([[1, 0, 0], [0, 2, 0], [0, 0, 1]])

substitutor = Substitutor(structure, sample=50)
os.makedirs("output", exist_ok=True)
# A generator of Pymatgen's CifWriters
for i, cifwriter in enumerate(substitutor.cifwriters()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"

    cifwriter.write_file(filename=output_filename)
