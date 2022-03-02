import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe1': 'Fe3Ti'})
structure.replace_species({'Fe2': 'Fe3Ti', 'Fe3': 'Fe3Ti'})

substitutor = Substitutor(structure)
substitutor.get_ap()
os.makedirs("output", exist_ok=True)
# A generator for Pymatgen's CifWriters
for i, cifwriter in enumerate(substitutor.cifwriters()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"

    cifwriter.write_file(filename=output_filename)
