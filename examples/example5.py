import os

from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure *= [1, 2, 1]
structure1 = structure.copy()
structure2 = structure.copy()
structure3 = structure.copy()
structure4 = structure.copy()
structure1.replace_species({'Fe1': 'Fe7Ti'})
structure2.replace_species({'Fe2': 'Fe6Ti2'})
structure3.replace_species({'Fe3': 'FeTi'})
structure4.replace_species({'Fe2': 'FeTi3'})

os.makedirs("output1", exist_ok=True)
os.makedirs("output2", exist_ok=True)
os.makedirs("output3", exist_ok=True)
os.makedirs("output4", exist_ok=True)

substitutor = Substitutor(structure1)
for i, _structure in enumerate(substitutor.order_structure()):
    output_filename = f"output1/{i}.cif"
    _structure.to(filename=output_filename)

substitutor.structure = structure2
for i, _structure in enumerate(substitutor.order_structure()):
    output_filename = f"output2/{i}.cif"
    _structure.to(filename=output_filename)

substitutor.structure = structure3
substitutor.sampled_indices = 50
for i, _structure in enumerate(substitutor.order_structure()):
    output_filename = f"output3/{i}.cif"
    _structure.to(filename=output_filename)

substitutor.structure = structure4
for i, _structure in enumerate(substitutor.order_structure()):
    output_filename = f"output4/{i}.cif"
    _structure.to(filename=output_filename)
