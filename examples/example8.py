import os

from pymatgen.core import Structure
from shry import Substitutor

cif_file = 'SmFe7Ti.cif'
structure = Structure.from_file(cif_file)

# Loops over Structure.sites
def groupby(site):
    return site.species

os.makedirs("output1", exist_ok=True)
os.makedirs("output2", exist_ok=True)
substitutor = Substitutor(structure)
print(structure)

# Equivalent to Fe1 -> Fe7Ti, Fe2 -> Fe7Ti, Fe3 -> Fe7Ti
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output1/{i}.cif"
    cifwriter.write_file(filename=output_filename)

# Equivalent to Fe -> Fe7Ti
substitutor.groupby = groupby
for i, cifwriter in enumerate(substitutor.cifwriters()):
    output_filename = f"output2/{i}.cif"
    cifwriter.write_file(filename=output_filename)
