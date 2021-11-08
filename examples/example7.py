import io
import os
import shutil

from shry import LabeledStructure, Substitutor
from pymatgen.analysis.ewald import EwaldSummation

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe2+7Ti3+'})
structure *= [1, 2, 1]

substitutor = Substitutor(structure)
substitutor.sampled_indices = 50
os.makedirs("output", exist_ok=True)

logio = io.StringIO()

# A generator of Structure instances.
for i, _structure in enumerate(substitutor.order_structure()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"
    
    ewald_sum = EwaldSummation(_structure).total_energy
    print(f"{i} {ewald_sum}", file=logio)
    _structure.to(filename=output_filename)

with open("output/log.txt", "w") as f:
    logio.seek(0)
    shutil.copyfileobj(logio, f)
