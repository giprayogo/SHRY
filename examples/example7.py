import io
import os
import shutil

from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.io.cif import CifParser
from shry import LabeledStructure, Substitutor

cif_file = 'SmFe12.cif'
structure = LabeledStructure.from_file(cif_file)
structure.replace_species({'Fe': 'Fe2+7Ti3+'})
# Also add oxidation state to the Sm (any)
structure.replace_species({'Sm': 'Sm+'})

substitutor = Substitutor(structure, sample=50)
substitutor.make_patterns()
os.makedirs("output", exist_ok=True)

logio = io.StringIO()

# A generator for Pymatgen's CifWriters
for i, cifwriter in enumerate(substitutor.cifwriters()):
    # Some naming logic.
    output_filename = f"output/{i}.cif"

    # Create Structure from the CifWriter
    cifparser = CifParser.from_string(str(cifwriter))
    structure = cifparser.get_structures(primitive=False)[0]

    ewald_sum = EwaldSummation(structure).total_energy
    print(f"{i} {ewald_sum}", file=logio)
    structure.to(filename=output_filename)

# Save the Ewald summations
with open("output/log.txt", "w") as f:
    logio.seek(0)
    shutil.copyfileobj(logio, f)
