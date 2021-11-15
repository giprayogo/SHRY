import warnings
warnings.simplefilter('ignore')

from ase import Atoms
from ase.io import read, write
from ase.spacegroup import Spacegroup
import pymatgen
from pymatgen.io.cif import CifParser
import os, sys
import glob

root_dir=os.getcwd()
sum_cif_files=0

for i in range(1,231):
    sg = Spacegroup(i)
    print('------------------------------')
    print('Space group', sg.no, sg.symbol)
    print('------------------------------')
    
    sg_dir=os.path.join(root_dir, f"SG{str(sg.no)}")
    
    cif_files=glob.glob(os.path.join(sg_dir, '*[0-9].cif'))
    
    print(f"The number of cif files is {len(cif_files)}")
    sum_cif_files+=len(cif_files)
    
    for cif in cif_files:
        parser = CifParser(cif)
        cif_instance = parser.get_structures()[0]
        reduced_str = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(cif_instance)
        #reduced_str = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(cif, symprec=0.01, angle_tolerance=1.0)
        #print("print structure")
        #print(cif)
        #print("------------------------------------------")
        #print("This is the space group symmbol")
        print(f"{os.path.basename(cif)}: Space group= {reduced_str.get_space_group_symbol()}, Space group number= {reduced_str.get_space_group_number()}")
        if reduced_str.get_space_group_number() != i:
            print("Error!! SGs are not consistent.")
            sys.exit()
        #sym_str = reduced_str.get_symmetrized_structure()
        #print("This is the symmetrized structure")
        #print(sym_str)
        #print(f"=========================================")
        #print(f" ")
        #print(f" ")
    print(f" ")

print("All the symmetry check has passed.")

print(f"Total number of cif files is {sum_cif_files}")