#!/usr/bin/env python

import warnings
warnings.simplefilter('ignore')

import glob
import itertools
import math
import os, time
import re
import tqdm
import sys
from datetime import datetime
import zipfile
import shutil

import numpy as np
import pandas as pd
from pymatgen.io.cif import CifParser, CifFile, str2float
from shry import LabeledStructure
from shry.core import NeedSupercellError, Substitutor, PatchedSpacegroupAnalyzer
from pymatgen.util.coord import find_in_coord_list_pbc
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SpacegroupOperations
from pymatgen.symmetry.structure import SymmetrizedStructure

from ase import atoms
from ase.io import read, write

#TOLERANCE
SUPERCELL_TOLERANCE=1.0e-2 #angstrom
SHRY_TOLERANCE=1.0e-2      #angstrom
SHRY_ANGLE_TOLERANCE=5.0   #degree

# FLAGS
check_internal_redundancy=True
check_consistency=True
    
# Display formatting
LINEWIDTH = 88
HLINE = "-" * LINEWIDTH
TQDM_CONF = {
	"ncols": LINEWIDTH
}
DISABLE_PROGRESSBAR = False

def load_cif_dict(cif_filename):
    cif_file = CifFile.from_file(cif_filename)
    # assume single fragment
    key = list(cif_file.data.keys())[0]
    return cif_file.data[key].data

def get_coordinate_and_symbols(cif_filename, round=4):
    cif_dict = load_cif_dict(cif_filename)
    # Get sorting key from the first data
    symbols = [str(re.sub('[^a-zA-Z]+', '', i)) for i in cif_dict["_atom_site_type_symbol"]]
    x = [str2float(i) for i in cif_dict["_atom_site_fract_x"]]
    y = [str2float(i) for i in cif_dict["_atom_site_fract_y"]]
    z = [str2float(i) for i in cif_dict["_atom_site_fract_z"]]
    coords = [[float(x),float(y),float(z)] for x,y,z in zip(x, y, z)]
    coords = np.round(coords, round)
    condition = coords == 1.0000000
    coords[condition] = 0.0000000
    return {"symbols":symbols, "coords":coords}

def check_structure_consistency(str_cif_1, str_cif_2, symmops, write_matched_str=False, round=4):
    
    str_ref_1=get_coordinate_and_symbols(str_cif_1, round)
    str_ref_2=get_coordinate_and_symbols(str_cif_2, round)
    perms=[]
    for isym,symmop in enumerate(symmops):
        coords=str_ref_1["coords"]
        o_coords=symmop.operate_multi(coords)
        atol = 1e-8
        atol_coord = 1e-2
        stol_coord = 1e0
        step = 10
        is_up = True
        
        for _ in range(10):
            matches = [
                find_in_coord_list_pbc(coords, coord, atol=atol)
                for coord in o_coords
            ]
            #print (f"atol = {atol}")
            
            match_lengths = [len(x) for x in matches]
            if all(x == 1 for x in match_lengths):
                #print (f"Final atol = {atol}")
                break
            elif any(not x for x in match_lengths):
                #print ("elif any(not x for x in match_lengths)")
                # Too strict
                if not is_up:
                    step /= 2
                    is_up = True
                atol *= step
                #print (f"atol = {atol}")
            elif any(x > 1 for x in match_lengths):
                #print ("elif any(x > 1 for x in match_lengths)")
                # Too loose
                if is_up:
                    step /= 2
                    is_up = False
                atol /= step
                #print (f"atol = {atol}")
        else:
            raise RuntimeError("Failed to build symmetry list.")

        indices = [x[0] for x in matches]
        perms.append(indices)
    
    # check matching!!!!
    ref_coords=np.array(str_ref_2["coords"])
    ref_symbols=str_ref_2["symbols"]
    
    match_symobols=[]
    for iperm, perm in enumerate(perms):
        coords=np.array(str_ref_1["coords"])
        symbols=str_ref_1["symbols"]
        #perm_coords=np.array([coords[i] for i in perm])
        perm_coords=coords
        perm_symbols=[symbols[i] for i in perm]
        
        #lexsort!! ref1
        lex_perm_coords=np.array([perm_coords[i] for i in np.lexsort(perm_coords.T)])
        lex_perm_symbols=[perm_symbols[i] for i in np.lexsort(perm_coords.T)]

        #lexsort!! ref2
        lex_ref_coords=np.array([ref_coords[i] for i in np.lexsort(ref_coords.T)])
        lex_ref_symbols=[ref_symbols[i] for i in np.lexsort(ref_coords.T)]
        
        coords_sum=np.sum(np.abs(np.array(lex_ref_coords)-np.array(lex_perm_coords)))
        coords_sum_each=np.abs([np.array(ref_coord) - np.array(perm_coord) for ref_coord, perm_coord in zip(lex_ref_coords,lex_perm_coords)])
        if not coords_sum < stol_coord:
            #print("Something wrong in the applied symmetry operations!!!")
            #print("This usually happens due to inconsistency in the significant digits")
            #print("Try to charge round in np.round(coords, round)")
            #print(f"str_cif_1 = {str_cif_1}")
            #print(f"str_cif_2 = {str_cif_2}")
            #print("coords_sum is not zero")
            #print(f"coords_sum={coords_sum}")
            #print("------------------------------")
            #print(lex_perm_symbols)
            #print(lex_ref_symbols)
            #print("------------------------------")
            #print(lex_perm_coords)
            #print(lex_ref_coords)
            #print("------------------------------")
            #print(coords_sum_each)
            raise ValueError
        if not any([ np.sum(coord) < atol_coord for coord in coords_sum_each]):
            #print("Something wrong in the applied symmetry operations!!!")
            #print("This usually happens due to inconsistency in the significant digits")
            #print("Try to charge round in np.round(coords, round)")
            #print(f"str_cif_1 = {str_cif_1}")
            #print(f"str_cif_2 = {str_cif_2}")
            #print("There is a finte element in coords_sum_each")
            #print(f"coords_sum_each={coords_sum_each}")
            raise ValueError
        match_symobols.append(lex_perm_symbols==lex_ref_symbols)
        
        if lex_perm_symbols==lex_ref_symbols and write_matched_str:
            print("A redundant cif file pair has been found!!")
            #print(f"num symop. {len(symmops)}")
            #print(f"coords_sum={coords_sum}")
            #print(f"coords_sum_each={coords_sum_each}")
            #print(f"symmops[{iperm}]")
            print(f"str_cif_1 = {str_cif_1}")
            print(f"str_cif_2 = {str_cif_2}")
            print(f"For the sym. op.")
            print(f"rot:")
            print(symmops[iperm].rotation_matrix)
            print(f"tra:")
            print(symmops[iperm].translation_vector)
            
            #ASE read and write
            atom_1=read(str_cif_1)
            atom_2=read(str_cif_2)
            
            atom_1_basename = os.path.basename(str_cif_1).replace(".cif","")
            atom_2_basename = os.path.basename(str_cif_2).replace(".cif","")
            atom_1_dirname = os.path.dirname(str_cif_1)
            cif_name=os.path.join(atom_1_dirname, f"v_{atom_1_basename}_to_{atom_2_basename}.cif")
            
            coord1=atom_1.get_scaled_positions()
            o_coords=symmops[iperm].operate_multi(coords)
            atom_1.set_scaled_positions(o_coords)
            write(cif_name, atom_1)
            print(f"Written as {cif_name}")
    
    return any(match_symobols)

def main():
    
    timestamp = int(datetime.now().timestamp())
    df = pd.read_excel(sys.argv[1])
    rows, _ = df.shape

    # Things for string matching
    coll_code_regex = re.compile(r"cod.*[0-9]+")
    number_sign_regex = re.compile(r"[0-9\+-]+")
    supercell_cif_glob = "supercell-*{}*/supercell*.cif"
    shry_cif_glob = "shry-*{}*/slice*/*.cif"
    np_round_list=[4,3] # for lexsort!!
    
    #df["Checked"] = False
    #df["Note"] = ""
    df.drop(columns="T_supercell_real", inplace=True, errors="ignore")
    df.drop(columns="T_supercell_user", inplace=True, errors="ignore")
    df.drop(columns="T_supercell_sys", inplace=True, errors="ignore")
    df.drop(columns="T_shry_real", inplace=True, errors="ignore")
    df.drop(columns="T_shry_user", inplace=True, errors="ignore")
    df.drop(columns="T_shry_sys", inplace=True, errors="ignore")
    
    try:
        # About this try block: save DataFrame no matter what
        for row, zipped in enumerate(
            zip(
                df["Supercell"],
                df["File"],
                df["Checked"],
            )
        ):
        
            supercell, reference_filename, checked = zipped
            # Skip checked structures.
            if checked:
                continue
            print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"Working on {reference_filename} ({row+1}/{rows})")
    
            # Every required information from the DataFrame
            reference_filename = reference_filename.replace(".cif", "_partial.cif")
            reference_dirname = os.path.dirname(reference_filename)
            reference_basename = os.path.basename(reference_filename)
            supercell_array = np.array(list(map(int, supercell.split("x"))))
            #print(f"supercell = {supercell_array}")
            
            reference_structure = LabeledStructure.from_file(reference_filename)
            sga=PatchedSpacegroupAnalyzer(reference_structure, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE)
            symmops = sga.get_symmetry_operations()
            #print("---------------------------------")
            #print("symmetry/unitcell")
            #print(symmops)
            #print(len(symmops))
            #print(reference_structure)
            #time.sleep(10)
            
            supercell_array = np.array(list(map(int, supercell.split("x"))))
            reference_structure *= supercell_array
            sga=PatchedSpacegroupAnalyzer(reference_structure, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE)
            symmops = sga.get_symmetry_operations()
            #print("---------------------------------")
            #print("symmetry/super cell")
            #print(symmops)
            #print(f"len(symmops) = {len(symmops)}")
            #print(reference_structure)
            #sys.exit()
            
            # NOTE: All CIFs should contain CollCode string
            coll_codes = coll_code_regex.findall(reference_basename)
            if len(coll_codes) != 1:
                df.at[row, "Note"] = (
                    df.at[row, "Note"] + ", invalid reference CIF filename"
                )
                continue
            coll_code = coll_codes[0]  # long
            
            # extract zips!! SHRY and Supercell
            #print("Extracting zip dirs")
            zip_files=glob.glob(
                os.path.join(reference_dirname, "*_{}_*.zip".format(coll_code))
            )
            
            for zip_file in zip_files:
                filebase = os.path.basename(zip_file)
                filebasenozip = filebase.replace(".zip", "")
                with zipfile.ZipFile(zip_file, 'r') as f:
                    f.extractall(os.path.join(reference_dirname, filebasenozip))
            
            # set SHRY and Supercell file names
            supercell_cifs = glob.glob(
                os.path.join(reference_dirname, supercell_cif_glob.format(coll_code))
            )
            shry_cifs = glob.glob(
                os.path.join(reference_dirname, shry_cif_glob.format(coll_code))
            )
            
            
            # Store all structures
            #print(supercell_cifs[0])
            #print(shry_cifs[0])
            
            # remove!!!
            #supercell_structures=[get_coordinate_and_symbols(cif) for cif in supercell_cifs]
            #shry_structures=[get_coordinate_and_symbols(cif) for cif in shry_cifs]
            
            if check_internal_redundancy:
                # Supercell_structure check:
                print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print("Supercell: Check if all the generated structures are symmetrically inequivalent.")
                match_flag_list=[]
                match_cif_1=[]
                match_cif_2=[]
                all_combinations=list(itertools.combinations(range(len(supercell_cifs)), 2))
                for cif_set in tqdm.tqdm(
                    all_combinations,
                    desc="Supercell: all str. combinations",
                    **TQDM_CONF,
                    disable=DISABLE_PROGRESSBAR,
                ):
                    str_cif_1=supercell_cifs[cif_set[0]]
                    str_cif_2=supercell_cifs[cif_set[1]]
                    
                    match_flag=check_structure_consistency(str_cif_1, str_cif_2, symmops, True)
                    match_flag_list.append(match_flag)
                    if match_flag:
                        match_cif_1.append(supercell_cifs[cif_set[0]])
                        match_cif_2.append(supercell_cifs[cif_set[1]])
                    #print(f"Str:{str_set[0]} v.s. Str.:{str_set[1]} = {match_flag}")
                    #sys.exit()
                
                if any(match_flag_list):
                    print(f"Supercell, {reference_filename}, the symmetrically-inequivalent test failed.")
                    pairs_pd=pd.DataFrame(np.array([match_cif_1, match_cif_2]).T, columns=["supercell1_cif","supercell2_cif"])
                    pairs_pd.to_csv(os.path.join(reference_dirname, f"{reference_basename}_supercell_matched.csv"))
                
                else:
                    print(f"Supercell, {reference_filename}, the symmetrically-inequivalent test succeeded.")
        
                # SHRY_structure check:
                print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print("SHRY: Check if all the generated structures are symmetrically inequivalent.")
                match_flag_list=[]
                all_combinations=list(itertools.combinations(range(len(shry_cifs)), 2))
                for cif_set in tqdm.tqdm(
                    all_combinations,
                    desc="SHRY: all str. combinations",
                    **TQDM_CONF,
                    disable=DISABLE_PROGRESSBAR,
                ):
                    str_cif_1=shry_cifs[cif_set[0]]
                    str_cif_2=shry_cifs[cif_set[1]]
                    
                    match_flag=check_structure_consistency(str_cif_1, str_cif_2, symmops, True)
                    match_flag_list.append(match_flag)
                    #print(f"Str:{str_set[0]} v.s. Str.:{str_set[1]} = {match_flag}")
                    #sys.exit()
                
                if any(match_flag_list):
                    print(f"SHRY, {reference_filename}, the symmetrically-inequivalent test failed.")
                else:
                    print(f"SHRY, {reference_filename}, the symmetrically-inequivalent test succeeded.")
    
            if check_consistency:
                
                """
                # Noted already by previous scripts;
                if len(supercell_cifs) != len(shry_cifs):
                    print("Mismatch in number of irreducible structure")
                    df.at[row, "Note"] = df.at[row, "Note"] + ", cif number mismatch"
                    sys.exit()
                    continue
                """
                
                # Consistency check between SHRY and Supercell!!:
                print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print("SHRY v.s. Supercell: Check if all the generated structures are consistent between SHRY and Supercell!")
                shry_match_pair_list=[]
                supercell_match_pair_list=[]
                non_match_shry_index_list=[]
                supercell_index_list=[i for i in range(len(supercell_cifs))]
                shry_index_list=[i for i in range(len(shry_cifs))]
                
                for shry_index in tqdm.tqdm(
                    shry_index_list,
                    desc="SHRY v.s. Supercell",
                    **TQDM_CONF,
                    disable=DISABLE_PROGRESSBAR,
                ):
                
                    match_flag=False
                    #print(supercell_index_list)
                    for supercell_index in supercell_index_list:
                        for np_round in np_round_list:
                            try:
                                match_flag=check_structure_consistency(
                                shry_cifs[shry_index], 
                                supercell_cifs[supercell_index], 
                                symmops,
                                False,
                                np_round)
                                #print(f"np_round = {np_round} is succesfull")
                                break
                            except ValueError:
                                #print(f"np_round = {np_round} fails")
                                continue
                        #print(match_flag)
                        if match_flag:
                            supercell_index_list.remove(supercell_index)
                            shry_match_pair_list.append(shry_index)
                            supercell_match_pair_list.append(supercell_index)
                            break
                    if not match_flag:
                        print("For SHRY output: {shry_cifs[shry_index]}")
                        print("There is no equivalent structure in the supercell outputs!!")
                        non_match_shry_index_list.append(shry_index)
                        raise
                
                pairs_pd=pd.DataFrame(np.array([[shry_cifs[i] for i in shry_match_pair_list], [supercell_cifs[i] for i in supercell_match_pair_list]]).T, columns=["shry_cif","supercell_cif"])
                pairs_pd.to_csv(os.path.join(reference_dirname, f"{reference_basename}_pairs.csv"))
                if len(non_match_shry_index_list) !=0 or len(supercell_index_list) !=0:
                    df.at[row, "Note"] = (df.at[row, "Note"] + f"non_match_shry_cifs: {[shry_cifs[i] for i in non_match_shry_index_list]}, non_match_supercell_cifs:{[supercell_cifs[i] for i in supercell_index_list]}")
            
            
            #Finally
            df.at[row, "Checked"] = True
            
            #finally, delete the extracted dirs
            #print("Delete zip-extracted dirs")
            for zip_file in zip_files:
                filebase = os.path.basename(zip_file)
                filebasenozip = filebase.replace(".zip", "")
                shutil.rmtree(os.path.join(reference_dirname,filebasenozip))

    except Exception as e:
        #print(f"type:{str(type(e))}")
        #print(e.args)
        import traceback
        print(traceback.format_exc())
        
    finally:
        #print(df)
        df.to_excel(f"fullcheck_{timestamp}.xls", index=False)
        
if __name__ == "__main__":
    main()
