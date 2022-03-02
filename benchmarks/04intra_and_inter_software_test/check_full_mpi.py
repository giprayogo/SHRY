#!/usr/bin/env python

import warnings
warnings.simplefilter('ignore')

import os
os.environ["OMP_NUM_THREADS"] = "1" # setenv OMP_NUM_THREADS 1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # setenv OPENBLAS_NUM_THREADS 1
os.environ["MKL_NUM_THREADS"] = "1" # setenv MKL_NUM_THREADS 1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # setenv VECLIB_MAXIMUM_THREADS 1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # setenv NUMEXPR_NUM_THREADS 1

import time
import glob
import itertools
import math
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

# Display formatting
LINEWIDTH = 88
HLINE = "-" * LINEWIDTH
TQDM_CONF = {
	"ncols": LINEWIDTH
}
DISABLE_PROGRESSBAR = False

# set mpi env
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
from mpi4py import MPI

#MPI init
MPI.Init()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
root = 0
    
def mpiabort_excepthook(type, value, traceback):
    comm.Abort()
    sys.__excepthook__(type, value, traceback)
    
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
            match_lengths = [len(x) for x in matches]
            if all(x == 1 for x in match_lengths):
                break
            elif any(not x for x in match_lengths):
                # Too strict
                if not is_up:
                    step /= 2
                    is_up = True
                atol *= step
            elif any(x > 1 for x in match_lengths):
                # Too loose
                if is_up:
                    step /= 2
                    is_up = False
                atol /= step
        else:
            raise RuntimeError(f"rank={rank}: Failed to build symmetry list.")

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
            #print("coords_sum is not zero")
            #print(f"coords_sum={coords_sum}")
            raise AssertionError(f"rank={rank}: Something wrong in the applied symmetry operations!!!")
        if not any([ np.sum(coord) < atol_coord for coord in coords_sum_each]):
            #print("Something wrong in the applied symmetry operations!!!")
            #print("This usually happens due to inconsistency in the significant digits")
            #print("Try to charge round in np.round(coords, round)")
            #print("there is a finte element in coords_sum_each")
            #print(f"coords_sum_each={coords_sum_each}")
            raise AssertionError(f"rank={rank}: Something wrong in the applied symmetry operations!!!")
        match_symobols.append(lex_perm_symbols==lex_ref_symbols)

        if lex_perm_symbols==lex_ref_symbols and write_matched_str:
            print(f"rank={rank}: A redundant cif file pair has been found!!")
            print(f"rank={rank}: str_cif_1 = {str_cif_1}")
            print(f"rank={rank}: str_cif_2 = {str_cif_2}")
            print(f"rank={rank}: For the sym. op.")
            print(f"rank={rank}: rot:")
            print(symmops[iperm].rotation_matrix)
            print(f"rank={rank}: tra:")
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
            print(f"rank={rank}: Written as {cif_name}")

    return any(match_symobols)

##################################
# main():
##################################

# FLAGS
check_internal_redundancy=True
check_consistency=True
#print(check_internal_redundancy)
#print(check_consistency)

#MPI read only for root
if rank==root:
    timestamp = int(datetime.now().timestamp())
    df = pd.read_excel(sys.argv[1])
    rows, _ = df.shape
    #df["Checked"] = False
    #df["Note"] = ""
    df.drop(columns="T_supercell_real", inplace=True, errors="ignore")
    df.drop(columns="T_supercell_user", inplace=True, errors="ignore")
    df.drop(columns="T_supercell_sys", inplace=True, errors="ignore")
    df.drop(columns="T_shry_real", inplace=True, errors="ignore")
    df.drop(columns="T_shry_user", inplace=True, errors="ignore")
    df.drop(columns="T_shry_sys", inplace=True, errors="ignore")
else:
    timestamp=None
    df=None

#mpi bcast
timestamp = comm.bcast(timestamp, root=root)
df = comm.bcast(df, root=root)

# Things for string matching
coll_code_regex = re.compile(r"cod.*[0-9]+")
number_sign_regex = re.compile(r"[0-9\+-]+")
supercell_cif_glob = "supercell-*{}*/supercell*.cif"
shry_cif_glob = "shry-*{}*/slice*/*.cif"
np_round_list=[4,3] # for lexsort!!

#####################################################
# Big for loop:
#####################################################

# About this try block: save DataFrame no matter what
for row, zipped in enumerate(
    zip(
        df["Supercell"],
        df["File"],
        df["Checked"],
    )
):

    try:
    
        # mpiabort test!
        #if rank == 1:
        #    print("rank1 raise error!!")
        #    raise AssertionError(f"rank={rank}: Failure")
        
        if rank==root:
    
            supercell, reference_filename, checked = zipped
            # Skip checked structures.
            if checked:
                continue
            print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print(f"Working on {reference_filename} ({row+1}/{rows})",flush=True)
    
            # Every required information from the DataFrame
            reference_filename = reference_filename.replace(".cif", "_partial.cif")
            reference_dirname = os.path.dirname(reference_filename)
            reference_basename = os.path.basename(reference_filename)
            supercell_array = np.array(list(map(int, supercell.split("x"))))
            #print(f"supercell = {supercell_array}")
            
            #print("LabeledStructure: start")
            reference_structure = LabeledStructure.from_file(reference_filename)
            #print("LabeledStructure: end")
            #sga=PatchedSpacegroupAnalyzer(reference_structure, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE)
            #symmops = sga.get_symmetry_operations()
            #print("---------------------------------")
            #print(symmops)
            #print(len(symmops))
            
            supercell_array = np.array(list(map(int, supercell.split("x"))))
            #print("supercell_array: start")
            reference_structure *= supercell_array
            #print("supercell_array: end")
            #print("PatchedSpacegrouptAnalyzer: start")
            sga=PatchedSpacegroupAnalyzer(reference_structure, symprec=SHRY_TOLERANCE, angle_tolerance=SHRY_ANGLE_TOLERANCE)
            symmops = sga.get_symmetry_operations()
            #print("PatchedSpacegrouptAnalyzer: end")
            #print("---------------------------------")
            #print(symmops)
            #print(len(symmops))
            #print(f"sys.getsizeof(reference_structure) = {sys.getsizeof(reference_structure)}")
            #print(f"The symmetry operations of the ref. structure are retrieved.")
            #sys.exit()
            
            # NOTE: All CIFs should contain CollCode string
            coll_codes = coll_code_regex.findall(reference_basename)
            if len(coll_codes) != 1:
                df.at[row, "Note"] = (
                    df.at[row, "Note"] + ", invalid reference CIF filename"
                )
                raise AssertionError(f"rank={rank}: invalid reference CIF filename")
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
            
            # Noted already by previous scripts;
            if len(supercell_cifs) != len(shry_cifs):
                print("Mismatch in number of irreducible structure")
                df.at[row, "Note"] = df.at[row, "Note"] + ", cif number mismatch"
                raise AssertionError(f"rank={rank}: Mismatch in number of irreducible structure")
            
            # Store all structures
            #print(supercell_cifs[0])
            #print(shry_cifs[0])
            
            # !! should be a problem in terms of memory.
            #supercell_structures=[get_coordinate_and_symbols(cif) for cif in supercell_cifs]
            #shry_structures=[get_coordinate_and_symbols(cif) for cif in shry_cifs]
        
        else:
            symmops=None
            supercell_cifs=None
            shry_cifs=None
            #supercell_structures=None
            #shry_structures=None
        
        #mpi bcast
        symmops = comm.bcast(symmops, root=root)
        supercell_cifs = comm.bcast(supercell_cifs, root=root)
        shry_cifs = comm.bcast(shry_cifs, root=root)
        #supercell_structures = comm.bcast(supercell_structures, root=root)
        #shry_structures = comm.bcast(shry_structures, root=root)
        
        #supercel check MPI loop
        if check_internal_redundancy:
            # Supercell_structure check:
            if (rank == root): print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            if (rank == root): print("Supercell: Check if all the generated structures are symmetrically inequivalent.")
            match_flag_list=[]
            all_combinations=list(itertools.combinations(range(len(supercell_cifs)), 2))
            split_combinations=np.array_split(all_combinations, size)[rank]
            if (rank == root): print(f"MPI: Num. Core = {size}")
            if (rank == root): print(f"Total num:{len(all_combinations)}.",flush=True)

            #tdqm + MPI is simply achieved by all_reduce.
            #for str_set in tqdm.tqdm(
            #    split_combinations,
            #    desc="Supercell: all str. combinations",
            #    **TQDM_CONF,
            #    disable=DISABLE_PROGRESSBAR,
            #):
            
            for i, cif_set in enumerate(split_combinations):
                #if (i%100000==0): print(f"rank={rank}, {i}/{len(split_combinations)}.",flush=True)
                str_cif_1=supercell_cifs[cif_set[0]]
                str_cif_2=supercell_cifs[cif_set[1]]
                for np_round in np_round_list:
                    try:
                        match_flag=check_structure_consistency(str_cif_1, str_cif_2, symmops, True, np_round)
                        break
                    except AssertionError:
                        match_flag=False
                        continue
                match_flag_list.append(match_flag)
                #print(f"Str:{str_set[0]} v.s. Str.:{str_set[1]} = {match_flag}")
                del str_cif_1
                del str_cif_2
                #sys.exit()
            
            if (rank == root): print()
            #print(f"rank={rank}, match_flag_list={match_flag_list}")
            
            # MPI Gather here. Collect local array sizes using the high-level mpi4py gather
            sendbuf=np.array(match_flag_list)
            sendcounts = np.array(comm.gather(len(sendbuf), root))

            if rank == root:
                #print("sendcounts: {}, total: {}".format(sendcounts, sum(sendcounts)))
                recvbuf = np.empty(sum(sendcounts), dtype=bool)
            else:
                recvbuf = None
            
            comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)
            
            if rank == root:
                match_flag_list=recvbuf
                #print(match_flag_list)
                if any(match_flag_list):
                    print(f"Supercell, {reference_filename}, the symmetrically-inequivalent test failed.")
                    df.at[row, "Note"] = df.at[row, "Note"] + f"the symmetrically-inequivalent test: supercell failed."
                    
                else:
                    print(f"Supercell, {reference_filename}, the symmetrically-inequivalent test succeeded.")
        
            comm.Barrier()
        
            # Finally, update df
            if rank==root:
                df.to_excel(f"fullcheck_{timestamp}.xls", index=False)
            else:
                df=None
            #mpi bcast
            df = comm.bcast(df, root=root)
        
        
        #shry check MPI loop
        if check_internal_redundancy:
            # SHRY_structure check:
            if (rank == root): print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            if (rank == root): print("SHRY: Check if all the generated structures are symmetrically inequivalent.")
            match_flag_list=[]
            all_combinations=list(itertools.combinations(range(len(shry_cifs)), 2))
            split_combinations=np.array_split(all_combinations, size)[rank]
            if (rank == root): print(f"MPI: Num. Core = {size}")
            if (rank == root): print(f"Total num:{len(all_combinations)}.",flush=True)
            #for str_set in tqdm.tqdm(
            #    split_combinations,
            #    desc="SHRY: all str. combinations",
            #    **TQDM_CONF,
            #    disable=DISABLE_PROGRESSBAR,
            #):
            
            for i, cif_set in enumerate(split_combinations):
                #if (i%100000==0): print(f"rank={rank}, {i}/{len(split_combinations)}.",flush=True)
                str_cif_1=shry_cifs[cif_set[0]]
                str_cif_2=shry_cifs[cif_set[1]]
                
                for np_round in np_round_list:
                    try:
                        match_flag=check_structure_consistency(str_cif_1, str_cif_2, symmops, True, np_round)
                        break
                    except AssertionError:
                        match_flag=False
                        continue
                match_flag_list.append(match_flag)
                #print(f"Str:{str_set[0]} v.s. Str.:{str_set[1]} = {match_flag}")
                del str_cif_1
                del str_cif_2
                #sys.exit()
            
            if (rank == root): print()

            # MPI Gather here. Collect local array sizes using the high-level mpi4py gather
            sendbuf=np.array(match_flag_list)
            sendcounts = np.array(comm.gather(len(sendbuf), root))

            if rank == root:
                #print("sendcounts: {}, total: {}".format(sendcounts, sum(sendcounts)))
                recvbuf = np.empty(sum(sendcounts), dtype=bool)
            else:
                recvbuf = None
            
            comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)
            
            if rank == root:
                match_flag_list=recvbuf
                if any(match_flag_list):
                    print(f"SHRY, {reference_filename}, the symmetrically-inequivalent test failed.")
                    df.at[row, "Note"] = df.at[row, "Note"] + f"the symmetrically-inequivalent test: SHRY failed."
                    #raise RuntimeError(f"rank={rank}: Failure")
                else:
                    print(f"SHRY, {reference_filename}, the symmetrically-inequivalent test succeeded.")
            
            comm.Barrier()

            # Finally, update df
            if rank==root:
                df.to_excel(f"fullcheck_{timestamp}.xls", index=False)
            else:
                df=None
            #mpi bcast
            df = comm.bcast(df, root=root)
        
        
        # consistency check between SHRY and supercell
        if check_consistency:
            
            # Consistency check between SHRY and Supercell!!:
            if (rank == root): print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            if (rank == root): print("SHRY v.s. Supercell: Check if all the generated structures are consistent between SHRY and Supercell!")
            shry_match_pair_list=[]
            supercell_match_pair_list=[]
            non_match_shry_index_list=[]
            
            supercell_index_list=[i for i in range(len(supercell_cifs))]
            shry_index_list=[i for i in range(len(shry_cifs))]
            
            split_shry_index_list=np.array_split(shry_index_list, size)[rank]
            if (rank == root): print(f"MPI: Num. Core = {size}")
            if (rank == root): print(f"Total num:{len(shry_index_list)}.",flush=True) 
            #for shry_index in tqdm.tqdm(
            #    shry_index_list,
            #    desc="SHRY v.s. Supercell",
            #    **TQDM_CONF,
            #    disable=DISABLE_PROGRESSBAR,
            #):
            supercell_match_index=-1
            
            #check if mpi gather is possible inside the for loop!
            if size > len(shry_index_list):
                split_inside_loop=False
            else:
                split_inside_loop=True
            
            # arrange split_shry_index_list 
            # (i.e. the lengths should be the same for all the loop for the gather command):
            max_length=np.max([len(l) for l in np.array_split(shry_index_list, size)])
            if len(split_shry_index_list) < max_length:
                for _ in range(max_length-len(split_shry_index_list)):
                    split_shry_index_list=np.append(split_shry_index_list,-1)

            for i, shry_index in enumerate(split_shry_index_list):

                # mpiabort test!
                #if rank == 1:
                #    print("rank1 raise error!!")
                #    raise AssertionError(f"rank={rank}: Failure")
                #if (i%10000==0): print(f"rank={rank}, {i}/{len(split_shry_index_list)}.",flush=True)
                match_flag=False
                #print(f"rank={rank}, supercell_match_index {supercell_match_index}")
                #print(f"rank={rank}, supercell_index_list {supercell_index_list}")
                #print(supercell_index_list)
                #print(f"len(supercell_index_list) = {len(supercell_index_list)}")
                
                if shry_index != -1:
                
                    for supercell_index in supercell_index_list:
                        
                        for np_round in np_round_list:
                            try:
                                match_flag=check_structure_consistency(
                                shry_cifs[shry_index], 
                                supercell_cifs[supercell_index], 
                                symmops, 
                                False,
                                np_round)
                                break
                            except AssertionError:
                                match_flag=False
                                continue
                        
                        #print(match_flag)
                        if match_flag:
                            if not split_inside_loop:
                                supercell_index_list.remove(supercell_index)
                            else:
                                supercell_match_index=supercell_index
                            shry_match_pair_list.append(shry_index)
                            supercell_match_pair_list.append(supercell_index)
                            break
                    
                    if not match_flag:
                        #print("For SHRY output: {shry_cifs[shry_index]}")
                        #print("There is no equivalent structure in the supercell outputs!!")
                        non_match_shry_index_list.append(shry_index)
                        supercell_match_index=-1
                
                else:
                    supercell_match_index=-1
                
                if split_inside_loop:
                    #print(f"rank={rank}, supercell_match_index={supercell_match_index}",flush=True)
                    #finally, remove matched supercell for all rank!!
                    #print(f":rank={rank}, supercell_match_index={supercell_match_index}",flush=True)
                    supercell_match_index = comm.gather(supercell_match_index, root=root)
                    #print(f":rank={rank} gather")
                    
                    if rank==root:
                        #print(f"rank={rank}, supercell_match_index={supercell_match_index}",flush=True)
                        for supercell_index in supercell_match_index:
                            if supercell_index != -1:
                                supercell_index_list.remove(supercell_index)
                    else:
                        supercell_index_list=None
                
                    #reset supercell_match_index for all rank.
                    supercell_match_index=-1
                    
                    #mpi bcast
                    if rank==root:
                        pass
                        #print(f"rank={rank}, len(supercell_index_list)={len(supercell_index_list)}",flush=True)
                        #print(f"rank={rank}, supercell_index_list={supercell_index_list}",flush=True)
                    else:
                        pass
                        #print(f"rank={rank}, supercell_index_list={supercell_index_list}",flush=True)
                    supercell_index_list = comm.bcast(supercell_index_list, root=root)
                    #print(f"rank={rank}, supercell_index_list={supercell_index_list}",flush=True)
                    #print(f"rank={rank}, len(supercell_index_list)={len(supercell_index_list)}",flush=True)
                    
                    #comm.Barrier()
                    
                else:
                    pass # done nothing!!!
                
            # MPI supercell_index_list/supercell_match_pair_list
            sendbuf_shry=np.array(shry_match_pair_list, dtype=int)
            sendcounts_shry = np.array(comm.gather(len(sendbuf_shry), root))
            sendbuf_supercell=np.array(supercell_match_pair_list, dtype=int)
            sendcounts_supercell = np.array(comm.gather(len(sendbuf_supercell), root))

            if rank == root:
                recvbuf_shry = np.empty(sum(sendcounts_shry), dtype=int)
                recvbuf_supercell = np.empty(sum(sendcounts_supercell), dtype=int)
            else:
                recvbuf_shry = None
                recvbuf_supercell = None
            
            comm.Gatherv(sendbuf=sendbuf_shry, recvbuf=(recvbuf_shry, sendcounts_shry), root=root)
            comm.Gatherv(sendbuf=sendbuf_supercell, recvbuf=(recvbuf_supercell, sendcounts_supercell), root=root)
            
            if rank==root:
                if len(recvbuf_shry) != len(recvbuf_supercell):
                    print("recvbuf_shry/recvbuf_supercell does not match")
                    raise AssertionError(f"rank={rank}: Failure")
                
                pairs_pd=pd.DataFrame(np.array([[shry_cifs[i] for i in recvbuf_shry], [supercell_cifs[i] for i in recvbuf_supercell]]).T, columns=["shry_cif","supercell_cif"])
                pairs_pd.to_csv(os.path.join(reference_dirname, f"{reference_basename}_pairs.csv"))
            
            # MPI Gather non_match_shry_index_list
            sendbuf=np.array(non_match_shry_index_list)
            sendcounts = np.array(comm.gather(len(sendbuf), root))

            if rank == root:
                recvbuf = np.empty(sum(sendcounts), dtype=int)
            else:
                recvbuf = None
            
            comm.Gatherv(sendbuf=sendbuf, recvbuf=(recvbuf, sendcounts), root=root)

            # write to df
            if rank==root:
                non_match_shry_index_list=recvbuf
                if len(non_match_shry_index_list) > 0:
                    df.at[row, "Note"] = df.at[row, "Note"] + f"No equiv. str in supercell:{[shry_cifs[i] for i in non_match_shry_index_list]}"
                    print("the consistency test was failure. See the xls file.")
                    #raise RuntimeError(f"rank={rank}: Failure")
                
            comm.Barrier()
            if (rank == root): print("SHRY v.s. Supercell: the consistency test has done.") 
        # fi
        
        # Finally, if all the test passes. / MPI=root
        #print(df)
        if rank==root:
            df.at[row, "Checked"] = True
            df.to_excel(f"fullcheck_{timestamp}.xls", index=False)
        else:
            df=None
        #mpi bcast
        df = comm.bcast(df, root=root)
        
        # Also, MPI=root, delete extracted files
        if rank==root:
            #finally, delete the extracted dirs
            #print("Delete zip-extracted dirs")
            for zip_file in zip_files:
                filebase = os.path.basename(zip_file)
                filebasenozip = filebase.replace(".zip", "")
                shutil.rmtree(os.path.join(reference_dirname,filebasenozip))

    except Exception as e:
        print(f"rank={rank}, MPI Abort call.")
        print(f"type:{str(type(e))}")
        print(e.args)
        df.to_excel(f"fullcheck_{timestamp}_rank_{rank}.xls", index=False)
        comm.Abort() # abort all MPI!!

if rank==root: 
    print("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print(f"Test Done: finalizing MPI...")

#MPI finalize
MPI.Finalize()
sys.exit(0)
