import warnings
warnings.simplefilter('ignore')

import os
os.environ["OMP_NUM_THREADS"] = "1" # setenv OMP_NUM_THREADS 1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # setenv OPENBLAS_NUM_THREADS 1
os.environ["MKL_NUM_THREADS"] = "1" # setenv MKL_NUM_THREADS 1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # setenv VECLIB_MAXIMUM_THREADS 1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # setenv NUMEXPR_NUM_THREADS 1

import math
import time
import subprocess
import sys
from datetime import datetime
from tqdm import tqdm

import pandas as pd
import numpy as np

TIMECOMMAND=["/usr/bin/time", "-p"]
LS = ["ls"]
REAL_TIME_COLUMN = 0
USER_TIME_COLUMN = 1
SYS_TIME_COLUMN =  2
TIME_LOOP=3

SUPERCELL_TOLERANCE="0.01" #angstrom
SHRY_TOLERANCE="0.01"      #angstrom
SHRY_ANGLE_TOLERANCE="5.0"

def timed_process_run(args, desc=None, average=TIME_LOOP):
    """
    Throws:
        - CalledProcessError: anything with the called program
        - AssertionError: weird time line
    """
    args = TIMECOMMAND + args
    r_list=[]; u_list=[]; s_list=[]
    
    for i in tqdm(range(average), desc=desc):
        try:
            p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            stdout_lines = p.stdout.decode().split("\n")
            stderr_lines = p.stderr.decode().split("\n")
        
            def has_real_time_output(string):
                return "real" in string
        
            def has_user_time_output(string):
                return "user" in string
        
            def has_sys_time_output(string):
                return "sys" in string
        
            time_line = [x for x in stderr_lines 
            if has_real_time_output(x) or
            has_user_time_output(x) or
            has_sys_time_output(x)
            ]
            assert len(time_line) == 3
        
            r = float(time_line[REAL_TIME_COLUMN].split()[1])
            u = float(time_line[USER_TIME_COLUMN].split()[1])
            s = float(time_line[SYS_TIME_COLUMN].split()[1])
            
            r_list.append(r); u_list.append(u); s_list.append(s)
        
        except:
            print(f"Subprocess error: time trial {i}/{len(range(average))}")
    
    try:
        assert len(r_list) != 0
        assert len(u_list) != 0
        assert len(s_list) != 0
        
        r=np.mean(r_list); u=np.mean(u_list); s=np.mean(s_list);
        return r, u, s, stdout_lines
            
    except AssertionError as e:
        import traceback
        print(traceback.format_exc())
        #print(e)
        raise AssertionError

# Supercell run
def has_supercell_substitutions(string):
    return "The total number of combinations is" in string

def has_supercell_equivalent(string):
    return "Combinations after merge" in string

def get_supercell_number(string):
    return int(string.split()[-1].split("(")[0])

# Shry run
def has_shry_equivalent(string):
    return "Expected total of" in string

def get_shry_number(string):
    return int(string.split()[3])
    
def main():
    
    filename="9011358_partial.cif"
    supercell_trial_list=["1x1x1", "1x1x2", "1x1x3", "1x2x2", "2x2x2"]
    
    for supercell in supercell_trial_list:
        supercell_list = supercell.split("x")
        supercell_args = ["supercell", "-m", "-d", "-s", supercell, "-i", filename, "-t", SUPERCELL_TOLERANCE]
        shry_args = ["shry","-s",supercell_list[0],supercell_list[1],supercell_list[2],"--no-write",
                    filename, "--symprec", SHRY_TOLERANCE, "--angle-tolerance", SHRY_ANGLE_TOLERANCE
                    ]
    
        r, u, s, stdout_lines = timed_process_run(supercell_args, desc="[time supercell]")
    
        supercell_substitutions = [
            x for x in stdout_lines if has_supercell_substitutions(x)
        ]
        assert len(supercell_substitutions) == 1
        supercell_substitutions = get_supercell_number(supercell_substitutions[0])
        
        supercell_equivalent = [
            x for x in stdout_lines if has_supercell_equivalent(x)
        ]
        assert len(supercell_equivalent) == 1
        # supercell_equivalent = int(supercell_equivalent[0].split()[-1])
        supercell_equivalent = get_supercell_number(supercell_equivalent[0])
    
        print("------------------------------------------------------")
        print(f"supercell = {supercell}")
        print("T_supercell: real, user, sys, {:.2f}, {:.2f}, {:.2f}".format(r,u,s))
        print(f"supercell_substitutions={supercell_substitutions}")
        print(f"supercell_equivalent={supercell_equivalent}")
        print("------------------------------------------------------")
        
        """
        r, u, s, stdout_lines = timed_process_run(shry_args, desc="[time shry]")
        
        shry_equivalent = [
            x for x in stdout_lines if has_shry_equivalent(x)
        ]
        assert len(shry_equivalent) == 1
        shry_equivalent = get_shry_number(shry_equivalent[0])
        
        print("------------------------------------------------------")
        print(f"supercell = {supercell}")
        print("T_shry: real, user, sys, {:.2f}, {:.2f}, {:.2f}".format(r,u,s))
        print(f"shry_equivalent={shry_equivalent}")
        print("------------------------------------------------------")
        """
        
if __name__ == "__main__":
    main()
