import warnings
warnings.simplefilter('ignore')

import os
import math
import time
import subprocess
import sys
from datetime import datetime
from tqdm import tqdm

import pandas as pd
import numpy as np

# NOTE: $GNU time command!! Not built-in.
# NOTE: units are seconds(time) and Kibibytes=KiB(memory).
TIMECOMMAND_WO_MEM=["/usr/bin/time", "-p"]
TIMECOMMAND_W_MEM=["/usr/bin/time", "-f", "real %e\nuser %U\nsys %S\nMemory %M"]
LS = ["ls"]
REAL_TIME_INDEX = 0
USER_TIME_INDEX = 1
SYS_TIME_INDEX =  2
MEMORY_INDEX =  3
TIME_LOOP=2

SUPERCELL_TOLERANCE="0.01" #angstrom
SHRY_TOLERANCE="0.01"      #angstrom
SHRY_ANGLE_TOLERANCE="5.0"

BABEL=False # True for supercell ver.1.2, False for that >= ver.2.0
SUPERCELL_PATH="/home/nkousuke/application/supercell/ver.2.0.2/supercell"
SHRY_PATH = "/home/nkousuke/application/shry/ver.1.1.0"

STDOUT_DIR="stdout"
STDERR_DIR="stderr"

SUPERCELL_FLAG=True
SHRY_FLAG=False
SHRY_NOCASH_FLAG=False
SHRY_DMAT_FLAG=False

def timed_process_run(args, desc=None, average=TIME_LOOP):
    """
    Throws:
        - CalledProcessError: anything with the called program
        - AssertionError: weird time line
    """
    args = TIMECOMMAND_W_MEM + args
    r_list=[]; u_list=[]; s_list=[]; m_list=[]
    
    for i in tqdm(range(average), desc=desc):
        try:
            #print(args)
            my_env = os.environ.copy()
            if BABEL: my_env["BABEL_DATADIR"] = SUPERCELL_PATH
            my_env["PYTHONPATH"] = SHRY_PATH
            p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, env=my_env)
            stdout_lines = p.stdout.decode().split("\n")
            stderr_lines = p.stderr.decode().split("\n")
            
            #print(p.stdout.decode())
            #print(p.stderr.decode())

            os.makedirs(STDOUT_DIR, exist_ok=True)
            os.makedirs(STDERR_DIR, exist_ok=True)

            with open(os.path.join(STDOUT_DIR, f"stdout_{desc}_{i}.out"), "w") as f:
                f.writelines(p.stdout.decode())
            with open(os.path.join(STDERR_DIR, f"stderr_{desc}_{i}.out"), "w") as f:
                f.writelines(p.stderr.decode())

            def has_real_time_output(string):
                return "real" in string
        
            def has_user_time_output(string):
                return "user" in string
        
            def has_sys_time_output(string):
                return "sys" in string
            
            def has_memory_output(string):
                return "Memory" in string

            time_line = [x for x in stderr_lines 
            if has_real_time_output(x) or
            has_user_time_output(x) or
            has_sys_time_output(x) or
            has_memory_output(x)
            ]
            assert len(time_line) == 4
        
            r = float(time_line[REAL_TIME_INDEX].split()[1])
            u = float(time_line[USER_TIME_INDEX].split()[1])
            s = float(time_line[SYS_TIME_INDEX].split()[1])
            m = float(time_line[MEMORY_INDEX].split()[1])
            
            r_list.append(r); u_list.append(u); s_list.append(s); m_list.append(m)
        
        except:
            print(f"Subprocess error: time trial {i+1}/{average}")
            import traceback
            print(traceback.format_exc())
    try:
        assert len(r_list) != 0
        assert len(u_list) != 0
        assert len(s_list) != 0
        assert len(m_list) != 0
        
        r=np.mean(r_list); u=np.mean(u_list); s=np.mean(s_list); m=np.mean(m_list)
        return r, u, s, m, stdout_lines
            
    except AssertionError as e:
        import traceback
        print(traceback.format_exc())
        #print(e)
        raise AssertionError


def main():
    timestamp = int(datetime.now().timestamp())
    df = pd.read_excel(sys.argv[1])
    prefix = sys.argv[2]
    rows, columns = df.shape
    df["Checked"] = df["Checked"].astype(float)
    df["Note"] = df["Note"].astype(str)
    # Add time columns
    if SUPERCELL_FLAG:
        df["T_supercell_real"] = np.nan
        df["T_supercell_user"] = np.nan
        df["T_supercell_sys"] = np.nan
        df["M_supercell_memory"] = np.nan
    if SHRY_FLAG:
        df["T_shry_real"] = np.nan
        df["T_shry_user"] = np.nan
        df["T_shry_sys"] = np.nan
        df["M_shry_memory"] = np.nan
        
    try:
        for row, zipped in enumerate(
            zip(
                df["Supercell"],
                df["File"],
                df["Substitutions"],
                df["Equivalent Structures"],
                df["Checked"],
            )
        ):
            supercell, filename, substitution, equivalent, checked = zipped
            assert isinstance(checked, float)
            if not math.isnan(checked):
                continue

            cif_basename = os.path.basename(filename).replace(".cif", "")
            filename = filename.replace(".cif", "_partial.cif")
            print(f"Working on {filename} ({row+1}/{rows})")

            substitution = int(substitution)
            equivalent = int(equivalent)
            
            # Supercell settings
            supercell_list = supercell.split("x")
            supercell_args = [SUPERCELL_PATH, "-m", "-d", "-s", supercell, "-i", filename, "-t", SUPERCELL_TOLERANCE]
            
            # SHRY settings
            shry_args = [
                "python",
                "-m",
                "shry.script",
                "-s",
                supercell_list[0],
                supercell_list[1],
                supercell_list[2],
                "--no-write",
                filename,
                "--symprec",
                SHRY_TOLERANCE,
                "--angle-tolerance",
                SHRY_ANGLE_TOLERANCE
            ]
            
            if SHRY_DMAT_FLAG:
                shry_args.append("--no-dmat")
            if SHRY_NOCASH_FLAG:
                shry_args.append("--no-cache")


            # Supercell run
            def has_supercell_substitutions(string):
                return "The total number of combinations is" in string

            def has_supercell_equivalent(string):
                return "Combinations after merge" in string
            
            def get_supercell_number(string):
                return int(string.split()[-1].split("(")[0])
    
            if SUPERCELL_FLAG:
                try:
                    r, u, s, m, stdout_lines = timed_process_run(supercell_args, desc=f"time_supercell_{cif_basename}")
                    supercell_substitutions = [
                        x for x in stdout_lines if has_supercell_substitutions(x)
                    ]
                    assert len(supercell_substitutions) == 1
                    try:
                        # supercell_substitutions = int(supercell_substitutions[0].split()[-1])
                        supercell_substitutions = get_supercell_number(supercell_substitutions[0])
                    except ValueError:
                        # keep as string
                        supercell_substitutions = supercell_substitutions[0].split()[-1]
    
                    supercell_equivalent = [
                        x for x in stdout_lines if has_supercell_equivalent(x)
                    ]
                    assert len(supercell_equivalent) == 1
                    try:
                        # supercell_equivalent = int(supercell_equivalent[0].split()[-1])
                        supercell_equivalent = get_supercell_number(supercell_equivalent[0])
                    except ValueError:
                        # keep as string
                        supercell_equivalent = supercell_equivalent[0].split()[-1]
    
                except (subprocess.CalledProcessError, AssertionError) as e:
                    print(e)
                    r, u, s, m= (0.0, 0.0, 0.0, 0.0)
                    supercell_equivalent = 0
                    supercell_substitutions = 0
                    df.at[row, "Note"] = df.at[row, "Note"] + ", supercell failed"
                df.at[row, "T_supercell_real"] = r
                df.at[row, "T_supercell_user"] = u
                df.at[row, "T_supercell_sys"] = s
                df.at[row, "M_supercell_memory"] = m
    
                try:
                    assert supercell_equivalent == equivalent
                    assert supercell_substitutions == substitution
                except AssertionError as e:
                    import traceback
                    print(traceback.format_exc())
                    #print(e)
                    df.at[row, "Note"] = df.at[row, "Note"] + ", supercell mismatch"
            
            else:
                print("Supercell test skipped.")
                
            time.sleep(1)
            
            # SHRY run
            def has_shry_substitutions(string):
                return "Total number of combinations is" in string 

            def has_shry_equivalent(string):
                return "Expected unique patterns is" in string
            
            def get_shry_number(string):
                return int(string.split()[-1])

            if SHRY_FLAG:
                try:
                    r, u, s, m, stdout_lines = timed_process_run(shry_args, desc=f"time_shry_{cif_basename}")
                    
                    
                    shry_substitutions = [
                        x for x in stdout_lines if has_shry_substitutions(x)
                    ]
                    assert len(shry_substitutions) == 1
                    try:
                        shry_substitutions = get_shry_number(shry_substitutions[0])
                    except ValueError:
                        # keep as string
                        raise ValueError
                        #shry_substitutions = supercell_substitutions[0].split()[-1]
    
                    shry_equivalent = [
                        x for x in stdout_lines if has_shry_equivalent(x)
                    ]
                    assert len(shry_equivalent) == 1
                    try:
                        shry_equivalent = get_shry_number(shry_equivalent[0])
                    except ValueError:
                        # keep as string
                        raise ValueError
                        #shry_equivalent = shry_equivalent[0].split()[-1]
                        
                except (subprocess.CalledProcessError, AssertionError) as e:
                    import traceback
                    print(traceback.format_exc())
                    #print(e)
                    stdout_lines = []
                    r, u, s, m = (0.0, 0.0, 0.0, 0.0)
                    df.at[row, "Note"] = df.at[row, "Note"] + ", shry failed"
                df.at[row, "T_shry_real"] = r
                df.at[row, "T_shry_user"] = u
                df.at[row, "T_shry_sys"] = s
                df.at[row, "M_shry_memory"] = m
    
                time.sleep(1)
            
            else:
                print("SHRY test skipped.")
                
            if SUPERCELL_FLAG:
                print(
                    f"Supercell/XSL file => "
                    f"Combinations:{supercell_equivalent}/{equivalent}; "
                    f"Substitutions:{supercell_substitutions}/{substitution}"
                )
                print("")
            
            if SHRY_FLAG:
                print(
                    f"SHRY/XSL file => "
                    f"Combinations:{shry_equivalent}/{equivalent}; "
                    f"Substitutions:{shry_substitutions}/{substitution}"
                )
                print("")
            
            print("============================================")
    finally:
        
        # write to excel
        df.to_excel(f"out_benchmark_{prefix}.xls", index=False)

if __name__ == "__main__":
    main()
