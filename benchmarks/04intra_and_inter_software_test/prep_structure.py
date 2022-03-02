import warnings
warnings.simplefilter('ignore')

import math
import os
import shutil
import subprocess
import sys
from datetime import datetime
import glob

import numpy as np
import pandas as pd

SUPERCELL_TOLERANCE="1.0e-2" #angstrom
SHRY_TOLERANCE="1.0e-2"      #angstrom
SHRY_ANGLE_TOLERANCE="5.0"   #degree

BABEL=False # True for supercell ver.1.2, False for that >= ver.2.0
SUPERCELL_PATH="/home/nkousuke/application/supercell/ver.2.0.2/supercell"
SHRY_PATH = "/home/nkousuke/application/shry/ver.1.1.0"

SUPERCELL_FLAG=True
SHRY_FLAG=True
SHRY_NOCASH_FLAG=False
SHRY_DMAT_FLAG=False

def process_run(args):
    
    my_env = os.environ.copy()
    if BABEL: my_env["BABEL_DATADIR"] = SUPERCELL_PATH
    my_env["PYTHONPATH"] = SHRY_PATH
    
    p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, env=my_env)
    stdout_lines = p.stdout.decode().split("\n")
    stderr_lines = p.stderr.decode().split("\n")
    
    return stdout_lines, stderr_lines

def main():
    timestamp = int(datetime.now().timestamp())
    df = pd.read_excel(sys.argv[1])
    rows, columns = df.shape
    #df["Checked"] = df["Checked"].astype(float)
    df["Checked"] = False
    df["Note"] = df["Note"].astype(str)
    
    success_row_list=[]
    
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
            #assert isinstance(checked, float)
            # Skip checked structures.
            if checked:
                continue

            filename = filename.replace(".cif", "_partial.cif")
            filedir = os.path.abspath(os.path.dirname(filename))
            filebase = os.path.basename(filename)
            filebasenocif = filebase.replace(".cif", "")
            supdir = f"supercell-{filebasenocif}"
            rootpath = os.getcwd()
            print(f"Working on {filename} ({row+1}/{rows})")

            supercell_list = supercell.split("x")
            supercell_args = [SUPERCELL_PATH, "-m", "-s", supercell, "-i", filebase, "-t", SUPERCELL_TOLERANCE]
            shry_args = [
                "python",
                "-m",
                "shry.script",
                filebase,
                "-s",
                supercell_list[0],
                supercell_list[1],
                supercell_list[2],
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
            if SUPERCELL_FLAG:
                try:
                    os.chdir(filedir)
                    os.makedirs(supdir, exist_ok=True)
                    shutil.copy(filebase, supdir)
                    os.chdir(supdir)
                    stdout_lines, stderr_lines = process_run(supercell_args)
                    
                    #archive!!
                    os.chdir(filedir)
                    shutil.make_archive(supdir, 'zip', root_dir=supdir)
                    shutil.rmtree(supdir)
    
                except Exception as e:
                    print(f"type:{str(type(e))}")
                    print(e.args)
    
                finally:
                    os.chdir(rootpath)
            else:
                print("Supercell skipped.")
            
            if SHRY_FLAG:
                # SHRY run
                try:
                    os.chdir(filedir)
                    stdout_lines, stderr_lines = process_run(shry_args)
                    
                    #archive!!
                    shry_glob = "shry-*{}*".format(filebasenocif)
                    shry_dirs = [shry_dir for shry_dir in glob.glob(shry_glob) if os.path.isdir(shry_dir)]
                    for shry_dir in shry_dirs:
                        shutil.make_archive(shry_dir, 'zip', root_dir=shry_dir)
                        shutil.rmtree(shry_dir)
                    
                except Exception as e:
                    print(f"type:{str(type(e))}")
                    print(e.args)
                finally:
                    os.chdir(rootpath)
            else:
                print("SHRY skipped.")
                
            print("============================================")
            
            # finally if both shry and supercell are successful,
            df.at[row, "Checked"] = True
            df.at[row, "Note"] = df.at[row, "Note"] + ", prep-ok"
            df.to_excel(f"out_prepstruct.xls", index=False)
    
    except Exception as e:
        print(f"type:{str(type(e))}")
        print(e.args)
    
    finally:
        print(df)

if __name__ == "__main__":
    main()
