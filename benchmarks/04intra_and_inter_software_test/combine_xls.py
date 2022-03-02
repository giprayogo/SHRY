import warnings
warnings.simplefilter('ignore')

import shutil
import os
import time
import subprocess
import sys
import pandas as pd

def main():
    df_list=[]
    sg_list=[f"SG{i}" for i in range(1,231)]
    
    for sg in sg_list:
        if not os.path.exists(f"autosub_{sg}.xls"):
            continue
        df_list.append(pd.read_excel(f"autosub_{sg}.xls", index_col=0))

    for i, df_pd in enumerate(df_list):
        if i == 0:
            df_all= pd.DataFrame(columns=df_pd.columns)
            for j in range(len(df_pd)):
                df_all=df_all.append(df_pd.loc[j], ignore_index=True)
        else:
            for j in range(len(df_pd)):
                df_all=df_all.append(df_pd.loc[j], ignore_index=True)
    
    #df_all=df_all.reset_index(drop=True)
    df_all.to_excel(f"autosub_SG_all.xls")
    
if __name__ == "__main__":
    main()
