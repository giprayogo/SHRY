import warnings
warnings.simplefilter('ignore')

import sys
import pandas as pd

def main():
    df_all = pd.read_excel("benchmark_SG_all.xls", index_col=0)
    
    sg_list=set(list(df_all["SpaceGroup_No"]))
    
    for sg in sg_list:
        df_all[df_all["SpaceGroup_No"] == sg].to_excel(f"benchmark_SG_{sg}.xls")
    
if __name__ == "__main__":
    main()
