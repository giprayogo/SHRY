#!/bin/bash

for no in `seq 5 230`
    do
    if [ -e SG${no} ]; then

        file_list=`ls SG${no} | grep "SG${no}.*\d+*\.cif"`
    
        for file in $file_list
            do
            prefix=`echo $file | rev | cut -c 5- | rev`
            
            cif_org=${prefix}.cif
            cif_sub=${prefix}_partial.cif
            echo ==$cif_org==
            
            HM=`grep "_space_group_name_H-M" SG${no}/$cif_org | tail -n 1 | awk -F "\'" '{print $2}'`
            echo HM=$HM
            HM_line="_symmetry_space_group_name_H-M  \'${HM}\'"
            
            line_num=`grep -n "_space_group_name_H-M" SG${no}/$cif_sub | awk -F ':' '{print $1}'`
            echo $line_num
            cp SG${no}/$cif_sub  SG${no}/bak_${cif_sub}
            #gsed -i "${line_num}d" SG${no}/bak_${cif_sub}
            #gsed -i -e "${line_num}i ${HM_line}" SG${no}/bak_${cif_sub}
            gsed -i "${line_num}d" SG${no}/${cif_sub}
            gsed -i -e "${line_num}i ${HM_line}" SG${no}/${cif_sub}
            echo ""
            done
    fi
done
