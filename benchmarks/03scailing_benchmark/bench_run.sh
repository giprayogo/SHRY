#!/bin/tcsh

#set start and end nos.
set start_sg=1
set end_sg=10

#bench
foreach sg ( `seq $start_sg $end_sg`)
    if ( -e "SG${sg}" ) then
        #on a local machine
        python bench.py benchmark_SG_${sg}.xls SG_${sg}
        
        #on a PBS system
        #qsub bench_sg$sg.sh
    endif
end

