#!/bin/tcsh

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

#exit

set start_sg=11
set end_sg=230

#bench
foreach sg ( `seq $start_sg $end_sg`)
    if ( -e "SG${sg}" ) then
	sed 's/sg_input/'$sg'/g' bench.sh > bench_sg$sg.sh
	qsub bench_sg$sg.sh
	rm bench_sg$sg.sh
    endif
end

