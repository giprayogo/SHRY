#!/bin/tcsh

#activate shry
setenv PATH "/home/nkousuke/application/anaconda3/bin:$PATH"
conda activate shry

set start_sg=1
set end_sg=230

#autosub
foreach sg ( `seq $start_sg $end_sg`)
	sed 's/sg_input/'$sg'/g' autosub.sh > autosub_sg$sg.sh
	qsub autosub_sg$sg.sh
	rm autosub_sg$sg.sh
end

