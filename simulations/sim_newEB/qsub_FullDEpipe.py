# coding: utf-8

import os

mean_seq = [1.5] #[1, 1.5, 2, 3]
disp_seq = [1, 2, 3, 4, 5]
cnfnd_seq = [0.5, 0.4, 0.3, 0.2]

for i in range(len(mean_seq)):
	for j in range(len(disp_seq)):
		for k in range(len(cnfnd_seq)):
			curr_exp_name = "M" + str(int(mean_seq[i]*10)) + "D" + str(disp_seq[j]) + "C" + str(int(cnfnd_seq[k]*10))
			os.system("qsub -P combat -l h_rt=96:00:00 -o ./DE_fullpipe/logs/" + curr_exp_name + ".txt" + \
    		          " -e ./DE_fullpipe/logs/err_" + curr_exp_name + ".txt" + " -N " + curr_exp_name + \
    		          " -cwd -b y -pe omp 8 /share/pkg/r/3.5.2/install/bin/Rscript sim_FullDEpipe.R " + \
    		          str(mean_seq[i]) + " " + str(disp_seq[j]) + " " + str(cnfnd_seq[k]) + " 20 5")
    	
#-l mem_per_core=16G 



