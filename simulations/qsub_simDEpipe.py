# coding: utf-8

import os

mean_seq = [1, 1.5, 2, 3]
disp_seq = [1, 2, 3, 4]

for i in range(len(mean_seq)):
	for j in range(len(disp_seq)):
		curr_exp_name = "largeM" + str(int(mean_seq[i]*10)) + "D" + str(disp_seq[j])
		os.system("qsub -P combat -l h_rt=160:00:00 -o ./DE_large/logs/" + curr_exp_name + ".txt" + \
			" -e ./DE_large/logs/err_" + curr_exp_name + ".txt" + " -N " + curr_exp_name + \
			" -cwd -b y -pe omp 8 /share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe.R " + \
			str(mean_seq[i]) + " " + str(disp_seq[j]) + " 20")
    	
 



