# coding: utf-8

import os
disp_seq = [1, 2, 3, 4, 5]
cnfnd_seq = [0.5, 0.4, 0.3, 0.2]

for i in range(len(disp_seq)):
    for j in range(len(cnfnd_seq)):
    	curr_exp_name = "EBout_Disp" + str(disp_seq[i]) + "_Cnfnd" + str(int(cnfnd_seq[j]*10))
    	os.system("qsub -P combat -l h_rt=96:00:00 -o ./DE_libcomp_EBcomplete_outliers/logs/" + curr_exp_name + ".txt" + \
    		      " -e ./DE_libcomp_EBcomplete_outliers/logs/err_" + curr_exp_name + ".txt" + " -N " + curr_exp_name + \
    		      " -cwd -b y -pe omp 8 /share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe_compFDR_EBcomplete_outliers.R " + \
    		      str(disp_seq[i]) + " " + str(cnfnd_seq[j]) + " 20 5")
    	




