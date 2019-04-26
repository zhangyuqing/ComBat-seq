# coding: utf-8
import os
import numpy as np
#alpha_fdr_seq = ["%.2f" % i for i in np.arange(0.01, 0.16, 0.01)]
alpha_fdr_seq = ["%.2f" % i for i in np.array([0.01, 0.06])]
for i in range(len(alpha_fdr_seq)):
	os.system("qsub -P combat -l h_rt=72:00:00 -o ./exprep_results/alpha" + str(i+1) + \
			  ".txt -e ./exprep_results/err_alpha" + str(i+1) + ".txt " +  \
              "-N ExpRep" + str(i+1) + " -cwd -b y -pe omp 8 " +  \
              "/share/pkg/r/3.5.2/install/bin/Rscript ../exprep_pipe.R " + alpha_fdr_seq[i])

