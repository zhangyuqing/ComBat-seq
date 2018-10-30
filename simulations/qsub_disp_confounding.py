#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.system("module load R/3.5.0")
os.system("module load gcc/7.2.0")



# In[2]:


disp_levels = [1, 2, 3, 4, 5]
confounding_levels = [0.5, 0.4, 0.3, 0.2, 0.1]


# In[5]:


for i in range(len(disp_levels)):
    for j in range(len(confounding_levels)):
        os.system("qsub -P combat -o logs/Disp" + str(disp_levels[i]) + "_Cnfnd" + str(int(confounding_levels[j]*10)) + ".txt" + \
        	" -e logs/err_Disp" + str(disp_levels[i]) + "_Cnfnd" + str(int(confounding_levels[j]*10)) + ".txt" +  \
        	" -N Disp" + str(disp_levels[i]) + "_Cnfnd" + str(int(confounding_levels[j]*10)) + \
        	" -cwd -b y -pe omp 8 Rscript sim_DEpipe_disp_confounding.R " + str(disp_levels[i]) + " " + str(confounding_levels[j]))


# In[ ]:




