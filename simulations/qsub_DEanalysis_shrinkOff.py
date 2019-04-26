# coding: utf-8

# In[ ]:


import os


# In[ ]:


fc_vec = [1, 1.5, 2, 3]
size_vec = [1, 10, 100, 1000]
Nsample_vec = [10, 20, 40]
design_vec = ['TRUE', 'FALSE']
coverage_vec = [1, 5, 10, 20]


# In[ ]:


### Fold changes (FC)
for i in range(len(fc_vec)):
    os.system("qsub -P combat -o logs/FC" + str(i+1) + ".txt -e logs/err_FC" + str(i+1) + ".txt " +  \
              "-N reFC" + str(i+1) + " -cwd -b y -pe omp 8 " +  \
              "/share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe_polyester.R " +  \
              "FC 1.3 " + str(fc_vec[i]) + " 100 10 20 TRUE 5")


# In[ ]:


## Dispersion differences (Disp)
for i in range(len(size_vec)):
    os.system("qsub -P combat -o logs/Disp" + str(i+1) + ".txt -e logs/err_Disp" + str(i+1) + ".txt " +  \
              "-N reDisp" + str(i+1) + " -cwd -b y -pe omp 8 " + \
              "/share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe_polyester.R " + \
              "Disp 1.3 1.5 " + str(size_vec[i]) + " 10 20 TRUE 5")


# In[ ]:


## Study design (Design)
for i in range(len(design_vec)):
    os.system("qsub -P combat -o logs/Design" + str(i+1) + ".txt -e logs/err_Design" + str(i+1) + ".txt " + \
              "-N reDesign" + str(i+1) + " -cwd -b y -pe omp 8 " + \
              "/share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe_polyester.R " +  \
              "Design 1.3 1.5 100 10 20 " + design_vec[i] + " 5")


# In[ ]:


## Number of samples (Nsample)
for i in range(len(Nsample_vec)):
    os.system("qsub -P combat -o logs/Nsample" + str(i+1) + ".txt -e logs/err_Nsample" + str(i+1) + ".txt " + \
              "-N reNsample" + str(i+1) + " -cwd -b y -pe omp 8 " +  \
              "/share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe_polyester.R " +  \
              "Nsample 1.3 1.5 100 10 " + str(Nsample_vec[i]) + " TRUE 5")


# In[ ]:


## Reads per gene (Depth)
for i in range(len(coverage_vec)):
    os.system("qsub -P combat -o logs/Depth" + str(i+1) + ".txt -e logs/err_Depth" + str(i+1) + ".txt " +  \
              "-N reDepth" + str(i+1) + " -cwd -b y -pe omp 8 " +  \
              "/share/pkg/r/3.5.2/install/bin/Rscript sim_DEpipe_polyester.R " +  \
              "Depth 1.3 1.5 100 10 20 TRUE " + str(coverage_vec[i]))


# In[ ]:




