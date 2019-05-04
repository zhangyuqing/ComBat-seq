# About ComBat-seq

<!--- 
[![Travis build status](https://travis-ci.com/zhangyuqing/ComBat-seq.svg?branch=master)](https://travis-ci.com/zhangyuqing/ComBat-seq)
--->

ComBat-seq is a batch effect adjustment (BEA) tool for RNA-seq read counts, using Negative Binomial regression. To use it, download the ComBat_seq.R and helper_seq.R scripts, and run the function with the syntax below:

```r
source("ComBat_seq.R")
source("helper_seq.R")
# include condition (group variable)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group, full_mod=TRUE)
# do not include condition
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=NULL, full_mod=FALSE)
```


