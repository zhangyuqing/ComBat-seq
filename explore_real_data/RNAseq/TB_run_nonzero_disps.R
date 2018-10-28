####  The dispersion estimates using non-zero portion of genes take too...ooo long to run. 
####  To save time, I'll run it once and save it as RData to avoid waiting too long each time I need to recreate the report
rm(list=ls())
setwd("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB")
sapply(c("SummarizedExperiment", "DESeq2", "edgeR", "dendextend", "ggplot2", "reshape2", "gridExtra", "scales", 
         "ggdendro", "MASS"), require, character.only=TRUE)

####  Load data
rds_obj <- readRDS("combined.rds")
cts <- assays(rds_obj)$counts
batch <- colData(rds_obj)$SequencingBatch
group <- colData(rds_obj)$Label
glm_ctrl <- glm.control(maxit=200)


####  Define function to compute dispersion
# calculate dispersion for one gene
genDisp_1gene <- function(curr_cts, curr_group, glm_ctrl){
  if(length(unique(curr_cts)) > 2 & var(curr_cts) > 0){
    if(length(unique(curr_group))<2){
      curr_f <- try(glm.nb(curr_cts ~ 1, control=glm_ctrl), silent=TRUE)
    }else{
      curr_f <- try(glm.nb(curr_cts ~ as.factor(curr_group), control=glm_ctrl), silent=TRUE)
    }
    if(class(curr_f)[1]!="try-error"){
      return(1/summary(curr_f)$theta)
    }
  }
  return(NA)
}

# use counts containing zeros
genDisp <- function(cts, group, batch, glm_ctrl, batch_name){
  disp_full_batch <- sapply(1:nrow(cts), function(ind){
    genDisp_1gene(curr_cts=cts[ind, batch==batch_name], curr_group=group[batch==batch_name], glm_ctrl=glm_ctrl)
  })
  disp_full_batch <- disp_full_batch[!is.na(disp_full_batch)]
  return(disp_full_batch)
}

# use counts without zeros
genDisp_Sep <- function(nonzero_cts, batch_sep, group_sep, glm_ctrl, batch_name){
  # take the batch specified by batch_name
  group_sep_thisbatch <- nonzero_cts_thisbatch <- list()
  for(i in 1:length(nonzero_cts)){
    group_sep_thisbatch[[i]] <- group_sep[[i]][batch_sep[[i]]==batch_name]
    nonzero_cts_thisbatch[[i]] <- nonzero_cts[[i]][batch_sep[[i]]==batch_name]
  }
  names(nonzero_cts_thisbatch) <- names(group_sep_thisbatch) <- names(nonzero_cts)
  
  # calculate dispersions
  disp_nonzero_thisbatch <- sapply(1:length(nonzero_cts_thisbatch), function(ii){
    curr_cts <- nonzero_cts_thisbatch[[ii]] 
    curr_group <- group_sep_thisbatch[[ii]]
    return(genDisp_1gene(curr_cts, curr_group, glm_ctrl))
  })
  disp_nonzero_thisbatch <- disp_nonzero_thisbatch[!is.na(disp_nonzero_thisbatch)]
  return(disp_nonzero_thisbatch)
}



####  Use counts WITH zeros to estimate dispersion
disp_full_africa <- genDisp(cts=cts, group=group, batch=batch, glm_ctrl=glm_ctrl, batch_name="Africa")
disp_full_G6 <- genDisp(cts=cts, group=group, batch=batch, glm_ctrl=glm_ctrl, batch_name="G6")
disp_full_india <- genDisp(cts=cts, group=group, batch=batch, glm_ctrl=glm_ctrl, batch_name="India")
disp_full_bzl1 <- genDisp(cts=cts, group=group, batch=batch, glm_ctrl=glm_ctrl, batch_name="Brazil_1")
disp_full_bzl2 <- genDisp(cts=cts, group=group, batch=batch, glm_ctrl=glm_ctrl, batch_name="Brazil_2")



####  Use only non-zero portions of each gene to estimate dispersion
##  Take the non-zero portion in each gene
batch_sep <- group_sep <- nonzero_cts <- list()
for(i in 1:nrow(cts)){
  nonzero_ind <- which(cts[i, ]!=0)
  batch_sep[[i]] <- as.character(batch[nonzero_ind])
  group_sep[[i]] <- as.character(group[nonzero_ind])
  nonzero_cts[[i]] <- cts[i, nonzero_ind]
}
# sanity check
if(!identical(sapply(nonzero_cts, length), as.integer(rowSums(cts!=0)))){stop("Error in taking non-zero portions!")}
names(nonzero_cts) <- rownames(cts); names(batch_sep) <- rownames(cts); names(group_sep) <- rownames(cts)

## Calculate dispersion in each batch
disp_nonzero_africa <- genDisp_Sep(nonzero_cts=nonzero_cts, batch_sep=batch_sep, group_sep=group_sep, glm_ctrl=glm_ctrl, batch_name="Africa")
disp_nonzero_G6 <- genDisp_Sep(nonzero_cts=nonzero_cts, batch_sep=batch_sep, group_sep=group_sep, glm_ctrl=glm_ctrl, batch_name="G6")
disp_nonzero_india <- genDisp_Sep(nonzero_cts=nonzero_cts, batch_sep=batch_sep, group_sep=group_sep, glm_ctrl=glm_ctrl, batch_name="India")
disp_nonzero_bzl1 <- genDisp_Sep(nonzero_cts=nonzero_cts, batch_sep=batch_sep, group_sep=group_sep, glm_ctrl=glm_ctrl, batch_name="Brazil_1")
disp_nonzero_bzl2 <- genDisp_Sep(nonzero_cts=nonzero_cts, batch_sep=batch_sep, group_sep=group_sep, glm_ctrl=glm_ctrl, batch_name="Brazil_2")



####  Compare
disp_full_lst <- list(Africa=disp_full_africa, G6=disp_full_G6, India=disp_full_india,
                      Brazil_1=disp_full_bzl1, Brazil_2=disp_full_bzl2)
disp_nonzero_lst <- list(Africa=disp_nonzero_africa, G6=disp_nonzero_G6, India=disp_nonzero_india,
                         Brazil_1=disp_nonzero_bzl1, Brazil_2=disp_nonzero_bzl2)

disp_full_stats <- lapply(c(min, median, mean, max), function(ff){sapply(disp_full_lst, ff)})
disp_full_stats <- do.call(rbind, disp_full_stats)
disp_nonzero_stats <- lapply(c(min, median, mean, max), function(ff){sapply(disp_nonzero_lst, ff)})
disp_nonzero_stats <- do.call(rbind, disp_nonzero_stats)
rownames(disp_full_stats) <- rownames(disp_nonzero_stats) <- c("min disp.", "median disp.", "mean disp.", "max disp.")

print(round(disp_full_stats, 4))
print(round(disp_nonzero_stats, 4))


####  Save result
save(disp_full_africa, disp_full_G6, disp_full_india, disp_full_bzl1, disp_full_bzl2,
     disp_nonzero_africa, disp_nonzero_G6, disp_nonzero_india, disp_nonzero_bzl1, disp_nonzero_bzl2,
     file="disps_nonzero.RData")
