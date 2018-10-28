####  The dispersion estimates using non-zero portion of genes take too...ooo long to run. 
####  To save time, I'll run it once and save it as RData to avoid waiting too long each time I need to recreate the report
rm(list=ls())
setwd("~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB")
sapply(c("SummarizedExperiment", "DESeq2", "edgeR", "dendextend", "ggplot2", "reshape2", "gridExtra", "scales", 
         "ggdendro", "MASS"), require, character.only=TRUE)

####  Load data
rds_obj <- readRDS(file.path(output_dir, "combined.rds"))
cts <- assays(rds_obj)$counts
batch <- colData(rds_obj)$SequencingBatch
group <- colData(rds_obj)$Label
glm_ctrl <- glm.control(maxit=200)


####  Define function to compute dispersion


####  Use counts WITH zeros to estimate dispersion
# Africa
disp_full_batch1 <- sapply(1:nrow(cts), function(kk){
  curr_cts <- cts[kk, batch==1] 
  if(length(curr_cts) > 2 & var(curr_cts) > 0){
    curr_f <- try(glm.nb(curr_cts ~ as.factor(group_num[batch==1]), control=glm_ctrl), silent=TRUE)
    if(class(curr_f)[1]!="try-error"){
      return(1/summary(curr_f)$theta)
    }
  }
  return(NA)
})
disp_full_batch1 <- disp_full_batch1[!is.na(disp_full_batch1)]

# G6
disp_full_batch2 <- sapply(1:nrow(cts), function(ll){
  curr_cts <- cts[ll, batch==2] 
  if(length(curr_cts) > 2 & var(curr_cts) > 0){
    curr_f <- try(glm.nb(curr_cts ~ as.factor(group_num[batch==2]), control=glm_ctrl), silent=TRUE)
    if(class(curr_f)[1]!="try-error"){
      return(1/summary(curr_f)$theta)
    }
  }
  return(NA)
})
disp_full_batch2 <- disp_full_batch2[!is.na(disp_full_batch2)]




####  Use only non-zero portions of each gene to estimate dispersion
##  Take the non-zero portion in each gene
batch_sep <- group_sep <- nonzero_cts <- list()
for(i in 1:nrow(cts)){
  nonzero_ind <- which(cts[i, ]!=0)
  batch_sep[[i]] <- batch[nonzero_ind]
  group_sep[[i]] <- group_num[nonzero_ind]
  nonzero_cts[[i]] <- cts[i, nonzero_ind]
}
# sanity check
if(!identical(sapply(nonzero_cts, length), as.integer(rowSums(cts!=0)))){stop("Error in taking non-zero portions!")}
names(nonzero_cts) <- rownames(cts); names(batch_sep) <- rownames(cts); names(group_sep) <- rownames(cts)

##  Dispersion in batch 1
# take batch 1 data
group_sep_batch1 <- nonzero_cts_batch1 <- list()
for(i in 1:length(nonzero_cts)){
  group_sep_batch1[[i]] <- group_sep[[i]][batch_sep[[i]]==1]
  nonzero_cts_batch1[[i]] <- nonzero_cts[[i]][batch_sep[[i]]==1]
}
names(nonzero_cts_batch1) <- names(group_sep_batch1) <- rownames(cts)
# calculate dispersions
disp_nonzero_batch1 <- sapply(1:length(nonzero_cts_batch1), function(ii){
  curr_cts <- nonzero_cts_batch1[[ii]] 
  curr_group <- group_sep_batch1[[ii]]
  if(length(curr_cts) > 2 & var(curr_cts) > 0){
    curr_f <- try(glm.nb(curr_cts ~ as.factor(curr_group), control=glm_ctrl), silent=TRUE)
    if(class(curr_f)[1]!="try-error"){
      return(1/summary(curr_f)$theta)
    }
  }
  # curr_y <- DGEList(counts=t(as.matrix(curr_cts)), group=as.factor(curr_group))
  # if(length(curr_cts) > 3 & var(curr_cts) > 0){
  #   if(length(unique(curr_group))>1){
  #     curr_design <- model.matrix(~as.factor(curr_group))
  #   }else{
  #     curr_design <- NULL
  #   }
  #   curr_y <- try(estimateDisp(curr_y, curr_design), silent=TRUE)
  #   if(class(curr_y)[1]!="try-error"){return(curr_y$tagwise.dispersion)}
  # }
  return(NA)
})
disp_nonzero_batch1 <- disp_nonzero_batch1[!is.na(disp_nonzero_batch1)]

##  Dispersion in batch 2
#take batch 2 data
group_sep_batch2 <- nonzero_cts_batch2 <- list()
for(j in 1:length(nonzero_cts)){
  group_sep_batch2[[j]] <- group_sep[[j]][batch_sep[[j]]==2]
  nonzero_cts_batch2[[j]] <- nonzero_cts[[j]][batch_sep[[j]]==2]
}
names(nonzero_cts_batch2) <- names(group_sep_batch2) <- rownames(cts)
# calculate dispersions
disp_nonzero_batch2 <- sapply(1:length(nonzero_cts_batch2), function(jj){
  curr_cts <- nonzero_cts_batch2[[jj]] 
  curr_group <- group_sep_batch2[[jj]]
  if(length(curr_cts) > 2 & var(curr_cts) > 0){
    curr_f <- try(glm.nb(curr_cts ~ as.factor(curr_group), control=glm_ctrl), silent=TRUE)
    if(class(curr_f)[1]!="try-error"){
      return(1/summary(curr_f)$theta)
    }
  }
  return(NA)
})
disp_nonzero_batch2 <- disp_nonzero_batch2[!is.na(disp_nonzero_batch2)]


####  Compare
disp_full_lst <- list(Batch1=disp_full_batch1, Batch2=disp_full_batch2)
disp_nonzero_lst <- list(Batch1=disp_nonzero_batch1, Batch2=disp_nonzero_batch2)

disp_full_stats <- lapply(c(min, median, mean, max), function(ff){sapply(disp_full_lst, ff)})
disp_full_stats <- do.call(rbind, disp_full_stats)
disp_nonzero_stats <- lapply(c(min, median, mean, max), function(ff){sapply(disp_nonzero_lst, ff)})
disp_nonzero_stats <- do.call(rbind, disp_nonzero_stats)
rownames(disp_full_stats) <- rownames(disp_nonzero_stats) <- c("min disp.", "median disp.", "mean disp.", "max disp.")

print(round(disp_full_stats, 4))
print(round(disp_nonzero_stats, 4))


####  Save result
save(disp_nonzero_batch1, disp_nonzero_batch2, disp_full_batch1, disp_full_batch2, file="disps_nonzero.RData")
