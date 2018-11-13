rm(list=ls())
setwd("~/Documents/ComBat_seq/real_data_app/")
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/TB"

# load data
rds_obj <- readRDS(file.path(data_dir, "combined.rds"))

# take brazil batch 1 & india
batch_filter <- rds_obj$SequencingBatch %in% c("Brazil_1", "India")
group_filter <- rds_obj$Label %in% c("Non-progressor", "Active")
rds_obj_subset <- rds_obj[, batch_filter & group_filter]

# count matrix
cts <- assays(rds_obj_subset)$counts

# batch indicator
batch_org <- as.character(colData(rds_obj_subset)$SequencingBatch)
batch <- rep(1, length(batch_org)); batch[batch_org=="India"] <- 2

# condition indicator
group_org <- as.character(colData(rds_obj_subset)$Label)
group <- rep(0, length(group_org)); group[group_org=="Active"] <- 1

# remove genes with only zero in any of batches
keep1 <- apply(cts[, batch==1], 1, function(x){!all(x==0)})
keep2 <- apply(cts[, batch==2], 1, function(y){!all(y==0)})
cts <- cts[keep1 & keep2, ]  

# get covariate for gender
covar <- rds_obj_subset$Sex
               
# save data
save(cts, batch, group, covar, file="TB_ExpRep.RData")
