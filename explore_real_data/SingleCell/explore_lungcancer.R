rm(list=ls())
sapply(c("Biobase", "dplyr"), require, character.only=TRUE)  # load libraries
setwd("/restricted/projectnb/combat/work/yuqingz/ComBat-seq")  # set working directory, change into your own path

# paths for the data
# batch1_dir <- "/restricted/projectnb/pulmseq/zhe/CEL_Seq_2_20170325/Count_results"
# batch2_dir <- "/restricted/project/pulmseq/Single_Cell_Test/Experiments/2017_06_29_IDA_2/Results/IDA2_Single_Cell_ExpressionSet_V2.rds"
data_dir <- "/restricted/projectnb/combat/data/single_cell_lung_cancer_IDA_Jen"
batch1_dir <- file.path(data_dir, "Count_results")
batch2_dir <- file.path(data_dir, "IDA2_Single_Cell_ExpressionSet_V2.rds")

# load batch 1
flist <- dir(batch1_dir)  # data files in batch 1 directory
dat1_lst <- lapply(flist, function(fname){
  tmp <- read.delim(file.path(batch1_dir, fname))
  rownames(tmp) <- tmp[,1]; tmp <- tmp[, -1]; return(tmp)
})  # loop through files in batch 1 directory and read data one by one
cts1 <- as.matrix(do.call(cbind, dat1_lst))  # bind all sample files in batch 1 into one count matrix

# load batch 2
dat2 <- readRDS(batch2_dir)  # read .rds file for batch 2
cts2 <- exprs(dat2)  # get count matrix for batch 2

# information for all samples
meta_df <- pData(dat2)

# take common genes of two batches
common_genes <- intersect(rownames(cts1), rownames(cts2))  # take overlap of genes across 2 batches
cts1 <- cts1[common_genes, ]
cts2 <- cts2[common_genes, ]

# combine the two batches
cts <- cbind(cts1, cts2)  

# get batch information and sample ID from column names
sample_colnames <- colnames(cts)  # column names of count matrix contain batch variable & sample ID   
sample_colnames_splt <- regmatches(sample_colnames, regexpr("_", sample_colnames), invert=TRUE)  # split names by first "_", into batch and sample IDs
batch <- sapply(sample_colnames_splt, function(x){x[[1]]})
sample_ID <- sapply(sample_colnames_splt, function(x){x[[2]]})  # all sample IDs in count matrix

# subsetting - get 3 former smokers no cancer & 3 former smokers with ADC & only cell samples
meta_sub <- filter(meta_df, status_today=="former", cancer=="N" | cancer_subtype=="ADC", sample.type=="Cell")  # meta information for subsetted samples
selected_sample_id <- gsub("IDA2_", "", meta_sub$ID)  # IDs for selected samples

# take subset of selected samples
batch <- batch[which(sample_ID %in% selected_sample_id)]
subset_id <- match(sample_ID, selected_sample_id); subset_id <- subset_id[!is.na(subset_id)]  
meta_sub_matched <- meta_sub[subset_id, ]
cts <- cts[, paste(batch, gsub("IDA2_", "", meta_sub_matched$ID), sep="_")]
meta_sub_matched$batch <- batch 

# save the R data file (.RData) for download
save(cts, meta_sub_matched, file="single_cell_lung_cancer.RData")

# download and organize into a SCE object locally
rm(list=ls())
load("~/Downloads/single_cell_lung_cancer.RData")
library(singleCellTK)
rownames(meta_sub_matched) <- paste(meta_sub_matched$batch, gsub("IDA2_", "", meta_sub_matched$ID), sep="_")
sce_obj <- createSCE(assayFile = cts, assayName = "counts",
                     annotFile = meta_sub_matched, 
                     inputDataFrames = TRUE, createLogCounts = TRUE)
saveRDS(sce_obj, file="single_cell_lung_cancer.rds")
