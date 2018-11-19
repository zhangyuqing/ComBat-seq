rm(list=ls())
setwd("~/Dropbox/Work/ComBat_Seq/ComBat-Seq/explore_real_data/SingleCell/")
sapply(c("ggplot2", "reshape2", "gridExtra", "scales"), require, character.only=TRUE)
data_dir <- "~/Google Drive/ComBat_seq/real_data_example/SingleCell/Tung2017Batch_ScientificReports"


cts <- read.table(file.path(data_dir, "molecules.txt"), sep = "\t")  # use molecule counts rather than read counts
anno <- read.table(file.path(data_dir, "annotation.txt"), sep = "\t", header = TRUE)
# ckeep <- anno[, "individual"] == "NA19101" | anno[, "individual"] == "NA19239"
# cts <- as.matrix(cts[, ckeep])
# group <- as.factor(as.character(anno[ckeep, "individual"]))
# batch <- as.factor(as.character(anno[ckeep, "batch"]))
#batch <- anno[ckeep, "replicate"]
# remove genes that aren't expressed in at least 6 cells


n_zeros_in_genes <- apply(cts, 1, function(x){length(which(x==0))})
percent_zeros_in_genes <- n_zeros_in_genes / ncol(cts)
percent_zeros_df <- data.frame(genes=rownames(cts), value=percent_zeros_in_genes)

binned_cts <- hist(percent_zeros_in_genes, breaks=10, plot=FALSE)$counts 
binned_proportions <- binned_cts / sum(binned_cts)
names(binned_proportions) <- paste0("<", percent(seq(from=0.1, to=1, by=0.1)))  
cat("######  Proportion of genes with X% zeros across samples  #######\n")
print(round(binned_proportions,3))

ggplot(percent_zeros_df, aes(0, value)) + 
  geom_violin() +
  geom_boxplot(aes(0.7, value), width=0.15) +
  coord_flip() +
  annotate(geom="text", label=percent(median(percent_zeros_df$value)), x=0.85, y=median(percent_zeros_df$value)) +
  labs(y="% Zeros", title="Proportion of zeros across samples for the remaining genes") +
  theme(axis.title.y=element_blank())#, axis.text.y=element_blank())

