rm(list=ls())
setwd("C:/Users/zhang/Google Drive/ComBat_seq/real_data_example/")
pkg_vec <- c("recount")
sapply(pkg_vec, require, character.only=TRUE)
print(sapply(pkg_vec, function(pkg_name){as.character(packageVersion(pkg_name))}))

##  Load data
# abstract_search('batch')$project
study_SRP <- "SRP047233"
download_study(study_SRP)
load(file.path(study_SRP, 'rse_gene.Rdata'))
meta_info <- as.data.frame(colData(rse_gene))

## Obtain read counts from provided base pair coverage counts
rse_gene <- read_counts(rse_gene, round=TRUE) #scale_counts(rse_gene)  

## Remove genes with all 0 counts
rse_gene <- rse_gene[apply(assays(rse_gene)$counts, 1, function(x){!all(x==0)}), ]

## Batch indicator
# split_info <- strsplit(meta_info$characteristics," ")
# batch <- rep(NA, length(split_info))
# for(i in 1:length(split_info)){
#   batch[i] <- split_info[[i]][grep("batch:", split_info[[i]])+1]
# }
# batch <- as.numeric(gsub(')', '', batch))
batch <- sapply(meta_info$characteristics, function(s){s[grep("batch",s)]})
batch <- as.factor(batch)
batch <- sapply(batch, function(b){return(which(levels(batch)==b))})

## Condition indicator
group <- rep(NA, nrow(meta_info))
group[grep("CHD8", meta_info$title)] <- "CHD8"
group[grep("GFP", meta_info$title)] <- "GFP"
group[grep("LacZ", meta_info$title)] <- "LacZ"
group_num <- rep(0, nrow(meta_info)); group_num[group=="CHD8"] <- 1
#group <- sapply(meta_info$characteristics, function(s){s[grep("shRNA target gene",s)]})
  
## Add batch and condition in meta info
colData(rse_gene)$batch <- as.factor(batch)
colData(rse_gene)$condition <- as.factor(group)
colData(rse_gene)$condition_bi <- as.factor(group_num)

##  Save data in desired format
save(rse_gene, file=file.path(study_SRP, 'rse_gene.Rdata'))
