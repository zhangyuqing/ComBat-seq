rm(list=ls())
sapply(c("plyr", "SummarizedExperiment", "edgeR", "ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)
source("../../../ComBat_seq_returnParams.R")
source("../../../helper_seq.R")

## Load data
data_name <- "lung"

if(data_name=="lung"){
  data_dir <- "~/Google Drive/ComBat_seq/real_data_example/RNAseq/WRBU_lungcancer"
  cts <- readRDS(file.path(data_dir, "count_matrix_BU_WR.rds"))
  meta_info <- readRDS(file.path(data_dir, "BU_WR_annotation_file.rds"))
  # drop 4790-024 (single sample as its own batch)
  cts <- cts[, -grep("4790-024", colnames(cts))]
  meta_info <- filter(meta_info, kitnumber!="4790-024")
  batch <- factor(meta_info$site)
  group <- factor(meta_info$smoking_status)
  group <- as.factor(as.character(revalue(group, c("Former smoker"="0", "Current smoker"="1"))))
  # remove all 0 genes
  tmp <- t(apply(cts,1,floor)); colnames(tmp) <- colnames(cts); cts <- tmp; rm(tmp)
  cts <- cts[rowVars(cts)>0, ]
  plt_title <- "Lung cancer"
}else if(data_name=="CHD8"){
  load("~/Google Drive/ComBat_seq/real_data_example/RNAseq/CHD8/rse_gene.Rdata")
  cts <- assays(rse_gene)$counts
  cts <- cts[rowVars(cts)>0, ]
  batch <- colData(rse_gene)$batch
  group <- colData(rse_gene)$condition_bi
  plt_title <- "CHD8"
}


## Run combat-seq
combatseq_res <- ComBat_seq(counts=cts, batch=batch, group=group, shrink=FALSE, shrink.disp=FALSE)


## Genes strongly affected by batch
genes_bup_1 <- which(combatseq_res$gamma_hat[, 1] > 2.5); names(genes_bup_1) <- NULL
genes_bdown_1 <- which(combatseq_res$gamma_hat[, 1] < -2.5); names(genes_bdown_1) <- NULL
genes_bup_2 <- which(combatseq_res$gamma_hat[, 2] > 2.5); names(genes_bup_2) <- NULL
genes_bdown_2 <- which(combatseq_res$gamma_hat[, 2] < -2.5); names(genes_bdown_2) <- NULL

identical(intersect(genes_bup_1, genes_bdown_2), genes_bup_1)
identical(intersect(genes_bup_2, genes_bdown_1), genes_bup_2)

batch_gene_lst <- list(bup1=genes_bup_1, bup2=genes_bup_2)

cts_norm <- apply(cts, 2, function(x){x/sum(x)})

# gene-wise line plot
# cts_norm_tidy_lst <- lapply(1:length(batch_gene_lst), function(ii){
#   cts_norm_sub <- cts_norm[batch_gene_lst[[ii]], ]
#   cts_norm_df <- data.frame(Sample=colnames(cts_norm_sub), Batch=batch, t(cts_norm_sub))
#   cts_norm_tidy <- melt(cts_norm_df, id.vars=c("Sample", "Batch"), variable.name="Gene")
#   cts_norm_tidy$Type <- ii
#   return(cts_norm_tidy)
# })
# cts_norm_tidy <- data.frame(do.call(rbind, cts_norm_tidy_lst))
# cts_norm_tidy$Type <- plyr::revalue(factor(cts_norm_tidy$Type), c("1"="Up in batch 1", "2"="Up in batch 2"))
# ggplot(cts_norm_tidy, aes(x=Sample, y=value, group=Gene, color=Type))+
#   geom_line() +
#   #facet_wrap(~Batch) +
#   ylim(0, 1e-07) +
#   theme(axis.text.x=element_text(angle=45, hjust=1))


# average gene line plot
avg_cts_norm_tidy_lst <- lapply(1:length(batch_gene_lst), function(ii){
  avg_cts_norm_sub <- colMeans(cts_norm[batch_gene_lst[[ii]], ])
  avg_cts_norm_df <- data.frame(Sample=colnames(cts_norm), Batch=batch, Group=group,
                                AvgComp=avg_cts_norm_sub, Type=ii)
  return(avg_cts_norm_df)
})
avg_cts_norm_tidy <- data.frame(do.call(rbind, avg_cts_norm_tidy_lst))
avg_cts_norm_tidy$Type <- plyr::revalue(factor(avg_cts_norm_tidy$Type), c("1"="Up in batch 1", "2"="Up in batch 2"))

avg_cts_norm_tidy <- arrange(avg_cts_norm_tidy, Batch, Group)
avg_cts_norm_tidy$Sample <- factor(avg_cts_norm_tidy$Sample, levels=colnames(cts)[order(paste(batch,group))])

png(sprintf("figures/libcomp_%s.png", data_name), width=6, height=5, units="in", res=300)
ggplot(avg_cts_norm_tidy, aes(x=Sample, y=AvgComp, group=Type, color=Type, shape=Group))+
  geom_line() +
  geom_point() +
  #ylim(0, 1e-07) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(y="Average composition of genes with large mean batch effect",
       title=plt_title)
dev.off()



## Genes affected by batch
genes_bup_1 <- which(combatseq_res$gamma_hat[, 1] > 0); names(genes_bup_1) <- NULL
genes_bdown_1 <- which(combatseq_res$gamma_hat[, 1] < 0); names(genes_bdown_1) <- NULL
pars_df <- data.frame(gamma=combatseq_res$gamma_hat[,2], phi=combatseq_res$phi_hat[,2], b1.type="Null")
pars_df$b1.type <- as.character(pars_df$b1.type)
pars_df$b1.type[1:nrow(cts) %in% genes_bup_1] <- "gamma > 0 in batch 1"
pars_df$b1.type[1:nrow(cts) %in% genes_bdown_1] <- "gamma < 0 in batch 1"
pars_df$b1.type <- factor(pars_df$b1.type, levels=c("gamma > 0 in batch 1", "gamma < 0 in batch 1", "Null"))
png(sprintf("figures/libcomp_full_%s.png", data_name), width=6, height=5, units="in", res=300)
ggplot(pars_df, aes(x=gamma, y=phi, color=b1.type)) +
  geom_point() +
  labs(title="Batch 2 parameter estimates")
dev.off()
