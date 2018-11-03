rm(list=ls())
setwd("/restricted/projectnb/combat/data/TB_merged")
sapply(c("SummarizedExperiment", "DESeq2", "edgeR", "scales", "MASS"), require, character.only=TRUE)


####  Load data
rds_obj <- readRDS("combined.rds")
cts <- assays(rds_obj)$counts
# batch indicator
batch <- colData(rds_obj)$SequencingBatch
# condition indicator
group <- colData(rds_obj)$Label



####  Africa vs G6 
## take sub-studies
cts_AG <- cts[, rds_obj$Dataset %in% c("Africa","G6")]
group_AG <- group[rds_obj$Dataset %in% c("Africa","G6")]
batch_AG <- as.factor(as.character(batch[rds_obj$Dataset %in% c("Africa","G6")]))
# Batch 1 - Africa; Batch 2 - G6


## dispersion in batches & in whole study (including batch as covariate) (edgeR)
y_batch1 <- DGEList(counts=cts_AG[, batch_AG=="Africa"], group=group_AG[batch_AG=="Africa"])
y_batch1 <- calcNormFactors(y_batch1, method="TMM")
design_batch1 <- model.matrix(~group_AG[batch_AG=="Africa"])
y_batch1 <- estimateDisp(y_batch1, design_batch1)

y_batch2 <- DGEList(counts=cts_AG[, batch_AG=="G6"], group=as.factor(as.character(group_AG[batch_AG=="G6"])))
y_batch2 <- calcNormFactors(y_batch2, method="TMM")
design_batch2 <- model.matrix(~as.factor(as.character(group_AG[batch_AG=="G6"])))
y_batch2 <- estimateDisp(y_batch2, design_batch2)

y_whole <- DGEList(counts=cts_AG)
y_whole <- calcNormFactors(y_whole, method="TMM")
design_whole <- model.matrix(~ group_AG + batch_AG)
y_whole <- estimateDisp(y_whole, design_whole)

disp1_gene_edgeR <- y_batch1$tagwise.dispersion
disp2_gene_edgeR <- y_batch2$tagwise.dispersion
dispW_gene_edgeR <- y_whole$tagwise.dispersion

disp1_trend_edgeR <- y_batch1$trended.dispersion
disp2_trend_edgeR <- y_batch2$trended.dispersion
dispW_trend_edgeR <- y_whole$trended.dispersion


## dispersion in batches & in whole study (including batch as covariate) (DESeq2)
dds1 <- DESeqDataSetFromMatrix(countData=cts_AG[, batch_AG=="Africa"], design=~Group,
                               colData=data.frame(Group=group_AG[batch_AG=="Africa"]))
dds1 <- estimateSizeFactors(dds1)
dds1 <- estimateDispersions(dds1)

dds2 <- DESeqDataSetFromMatrix(countData=cts_AG[, batch_AG=="G6"], design=~Group,
                               colData=data.frame(Group=as.factor(as.character(group_AG[batch_AG=="G6"]))))
dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)

ddsW <- DESeqDataSetFromMatrix(countData=cts_AG, design=~Group + Batch,
                               colData=data.frame(Group=group_AG, Batch=batch_AG))
ddsW <- estimateSizeFactors(ddsW)
ddsW <- estimateDispersions(ddsW)

disp1_gene_DESeq2 <- dispersions(dds1)
disp2_gene_DESeq2 <- dispersions(dds2)
dispW_gene_DESeq2 <- dispersions(ddsW)


## summarize in data frame
disps_all <- data.frame(Batch1.edgeR=disp1_gene_edgeR, 
                        Batch2.edgeR=disp2_gene_edgeR, 
                        Whole.edgeR=dispW_gene_edgeR,
                        Batch1.edgeR.trend=disp1_trend_edgeR, 
                        Batch2.edgeR.trend=disp2_trend_edgeR, 
                        Whole.edgeR.trend=dispW_trend_edgeR,
                        Batch1.DESeq2=disp1_gene_DESeq2, 
                        Batch2.DESeq2=disp2_gene_DESeq2, 
                        Whole.DESeq2=dispW_gene_DESeq2)

saveRDS(disps_all, file="TB_disps_AG.rds")



####  Brazil 1 vs India
rm(disps_all)

## take sub-studies
cts_BI <- cts[, rds_obj$SequencingBatch %in% c("Brazil_1","India")]
group_BI <- group[rds_obj$SequencingBatch %in% c("Brazil_1","India")]
batch_BI <- as.factor(as.character(batch[rds_obj$SequencingBatch %in% c("Brazil_1","India")]))
# Batch 1 -  Brazil 1; Batch 2 - India
group_BI_india <- factor(as.character(group_BI[batch_BI=="India"]), levels=c("Non-progressor", "Active"))

## dispersion in batches & in whole study (including batch as covariate) (edgeR)
y_batch1 <- DGEList(counts=cts_BI[, batch_BI=="Brazil_1"], group=group_BI[batch_BI=="Brazil_1"])
y_batch1 <- calcNormFactors(y_batch1, method="TMM")
design_batch1 <- model.matrix(~group_BI[batch_BI=="Brazil_1"])
y_batch1 <- estimateDisp(y_batch1, design_batch1)

y_batch2 <- DGEList(counts=cts_BI[, batch_BI=="India"], group=group_BI_india)
y_batch2 <- calcNormFactors(y_batch2, method="TMM")
design_batch2 <- model.matrix(~group_BI_india)
y_batch2 <- estimateDisp(y_batch2, design_batch2)

y_whole <- DGEList(counts=cts_BI)
y_whole <- calcNormFactors(y_whole, method="TMM")
design_whole <- model.matrix(~ group_BI + batch_BI)
y_whole <- estimateDisp(y_whole, design_whole)

disp1_gene_edgeR <- y_batch1$tagwise.dispersion
disp2_gene_edgeR <- y_batch2$tagwise.dispersion
dispW_gene_edgeR <- y_whole$tagwise.dispersion

disp1_trend_edgeR <- y_batch1$trended.dispersion
disp2_trend_edgeR <- y_batch2$trended.dispersion
dispW_trend_edgeR <- y_whole$trended.dispersion


## dispersion in batches & in whole study (including batch as covariate) (DESeq2)
dds1 <- DESeqDataSetFromMatrix(countData=cts_BI[, batch_BI=="Brazil_1"], design=~Group,
                               colData=data.frame(Group=group_BI[batch_BI=="Brazil_1"]))
dds1 <- estimateSizeFactors(dds1)
dds1 <- estimateDispersions(dds1)

dds2 <- DESeqDataSetFromMatrix(countData=cts_BI[, batch_BI=="India"], design=~Group,
                               colData=data.frame(Group=group_BI_india))
dds2 <- estimateSizeFactors(dds2)
dds2 <- estimateDispersions(dds2)

ddsW <- DESeqDataSetFromMatrix(countData=cts_BI, design=~ Group + Batch,
                               colData=data.frame(Group=group_BI, Batch=batch_BI))
ddsW <- estimateSizeFactors(ddsW)
ddsW <- estimateDispersions(ddsW)

disp1_gene_DESeq2 <- dispersions(dds1)
disp2_gene_DESeq2 <- dispersions(dds2)
dispW_gene_DESeq2 <- dispersions(ddsW)

## summarize in data frame
disps_all <- data.frame(Batch1.edgeR=disp1_gene_edgeR, 
                        Batch2.edgeR=disp2_gene_edgeR, 
                        Whole.edgeR=dispW_gene_edgeR,
                        Batch1.edgeR.trend=disp1_trend_edgeR, 
                        Batch2.edgeR.trend=disp2_trend_edgeR, 
                        Whole.edgeR.trend=dispW_trend_edgeR,
                        Batch1.DESeq2=disp1_gene_DESeq2, 
                        Batch2.DESeq2=disp2_gene_DESeq2, 
                        Whole.DESeq2=dispW_gene_DESeq2)

saveRDS(disps_all, file="TB_disps_BI.rds")

