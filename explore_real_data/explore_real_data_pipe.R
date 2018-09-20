rm(list=ls())
setwd("C:/Users/zhang/Google Drive/ComBat_seq/real_data_example/")
pkg_vec <- c("DESeq2", "edgeR","dendextend", "ggplot2")
sapply(pkg_vec, require, character.only=TRUE)
print(sapply(pkg_vec, function(pkg_name){as.character(packageVersion(pkg_name))}))


####  Load desired data
command_args <- commandArgs(trailingOnly=TRUE)
study_SRP <- command_args[1]
test_batch <- as.numeric(command_args[2])
test_cond <- as.numeric(command_args[3])
# study_SRP <- "SRP047233"; test_batch <- 1; test_cond <- 0
load(file.path(study_SRP, 'rse_gene.Rdata'))  

batch <- as.numeric(as.character(rse_gene$batch))
cond_bi <- as.numeric(as.character(rse_gene$condition_bi))


####  Visualize to see if there is batch effect
p <- plotPCA(DESeqTransform(rse_gene), intgroup="batch")
ggsave(file.path(study_SRP, 'PCA_batch.png'), p, width=9, height=5, units="in")
dev.off(); rm(p)

p2 <- plotPCA(DESeqTransform(rse_gene), intgroup="condition_bi")
ggsave(file.path(study_SRP, 'PCA_cond.png'), p2, width=9, height=5, units="in")
dev.off(); rm(p2)

counts_norm <- scale(t(assays(rse_gene)$counts), center=TRUE, scale=TRUE)
hc <- hclust(dist(counts_norm))
dend <- as.dendrogram(hc)

png(file.path(study_SRP, 'hclust_batch.png'), width=7, height=7, units="in", res=300)
dend_batch <- color_branches(dend, col=batch[order.dendrogram(dend)]+1)
plot(dend_batch) 
dev.off(); rm(dend_batch)

png(file.path(study_SRP, 'hclust_cond.png'), width=7, height=7, units="in", res=300)
dend_cond <- color_branches(dend, col=cond_bi[order.dendrogram(dend)]+2)
plot(dend_cond) 
dev.off(); rm(dend_cond)


####  Distribution of library size
lib_sizes <- colSums(assays(rse_gene)$counts)
lib_sizes_df <- data.frame(lib_sizes=lib_sizes, batch=as.factor(batch))

# overall sample library sizes
png(file.path(study_SRP, 'lib_sizes_hist.png'), width=5.5, height=5.5, units="in", res=300)
ggplot(lib_sizes_df, aes(x=lib_sizes)) + 
  geom_histogram(bins=round(ncol(rse_gene)/5)) +
  labs(x="Library sizes in samples", y="Number of samples", 
       title="Distribution of observed library sizes") +
  geom_vline(aes(xintercept=mean(lib_sizes)), color="red", linetype="dashed", size=1) +
  annotate(geom="text", x=mean(lib_sizes), y=max(hist(lib_sizes)$counts)+1, 
           label=round(mean(lib_sizes),0), color="red")
dev.off()

print(round(range(lib_sizes),0))

# boxplot comparing distribution of library size across batches
png(file.path(study_SRP, 'lib_sizes_batch_boxplot.png'), width=5.5, height=5.5, units="in", res=300)
ggplot(lib_sizes_df, aes(x=batch, y=lib_sizes, fill=batch)) + 
  geom_boxplot() +
  labs(x="", y="Library sizes", title="Distribution of observed library sizes across batch") 
dev.off()  


####  Differential genes relative to biological condition (within batch 1)
counts_batch <- assays(rse_gene)$counts[, batch==test_batch]
cond_batch <- cond_bi[batch==test_batch]

y <- DGEList(counts=counts_batch, group=cond_batch)
y <- calcNormFactors(y)
design <- model.matrix(~as.factor(cond_batch))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
de_res <- topTags(qlf, n=nrow(rse_gene))

up_genes <- rownames(de_res)[de_res$table$logFC >= 0]  # all up-regulated genes
head_up_genes <- head(up_genes, n=50)  # top-50 up-regulated genes
tail_up_genes <- tail(up_genes, n=50)  # 50 least up-regulated genes

down_genes <- rownames(de_res)[de_res$table$logFC < 0]  # all down-regulated genes
head_down_genes <- head(down_genes, n=50)  # top-50 down-regulated genes
tail_down_genes <- tail(down_genes, n=50)  # 50 least down-regulated genes

# changes in terms of mean
print(mean(counts_batch[head_up_genes, cond_batch==0]));
print(mean(counts_batch[head_up_genes, cond_batch==1]))

print(mean(counts_batch[tail_up_genes, cond_batch==0]));
print(mean(counts_batch[tail_up_genes, cond_batch==1]))

print(mean(counts_batch[head_down_genes, cond_batch==0]));
print(mean(counts_batch[head_down_genes, cond_batch==1]))

print(mean(counts_batch[tail_down_genes, cond_batch==0]));
print(mean(counts_batch[tail_down_genes, cond_batch==1]))

# changes in terms of fold change
print(round(range(exp(de_res$table[head_up_genes, "logFC"])), 1));
print(round(range(exp(de_res$table[tail_up_genes, "logFC"])), 1))

print(round(range(exp(de_res$table[head_down_genes, "logFC"])), 1));
print(round(range(exp(de_res$table[tail_down_genes, "logFC"])), 1))

rm()


## Differential genes relative to batch (within condition 0)
counts_cond <- assays(rse_gene)$counts[, cond_bi==test_cond]
batch_cond <- batch[cond_bi==test_cond]

y <- DGEList(counts=counts_batch, group=batch_cond)
y <- calcNormFactors(y)
design <- model.matrix(~as.factor(batch_cond))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)
de_res <- topTags(qlf, n=nrow(rse_gene))

up_genes <- rownames(de_res)[de_res$table$logFC >= 0]  # all up-regulated genes
head_up_genes <- head(up_genes, n=50)  # top-50 up-regulated genes
tail_up_genes <- tail(up_genes, n=50)  # 50 least up-regulated genes

down_genes <- rownames(de_res)[de_res$table$logFC < 0]  # all down-regulated genes
head_down_genes <- head(down_genes, n=50)  # top-50 down-regulated genes
tail_down_genes <- tail(down_genes, n=50)  # 50 least down-regulated genes


## Dispersion differences across batch 1 and 2
y_batch1 <- DGEList(counts=assays(rse_gene)$counts[, batch==1], group=cond_bi[batch==1])
y_batch1 <- calcNormFactors(y_batch1)
design_batch1 <- model.matrix(~as.factor(cond_bi[batch==1]))
y_batch1 <- estimateDisp(y_batch1, design_batch1)

y_batch2 <- DGEList(counts=assays(rse_gene)$counts[, batch==2], group=cond_bi[batch==2])
y_batch2 <- calcNormFactors(y_batch2)
design_batch2 <- model.matrix(~as.factor(cond_bi[batch==2]))
y_batch2 <- estimateDisp(y_batch2, design_batch2)

