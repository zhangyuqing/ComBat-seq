sapply(c("ggplot2", "reshape2", "gridExtra", "dendextend", "edgeR", "DESeq2"), require, character.only=TRUE)

########### Control only
cts_norm <- apply(cts_ctrl, 2, function(x){x/sum(x)})
#combatseq_ctrl_f <- combatseq_ctrl[apply(combatseq_ctrl, 1, function(x){all(!is.na(x))}), ]
combatseq_norm <- apply(combatseq_ctrl, 2, function(x){x/sum(x)})

###  PCA
pcres1 <- prcomp(t(cts_norm))
pcres1_df <- data.frame(pcres1$x[,1:2])
pcres2 <- prcomp(t(combatseq_norm))
pcres2_df <- data.frame(pcres2$x[,1:2])

plt_df <- data.frame(rbind(pcres1_df, pcres2_df),
                     Batch=factor(rep(batch_ctrl,2)),
                     Method=factor(c(rep("Before",nrow(pcres1_df)), rep("After",nrow(pcres2_df))), 
                                   levels=c("Before", "After")))
ggplot(plt_df, aes(PC1, PC2, color=Batch)) + 
  geom_point() +
  facet_wrap(~Method)

###  clustering
hc_base <- hclust(dist(t(cts_norm)))
dend_base <- as.dendrogram(hc_base)
dend_base <- color_branches(dend_base, col=batch_ctrl[order.dendrogram(dend_base)]) 
plot(dend_base)  

adj_hc <- hclust(dist(t(combatseq_norm)))
adj_dend <- as.dendrogram(adj_hc)
adj_dend <- color_branches(adj_dend, col=batch_ctrl[order.dendrogram(adj_dend)])
plot(adj_dend) 


########### Whole dataset
cts_mat <- cts_mat[apply(cts_mat,1,function(x){!all(x==0)}), ]
group_norm <- as.character(group); group_norm[group_norm %in% c("gfp_for_egfr", "gfp18", "gfp30")] <- "gfp"
group_norm <- factor(group_norm, levels=c("gfp", "akt", "bad", "her2", "igf1r","raf", "egfr", "krasgv", "krasqh", "kraswt"))
start_time <- Sys.time()
combatseq_mat <- ComBat_seq(counts=cts_mat, batch=batch, group=group_norm, full_mod=TRUE, gene.subset.n=1000, Cpp=FALSE)
end_time <- Sys.time()
print(end_time - start_time)

