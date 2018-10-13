## de_called1: edgeR directly used on raw counts
## de_called5: edgeR on ComBat-seq adjusted counts
sapply(c("ggplot2", "reshape2", "gridExtra"), require, character.only=TRUE)

setdiff(de_called5, de_called1)
setdiff(de_called1, de_called5)

length(setdiff(de_called5, de_called1)) / 100
length(setdiff(de_called1, de_called5)) / 818

# violin plots of counts in FP genes (removed by ComBat-seq) before vs after adjustment
fp_genes <- setdiff(de_called1, de_called5)
for(i in 1:length(fp_genes)){
  cts_df <- data.frame(Before=c(counts_matrix[fp_genes[i], group==0], counts_matrix[fp_genes[i], group==1]),
                       After=c(adj_counts_combatseq[fp_genes[i], group==0], adj_counts_combatseq[fp_genes[i], group==1]),
                       Group=as.factor(c(rep(0, length(which(group==0))), rep(1, length(which(group==1))))))
  cts_df_mlt <- melt(cts_df, id.vars="Group")
  
  png(sprintf("fp_genes/fp_%s.png", fp_genes[i]), width=6, height=6, units="in", res=300)
  print(ggplot(cts_df_mlt, aes(x=Group, y=value)) + 
          facet_grid(~variable) +
          geom_violin() +
          geom_boxplot(width=0.1) +
          labs(x="Group", y="Counts", title=fp_genes[i]))
  dev.off()
}


# violin plots of counts in true DE genes (first 100) before vs after adjustment
DE_genes <- paste0("gene",1:100)
for(j in 1:length(DE_genes)){
  cts_df <- data.frame(Before=c(counts_matrix[DE_genes[j], group==0], counts_matrix[DE_genes[j], group==1]),
                       After=c(adj_counts_combatseq[DE_genes[j], group==0], adj_counts_combatseq[DE_genes[j], group==1]),
                       Group=as.factor(c(rep(0, length(which(group==0))), rep(1, length(which(group==1))))))
  cts_df_mlt <- melt(cts_df, id.vars="Group")
  
  png(sprintf("DE_genes/DE_%s.png", DE_genes[j]), width=6, height=6, units="in", res=300)
  print(ggplot(cts_df_mlt, aes(x=Group, y=value)) + 
          facet_grid(~variable) +
          geom_violin() +
          geom_boxplot(width=0.1) +
          labs(x="Group", y="Counts", title=DE_genes[j]))
  dev.off()
}

