rm(list=ls())
setwd("~/Documents/ComBat_seq/DE_analysis_tmp/")
sapply(c("ggplot2", "reshape2", "gridExtra", "plyr"), require, character.only=TRUE)

load("DE_results.RData")
iters_sel <- 1:10 #1:length(predDEgenes_edgeR_unadj)
sel_method_list <- list(edgeR=c("BaseIndi.edgeR", #"BaseQuant.edgeR", 
                                "Batch.edgeR", "OneStep.edgeR", "ComBat.lm", 
                                "RUVseq.edgeR", "SVAseq.edgeR", "ComBatseq.edgeR"),
                        DESeq2=c("BaseIndi.DESeq2", #"BaseQuant.DESeq2", 
                                 "Batch.DESeq2", "OneStep.DESeq2", "ComBat.lm", 
                                 "RUVseq.DESeq2", "SVAseq.DESeq2", "ComBatseq.DESeq2"))

psig_fracs <- psig_numbers <- list()
for(iter in iters_sel){
  res_lst <- lapply(names(predDE_edgeR_resdfs[[iter]]), function(nm){
    de_res <- predDE_edgeR_resdfs[[iter]][[nm]]
    p_ups <- de_res[true_ups, "PValue"]
    p_downs <- de_res[true_downs, "PValue"]
    p_nulls <- de_res[true_nulls, "PValue"]
    
    p_df <- melt(list(All=de_res$PValue, Up=p_ups, Down=p_downs, Null=p_nulls))
    p_df$L1 <- factor(p_df$L1, levels=c("All", "Up", "Down", "Null"))
    plt <- ggplot(p_df, aes(x=value)) +
      facet_grid(~L1) +
      geom_histogram(breaks=seq(0,1,0.05)) +
      geom_vline(xintercept=0.05, color="blue", linetype="dashed") +
      theme(axis.text.x = element_text(angle=45, hjust=1)) + #,
            #panel.grid.major = element_blank(), 
            #panel.grid.minor = element_blank()) +
      #theme_linedraw() +
      labs(x="Uncorrected P values", y="Count", title=nm) 
    return(list(p=list(up=p_ups, down=p_downs, null=p_nulls), plt=plt))
  })
  names(res_lst) <- names(predDE_edgeR_resdfs[[iter]])
  
  # visualize distribution of (un-corrected) P values
  plt_lst <- lapply(res_lst, function(res){return(res$plt)})
  png(sprintf("orgPValDistribs_iter%s.png", iter), height=22, width=9, units="in", res=300)
  print(do.call(grid.arrange, c(plt_lst, ncol=1)))
  dev.off()
  
  # proportion of p values
  p_lst <- lapply(res_lst, function(res){return(res$p)})
  psig_fracs[[iter]] <- do.call(rbind, lapply(p_lst, function(curr_p_lst){
    sapply(curr_p_lst, function(pvec){sum(pvec < 0.05)/length(pvec)})
  }))
  psig_numbers[[iter]] <- psig_fracs * matrix(rep(c(length(true_ups), length(true_downs), length(true_nulls)), 
                                          nrow(psig_fracs)), ncol=3, byrow=TRUE)
}
