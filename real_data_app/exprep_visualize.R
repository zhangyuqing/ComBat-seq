rm(list=ls())
setwd("~/Google Drive/ComBat_seq/")
sapply(c("ggplot2", "gridExtra", "reshape2", "ggpubr"), require, character.only=TRUE)
load("real_data_app/perf_stats.RData")
alpha_fdr_sel <- 0.1

methods_vec <- levels(stats_out$FDR$Verification)
#i=1;j=1
sens_sel <- list()
for(i in seq_along(methods_vec)){  # for each verification method
  for(j in seq_along(methods_vec)){  # for each evaluation method
    # for ground truth generated using the current method (veri_method)
    curr_fdr_df <- stats_out$FDR[stats_out$FDR$Verification==methods_vec[i] & stats_out$FDR$Evaluation==methods_vec[j], ]
    fdr_avg <- aggregate(curr_fdr_df$Stats, by=list(FDR.cutoff=curr_fdr_df$Cutoff), mean, na.rm=TRUE)
    cutoff_to_use <- fdr_avg[which.min(abs(fdr_avg[,2]-alpha_fdr_sel)), 1]
    sens_sel <- c(sens_sel, list(stats_out$SENS[stats_out$SENS$Verification==methods_vec[i] & 
                                                  stats_out$SENS$Evaluation==methods_vec[j] & 
                                                  stats_out$SENS$Cutoff==cutoff_to_use, ]))
  }
}
sens_sel <- do.call(rbind, sens_sel)

sens_sel$Evaluation <- factor(sens_sel$Evaluation, levels=c("raw.count", "one.step", "curr.combat", "ruvseq", "svaseq", "combatseq"))
sens_sel$Verification <- factor(sens_sel$Verification, levels=c("raw.count", "one.step", "curr.combat", "ruvseq", "svaseq", "combatseq"))

png("realdata_app_TB.png", width=7, height=6, units="in", res=300)
plt <- ggplot(sens_sel, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  facet_wrap(~Verification, ncol=2) +
  labs(y="Sensitivity") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.title.x=element_blank())
print(plt)
dev.off()



# png("realdata_app.png", width=13, height=6.5, units="in", res=300)
# tpr_plt <- ggplot(stats_out$SENS, aes(x=Evaluation, y=Stats, color=Evaluation)) +
#   geom_boxplot() +
#   facet_wrap(~Verification, ncol=2) +
#   labs(y="Sensitivity") +
#   theme(axis.text.x=element_text(angle=45, hjust=1), 
#         axis.title.x=element_blank())
# fpr_plt <- ggplot(fpr_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
#   geom_boxplot() +
#   geom_hline(yintercept=alpha, color="blue", linetype="dashed") +
#   facet_wrap(~Verification, ncol=2) +
#   labs(y="FPR") +
#   theme(axis.text.x=element_text(angle=45, hjust=1), 
#         axis.title.x=element_blank())
# grid.arrange(tpr_plt, fpr_plt, ncol=2)
# dev.off()
