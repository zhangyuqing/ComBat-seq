rm(list=ls())
setwd("")
sapply(c("ggplot2", ""), require, character.only=TRUE)
load("perf_stats.RData")


png("perf_stats.png", width=13, height=6.5, units="in", res=300)
tpr_plt <- ggplot(tpr_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  facet_wrap(~Verification, ncol=2) +
  labs(y="TPR") +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.title.x=element_blank())
fpr_plt <- ggplot(fpr_mat, aes(x=Evaluation, y=Stats, color=Evaluation)) +
  geom_boxplot() +
  geom_hline(yintercept=alpha, color="blue", linetype="dashed") +
  facet_wrap(~Verification, ncol=2) +
  labs(y="FPR") +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.title.x=element_blank())
grid.arrange(tpr_plt, fpr_plt, ncol=2)
dev.off()
