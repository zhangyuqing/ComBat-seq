rm(list=ls())
setwd("~/Documents/ComBat_seq/DE_analysis_tmp/")
sapply(c("ggplot2", "gridExtra", "reshape2"), require, character.only=TRUE)

tpr_res <- read.csv("tpr_simFC_bio2_batch2_sizes100_10_N20_B_depth5.csv")
fpr_res <- read.csv("fpr_simFC_bio2_batch2_sizes100_10_N20_B_depth5.csv")

p_tpr <- ggplot(melt(tpr_res), aes(x=variable, y=value)) +
  geom_boxplot() +
  labs(x="", y="TPR", title="Power") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
p_fpr <- ggplot(melt(fpr_res), aes(x=variable, y=value)) +
  geom_boxplot() +
  geom_hline(yintercept=0.05, col="red") +
  labs(x="", y="FPR", title="Type-I error rate") +
  theme(axis.text.x=element_text(angle=45, hjust=1))

png("sim_zero.png", width=7, height=7, units="in", res=300)
grid.arrange(p_fpr, p_tpr, ncol=2)
dev.off()
