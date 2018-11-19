ObsNominalFDR <- function(prec_res, plt_title){
  avg_fdr_obs <- aggregate(1-prec_res[,-1], by=list(FDR.cutoff=prec_res$FDR.cutoff), mean, na.rm=TRUE)
  avg_fdr_obs_df <- melt(avg_fdr_obs, id.vars="FDR.cutoff")
  avg_fdr_obs_df$FDR.cutoff <- as.numeric(as.character(avg_fdr_obs_df$FDR.cutoff))
  plt <- ggplot(avg_fdr_obs_df, aes(x=FDR.cutoff, y=value, group=variable, color=variable)) +
    geom_line() +
    geom_abline(slope=1, intercept=0, color="black") +
    labs(x="evaluation set FDR cutoff", y="observed FDR", title=plt_title) +
    scale_x_continuous(limits=c(0, 0.2)) +
    scale_y_continuous(limits=c(0, 0.3)) +
    theme_bw()
  return(list(plt=plt, fdr_df=avg_fdr_obs))
}


#prec_res=prec_res_edgeR; sens_res=sens_res_edgeR
PrecSensMapping <- function(prec_res, sens_res, alpha_fdr_sel){
  avg_prec_obs <- aggregate(prec_res[,-1], by=list(FDR.cutoff=prec_res$FDR.cutoff), mean, na.rm=TRUE)
  methods_vec <- colnames(prec_res)[-1]
  sens_sub <- list()
  for(i in seq_along(methods_vec)){
    cutoff_to_use <- avg_prec_obs[which.min(abs(avg_prec_obs[, methods_vec[i]]-(1-alpha_fdr_sel))), "FDR.cutoff"]
    sens_sub[[i]] <- sens_res[sens_res$FDR.cutoff==cutoff_to_use, methods_vec[i]]
  }
  names(sens_sub) <- methods_vec; sens_sub <- do.call(cbind, sens_sub)
  return(sens_sub)
}





