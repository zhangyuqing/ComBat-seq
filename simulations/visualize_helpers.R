ObsNominalFDR <- function(prec_df){   
  avg_fdr_obs <- aggregate(1-prec_df[,-1], by=list(FDR.cutoff=prec_df$FDR.cutoff), mean, na.rm=TRUE)
  avg_fdr_obs_df <- melt(avg_fdr_obs, id.vars="FDR.cutoff")
  return(list(fdr_df=avg_fdr_obs, fdr_df_mlt=avg_fdr_obs_df))
}


PrecSensMapping <- function(prec_df, sens_df, alpha_fdr_sel){
  avg_prec_obs <- aggregate(prec_df[,-1], by=list(FDR.cutoff=prec_df$FDR.cutoff), mean, na.rm=TRUE)
  
  methods_vec <- colnames(prec_df)[-1]
  sens_sub <- list()
  for(ii in seq_along(methods_vec)){
    cutoff_to_use <- avg_prec_obs[which.min(abs(avg_prec_obs[, methods_vec[ii]]-(1-alpha_fdr_sel))), "FDR.cutoff"]
    sens_sub[[ii]] <- sens_df[sens_df$FDR.cutoff==cutoff_to_use, methods_vec[ii]]
  }
  names(sens_sub) <- methods_vec
  return(sens_sub)
}


CleanSensOut <- function(sens_out, DE_method){
  sens_out_clean <- melt(sens_out) 
  colnames(sens_out_clean) <- c("Sensitivity", "Method", "Confounding")
  sens_out_clean$Method <- gsub(DE_method, "", sens_out_clean$Method)
  sens_out_clean$Method <- factor(sens_out_clean$Method, levels=unique(sens_out_clean$Method))
  sens_out_clean$Confounding <- factor(sens_out_clean$Confounding, levels=unique(sens_out_clean$Confounding))
  return(sens_out_clean)
}

CleanOut <- function(out_lst, stat_name, DE_method){
  out_merged <- melt(out_lst); 
  colnames(out_merged) <- c(stat_name, "Method", "Disp", "Mean")
  out_merged$Method <- gsub(DE_method, "", out_merged$Method)
  out_merged$Method <- factor(out_merged$Method, levels=unique(out_merged$Method))
  out_merged$Mean <- factor(out_merged$Mean, levels=unique(out_merged$Mean))
  return(out_merged)
}


CleanOut_Cnfnd <- function(out_lst, stat_name, DE_method){
  out_merged <- melt(out_lst); 
  colnames(out_merged) <- c(stat_name, "Method", "Confounding")
  out_merged$Method <- gsub(DE_method, "", out_merged$Method)
  out_merged$Method <- factor(out_merged$Method, levels=unique(out_merged$Method))
  out_merged$Confounding <- factor(out_merged$Confounding, levels=unique(out_merged$Confounding))
  return(out_merged)
}

