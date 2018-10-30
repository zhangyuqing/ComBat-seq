plotPerf <- function(tpr_res, fpr_res, fpr_cutoff=0.05){
  if(!identical(dim(tpr_res), dim(fpr_res))){stop("TPR and FPR result file should have matching dimensions!")}
  tpr_res <- as.matrix(tpr_res); fpr_res <- as.matrix(fpr_res)
  
  # only keep TPR where FPR is within the given cutoff
  fpr_mask <- fpr_res < fpr_cutoff
  tpr_res[!fpr_mask] <- NA
  
  tpr_lst <- as.list(as.data.frame(tpr_res))
  boxplot(tpr_lst)
}
