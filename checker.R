# design matrix
cat("\n########  Design matrix  ########\n")
cat("Partial design matrix:\n"); print(head(design))
# dispersion
cat("\n\n########  Dispersion  ########\n")
cat("Estimated common dispersion:\n"); print(round(disp_common,3))
cat("\nAverage estimated gene-wise dispersion:\n"); print(round(sapply(genewise_disp_lst, mean), 3))
# GLM model
cat("\n\n########  GLM model coefs  ########\n")
cat("Coefficients (model 1):\n"); print(head(glm_f$coefficients)); cat('...\n'); print(tail(glm_f$coefficients))
cat("\nSanity check: estimated batch mean in gene group 1:\n"); print(exp(colMeans(glm_f$coefficients[1:1000,1:n_batch]))*mean(colSums(counts)))
cat("\nAverage background count alpha:\n"); print(mean(alpha_g)); print(exp(mean(alpha_g))*mean(colSums(counts)))
cat("\nCoefficients (model 2):\n"); print(head(glm_f2$coefficients)); cat('...\n'); print(tail(glm_f2$coefficients))
# Posterior estimates of batch parameters
cat("\n\n########  Batch effect estimates  ########\n")
cat("Prior average exp(gamma) in gene group 1:\n")
cat("NOTE: adjusted mean = original mean / exp(gamma) \n")
print(exp(colMeans(gamma_hat[1:1000,])))
cat("\nPosterior average exp(gamma) in gene group 1:\n") 
print(exp(colMeans(gamma_star_mat[1:1000,])))
cat("\nPrior average dispersion:\n")
print(colMeans(phi_hat))
cat("\nPosterior average dispersion:\n")
print(colMeans(phi_star_mat))
cat("\nBatch-free mean mu:\n")
print(mean(mu_star[1001:2000,batch=="A"])); print(mean(mu_star[1001:2000,batch=="B"]))
cat("\nBatch-free dispersion phi:\n")
print(mean(phi_star))

