## Example 1
y1_ctrl <- matrix(rnbinom(1000*5,mu=10,size=10),ncol=5); mean(c(y1_ctrl)); var(c(y1_ctrl))
y1_case <- matrix(rnbinom(1000*5,mu=15,size=10),ncol=5); mean(c(y1_case)); var(c(y1_case))
y2_ctrl <- matrix(rnbinom(1000*3,mu=15,size=5), ncol=3); mean(c(y2_ctrl)); var(c(y2_ctrl))
y2_case <- matrix(rnbinom(1000*3,mu=20,size=5), ncol=3); mean(c(y2_case)); var(c(y2_case))
y3_ctrl <- matrix(rnbinom(1000*1,mu=5, size=20),ncol=1); mean(c(y3_ctrl)); var(c(y3_ctrl))
y3_case <- matrix(rnbinom(1000*1,mu=10,size=20),ncol=1); mean(c(y3_case)); var(c(y3_case))
counts <- cbind(y1_ctrl, y1_case, y2_ctrl, y2_case, y3_ctrl, y3_case)
batch <- c(rep("B", ncol(y1_ctrl)+ncol(y1_case)), 
           rep("A", ncol(y2_ctrl)+ncol(y2_case)), 
           rep("C", ncol(y3_ctrl)+ncol(y3_case))); table(batch)
group <- c(rep(0, ncol(y1_ctrl)), rep(1, ncol(y1_case)),
           rep(0, ncol(y2_ctrl)), rep(1, ncol(y2_case)),
           rep(0, ncol(y3_ctrl)), rep(1, ncol(y3_case)))
full_mod <- TRUE
source("helper_seq.R")



## Example 2 - a simpler example with equal sample size
N_genes <- 6
N_sample_vec <- rep(5, 4)
N_batch <- 2

lambda_1_ctrl <- matrix(rgamma(N_genes*N_sample_vec[1], shape=0.5, rate=1), nrow=N_genes)
lambda_1_case <- matrix(rgamma(N_genes*N_sample_vec[2], shape=1, rate=2), nrow=N_genes)
full_mod <- TRUE
source("helper_seq.R")
