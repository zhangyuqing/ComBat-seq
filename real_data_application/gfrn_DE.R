rm(list=ls())
sapply(c("sva", "DESeq2", "ggplot2", "reshape2", "gridExtra", "scales", "dendextend", 
         "RUVSeq", "ggpubr", "BatchQC", "clusterProfiler", "dplyr", "org.Hs.eg.db", "DOSE", "fgsea"), 
       require, character.only=TRUE)

## Parameters (change paths when necessary)
data_dir <- "~/Desktop/ComBat-seq/real_data_application"  # path to the signature data (.rds)
source("~/Desktop/ComBat-seq/real_data_application/gfrn_helpers.R")  # path to gfrn_helpers.R
source("~/Desktop/ComBat-seq/ComBat_seq.R"); source("~/Desktop/ComBat-seq/helper_seq.R")   
output_dir <- "~/Desktop/ComBat-seq/real_data_application"

pathway_regex <- c("her2", "^egfr", "kraswt")  #"^egfr" 
set.seed(1)



###########  Load data
sigdata <- readRDS(file.path(data_dir, "signature_data.rds"))
cts_mat <- assay(sigdata, "counts")   # count matrix 
batch <- colData(sigdata)$batch
group <- colData(sigdata)$group

# Take the subset (all controls & 1 condition specified by pathway_regex)
pathway_condition_ind <- grep(paste(pathway_regex, collapse="|"), group)
#pathway_batch <- unique(batch[pathway_condition_ind])  # which batch does this pathway belong to
ctrl_ind <- which(group %in% c("gfp_for_egfr", "gfp18", "gfp30"))
subset_ind <- sort(c(ctrl_ind, pathway_condition_ind))  
cts_sub <- cts_mat[, subset_ind]
batch_sub <- batch[subset_ind]
group_sub <- group[subset_ind]
# remove genes with only 0 counts in the subset & in any batch
#cts_sub <- cts_sub[apply(cts_sub,1,function(x){!all(x==0)}), ]
keep1 <- apply(cts_sub[, batch_sub==1],1,function(x){!all(x==0)})
keep2 <- apply(cts_sub[, batch_sub==2],1,function(x){!all(x==0)})
keep3 <- apply(cts_sub[, batch_sub==3],1,function(x){!all(x==0)})
cts_sub <- cts_sub[keep1 & keep2 & keep3, ]
dim(cts_sub)



###########  Batch correction
group_sub <- factor(as.character(group_sub), levels=c("gfp_for_egfr", "gfp18", "gfp30",  gsub("^", "", pathway_regex, fixed=T)))
group_sub <- plyr::revalue(group_sub, c("gfp_for_egfr"="gfp", "gfp18"="gfp", "gfp30"="gfp"))

## Use ComBatSeq to adjust data
combatseq_sub <- ComBat_seq(counts=cts_sub, batch=batch_sub, group=group_sub,
                            shrink=FALSE, shrink.disp=FALSE)

## Use original ComBat on logCPM
combat_sub <- ComBat(cpm(cts_sub, log=TRUE), batch=batch_sub, mod=model.matrix(~group_sub))

## RUVseq
group1 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==1])), c("gfp"="0", "her2"="1"))
deres1 <- edgeR_DEpipe(cts_sub[, batch_sub==1], batch=NULL, group=group1,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
group2 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==2])), c("gfp"="0", "egfr"="1"))
deres2 <- edgeR_DEpipe(cts_sub[, batch_sub==2], batch=NULL, group=group2,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
group3 <- plyr::revalue(as.factor(as.character(group_sub[batch_sub==3])), c("gfp"="0", "kraswt"="1"))
deres3 <- edgeR_DEpipe(cts_sub[, batch_sub==3], batch=NULL, group=group3,
                       include.batch=FALSE, alpha.unadj=1, alpha.fdr=1)
null_genes <- Reduce(intersect, list(which(deres1$de_res$FDR>0.95), which(deres2$de_res$FDR>0.95), which(deres3$de_res$FDR>0.95)))
group_obj <- makeGroups(group_sub)
ruvseq_sub <- RUVs(cts_sub, cIdx=null_genes, scIdx=group_obj, k=1)$normalizedCounts


## SVAseq
mod1 <- model.matrix(~group_sub); mod0 <- cbind(mod1[,1])
num.sv(cts_sub, mod=mod1)
svseq <- svaseq(cts_sub, mod=mod1, mod0=mod0, n.sv=1); cat("\n")



###########  Look at oncogene expression
path_mapping <- c("her2", "egfr", "kraswt")
names(path_mapping) <- c("ERBB2", "EGFR", "KRAS")

curr_path <- "EGFR"
for(curr_path in names(path_mapping)){
  curr_sample_ind <- which(group_sub %in% c("gfp", path_mapping[curr_path]))
  group_sub_rev <- plyr::revalue(group_sub, 
                                 c("gfp"="GFP (control)", "her2"="HER2", "egfr"="EGFR", "kraswt"="KRAS"))
  
  ## transform to logCPM
  cpm_df <- data.frame(Unadjusted=cpm(cts_sub, log=TRUE)[curr_path, curr_sample_ind],
                       `Original ComBat`=combat_sub[curr_path, curr_sample_ind],
                       `RUV-Seq`=cpm(ruvseq_sub, log=TRUE)[curr_path, curr_sample_ind],
                       `ComBat-Seq`=cpm(combatseq_sub, log=TRUE)[curr_path, curr_sample_ind],
                       Batch=as.factor(batch_sub[curr_sample_ind]),
                       Condition=as.character(group_sub_rev)[curr_sample_ind])
  cpm_df_mlt <- melt(cpm_df, variable.name="Method", id.vars=c("Batch", "Condition"))
  cpm_df_mlt$Condition <- factor(cpm_df_mlt$Condition, levels=c("GFP (control)", 
                                                                setdiff(levels(cpm_df_mlt$Condition), "GFP (control)")))
  
  png(file.path(output_dir, sprintf("gfrnEXPR_logCPM_%s.png", curr_path)), 
      width=8, height=5, units="in", res=300)
  print(ggplot(cpm_df_mlt, aes(x=Condition, y=value)) +
          facet_wrap(~Method, ncol=4)+ #, scales="free_y") +
          geom_boxplot(outlier.shape=NA) +
          geom_jitter(aes(color=Batch)) +
          labs(y=sprintf("%s Expression (logCPM)", curr_path)))
  dev.off()
}



###########  DE & enrichment analysis 
path_mapping <- c("her2", "egfr", "kraswt")
names(path_mapping) <- c("ERBB2", "EGFR", "KRAS")
#pathways.msigdb <- gmtPathways(file.path(data_dir, "msigdb.v7.1.symbols.gmt"))

curr_path <- "her2"
out_allgenes <- out_DEgenes <- num_DEgenes <- list()
for(curr_path in path_mapping){
  ####  Compare each condition with all controls pooled across batch
  ##  Subset of this pathway
  curr_sample_ind <- which(group_sub %in% c("gfp", curr_path))
  
  curr_batch <- batch_sub[curr_sample_ind]
  curr_group <- factor(as.character(group_sub[curr_sample_ind]), levels=c("gfp", curr_path))
  print(table(curr_batch))
  print(table(curr_group))
  
  curr_cts <- cts_sub[, curr_sample_ind]
  curr_combatseq <- combatseq_sub[, curr_sample_ind]
  curr_combat <- combat_sub[, curr_sample_ind]
  curr_ruvseq <- ruvseq_sub[, curr_sample_ind]
  curr_sv <- matrix(svseq$sv, ncol=1)[curr_sample_ind, ]
  
  ##  DE
  res_unadjusted <- edgeR_DEpipe(cts=curr_cts, batch=curr_batch, group=curr_group, 
                                 include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05) 
  res_onestep <- edgeR_DEpipe(cts=curr_cts, batch=curr_batch, group=curr_group, 
                              include.batch=TRUE, alpha.unadj=0.05, alpha.fdr=0.05) 
  res_combatseq <- edgeR_DEpipe(cts=curr_combatseq, batch=curr_batch, group=curr_group, 
                                include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05) 
  res_ruvseq <- edgeR_DEpipe(cts=curr_ruvseq, batch=curr_batch, group=curr_group, 
                             include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05) 
  res_svaseq <- edgeR_DEpipe(cts=curr_cts, batch=curr_batch, group=curr_group, covar=curr_sv, 
                             include.batch=FALSE, alpha.unadj=0.05, alpha.fdr=0.05)  
  fit <- lmFit(curr_combat, design=model.matrix(~curr_group))
  fit <- eBayes(fit)
  res_combat <- list(de_res=topTable(fit, n=nrow(curr_combat)), 
                     design=model.matrix(~curr_group))
  
  DE_objs <- list(Unadjusted=res_unadjusted, 
                  OneStep=res_onestep,
                  ComBat=res_combat, 
                  RUVseq=res_ruvseq,
                  SVAseq=res_svaseq,
                  ComBatseq=res_combatseq)
  DE_tables <- lapply(DE_objs, function(de_obj){de_obj$de_res})
  #tb_nm=names(DE_tables)[1]
  DEgene_list <- lapply(names(DE_tables), function(tb_nm){
    curr_tb <- DE_tables[[tb_nm]]
    # output DE gene table
    write.csv(curr_tb[1:100, c(1,3,4,5)], row.names=TRUE, quote=FALSE,
              file=file.path(output_dir, sprintf(sprintf("DEgenes_%s_%s.csv", curr_path, tb_nm))))

    if(tb_nm=="ComBat"){col_nm <- "adj.P.Val"}else{col_nm <- "FDR"}
    curr_tb_sig <- curr_tb[curr_tb[, col_nm] < 0.05, ]
    curr_gene_list <- curr_tb_sig$logFC
    names(curr_gene_list) <- rownames(curr_tb_sig)
    
    # gsea
    # png(file.path(output_dir, sprintf(sprintf("GSEA_%s_%s.png", curr_path, tb_nm))),
    #     width=10, height=7, units="in", res=300)
    # print(plotEnrichment(pathway=pathways.msigdb[["PID_RAS_PATHWAY"]], stats=curr_gene_list) +
    #         labs(title=sprintf("RAS PATHWAY, %s %s", curr_path, tb_nm)))
    # dev.off()
    
    # curr_tb_sig <- curr_tb[curr_tb[, col_nm] < 0.05, ]
    # curr_gene_list <- curr_tb_sig$logFC
    # curr_tb_sig <- curr_tb_sig %>% arrange(desc(abs(logFC)))
    # names(curr_gene_list) <- rownames(curr_tb_sig)
    # curr_gene_list <- sort(curr_gene_list, decreasing=TRUE)
    # gse <- gseGO(geneList=curr_gene_list, ont="ALL", keyType="SYMBOL", 
    #              OrgDb="org.Hs.eg.db", pvalueCutoff=0.2)
    # if(dim(gse)[1]!=0){
    #   png(file.path(output_dir, sprintf(sprintf("Enrich_%s_%s.png", curr_path, tb_nm))),
    #       width=10, height=7, units="in", res=300)
    #   print(dotplot(gse, showCategory=10, split=".sign") + 
    #     facet_grid(.~.sign))
    #   dev.off()
    # }
    
    return(curr_gene_list)
  })
  
  names(DEgene_list) <- names(DE_tables)
  num_DEgenes[[curr_path]] <- sapply(DEgene_list, length)
  out_DEgenes[[curr_path]] <- DEgene_list
  out_allgenes[[curr_path]] <- DE_tables
}
#save(out_DEgenes, out_allgenes, file="GFRN_DE.RData")


## Number of DE genes by each method for each pathway
do.call(rbind, num_DEgenes)


## Overlap between DE genes
degenes <- lapply(out_DEgenes$her2, function(x){names(x)})
ll <- combn(degenes, 2, simplify=FALSE)
out <- lapply(ll, function(x) length(intersect(x[[1]], x[[2]])))
do.call(c, out)
names(degenes)
sapply(degenes, length)


## Ranking and FDR of EGFR
sapply(out_DEgenes$egfr, function(x){
  which(names(x)=="EGFR")/length(x)
})
lapply(out_allgenes$egfr, function(x){x["EGFR",]})


## RAS pathway gene rankings
ras_pathway_genes <- read.csv(file.path(data_dir, "ras-pathway-gene-names.csv"), as.is=TRUE)[,1]
i=1;j=1
out2 <- list()
for(i in seq_along(out_allgenes)){
  out1 <- c()
  for(j in seq_along(out_allgenes[[i]])){
    curr_tb <- out_allgenes[[i]][[j]]
    overlap_genes <- intersect(rownames(curr_tb), ras_pathway_genes)
    out1[j] <- sum(which(rownames(curr_tb) %in% overlap_genes))
  }
  names(out1) <- names(out_allgenes[[i]])
  out2[[i]] <- out1
}
names(out2) <- names(out_allgenes)
apply(do.call(rbind, out2), 1, function(x){colnames(do.call(rbind, out2))[order(x)]})

curr_path=1; m=1
ras_in_dt <- intersect(rownames(out_allgenes[[1]][[1]]), ras_pathway_genes)
pcnt_cmp <- list()
for(curr_path in seq_along(out_allgenes)){
  curr_res <- out_allgenes[[curr_path]]
  pcnt_lst <- lapply(1:length(curr_res), function(m){
    curr_gene_tb <- curr_res[[m]] 
    percent_in_ras <- sapply(1:nrow(curr_gene_tb), function(r){
      length(intersect(rownames(curr_gene_tb)[1:r], ras_in_dt)) / length(ras_in_dt)
    })                                                         
  })
  names(pcnt_lst) <- names(curr_res)
  pcnt_dt <- melt(pcnt_lst)
  colnames(pcnt_dt)[2] <- c("Method")
  pcnt_dt$Method <- factor(pcnt_dt$Method, levels=c("Unadjusted", "OneStep", "ComBat", "SVAseq", "RUVseq", "ComBatseq"))
  pcnt_cmp[[curr_path]] <- do.call(rbind, pcnt_lst)
}

print(pcnt_cmp[[3]][order(pcnt_cmp[[3]][, 1000], decreasing=TRUE), 1000] * length(ras_in_dt))

# Fisher exact test
ras_res <- pcnt_cmp[[3]][order(pcnt_cmp[[3]][, 1000], decreasing=TRUE), 1000] * length(ras_in_dt)
sapply(ras_res, function(x){
  tb <- matrix(c(1000, nrow(cts_sub), x, length(ras_in_dt)), nrow=2, ncol=2)
  fisher.test(tb)$p.value
})

