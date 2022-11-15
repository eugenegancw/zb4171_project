#load the library
library(GSEABenchmarkeR)

#import output files from EA tools as DataFrame objects
#Note: need to change some of the column header names
#notably the gene set id column to GENE.SET 
#and the qvalue/FDR/adjusted pvalue to PVAL
setwd("C:/Users/Shu Yuan/Documents/Shi SY/NUS/y3s1/zb4171/project/datafiles")

#import bladder cancer enrichment analysis output
blca_clusterProfiler_bp <- DataFrame(read.csv("cluseterProfiler_BP_custom_noThreshold.csv"))
blca_clusterProfiler_mf <- DataFrame(read.csv("cluseterProfiler_MF_custom_noThreshold.csv"))
blca_gProfiler_bp <- DataFrame(read.csv("gProfiler_bp_custom_noThreshold.csv"))
blca_gProfiler_mf <- DataFrame(read.csv("gProfiler_mf_custom_noThreshold.csv"))
blca_WebGestalt_bp <- DataFrame(read.csv("WebGestalt_BP_custom_noThreshold.csv"))
blca_WebGestalt_mf <- DataFrame(read.csv("WebGestalt_MF_custom_noThreshold.csv"))
blca_enrichr_bp <- DataFrame(read.csv("enrichr_go-bp.tsv",sep = "\t"))
blca_enrichr_mf <- DataFrame(read.csv("enrichr_go-mf.tsv",sep = "\t"))
blca_panther_bp <- DataFrame(read.csv("panther_go-bp.tsv",sep = "\t"))
blca_panther_mf <- DataFrame(read.csv("panther_go-mf.tsv",sep = "\t"))
blca_david_bp <- DataFrame(read.csv("david_blca_go-bp.tsv",sep = "\t"))
blca_david_mf <- DataFrame(read.csv("david_blca_go-mf.tsv",sep = "\t"))

#import breast cancer enrichment analysis output
brca_clusterProfiler_bp <- DataFrame(read.csv("brca_cluseterProfiler_BP_custom_noThreshold.csv"))
brca_clusterProfiler_mf <- DataFrame(read.csv("brca_cluseterProfiler_MF_custom_noThreshold.csv"))
brca_gProfiler_bp <- DataFrame(read.csv("brca_gProfiler_bp_custom_noThreshold.csv"))
brca_gProfiler_mf <- DataFrame(read.csv("brca_gProfiler_mf_custom_noThreshold.csv"))
brca_WebGestalt_bp <- DataFrame(read.csv("brca_WebGestalt_BP_custom_noThreshold.csv"))
brca_WebGestalt_mf <- DataFrame(read.csv("brca_WebGestalt_MF_custom_noThreshold.csv"))
brca_enrichr_bp <- DataFrame(read.csv("enrichr_brca_go-bp.tsv",sep = "\t"))
brca_enrichr_mf <- DataFrame(read.csv("enrichr_brca_go-mf.tsv",sep = "\t"))
brca_panther_bp <- DataFrame(read.csv("panther_brca_go-bp.tsv",sep = "\t"))
brca_panther_mf <- DataFrame(read.csv("panther_brca_go-mf.tsv",sep = "\t"))
brca_david_bp <- DataFrame(read.csv("david_brca_go-bp.tsv",sep = "\t"))
brca_david_mf <- DataFrame(read.csv("david_brca_go-mf.tsv",sep = "\t"))

#load reference relevance rankings from package files
data.dir <- system.file("extdata", package="GSEABenchmarkeR")
mala.go.file <- file.path(data.dir, "malacards","GO_BP.rds")
mala.go.bp <- readRDS(mala.go.file)
mala.go.mf.file <- file.path(data.dir, "malacards","GO_MF.rds")
mala.go.mf <- readRDS(mala.go.mf.file)

#obtain phenotypic enrichment scores

relevancescores_blca_bp <- c("clusterProfiler"=evalRelevance(blca_clusterProfiler_bp,mala.go.bp$BLCA),"gProfiler"=evalRelevance(blca_gProfiler_bp,mala.go.bp$BLCA), "WebGestalt" = evalRelevance(blca_WebGestalt_bp,mala.go.bp$BLCA),"Enrichr"=evalRelevance(blca_enrichr_bp,mala.go.bp$BLCA),"PANTHER"=evalRelevance(blca_panther_bp,mala.go.bp$BLCA),"DAVID"=evalRelevance(blca_david_bp,mala.go.bp$BLCA))
relevancescores_blca_mf <- c("clusterProfiler"=evalRelevance(blca_clusterProfiler_mf,mala.go.mf$BLCA),"gProfiler"=evalRelevance(blca_gProfiler_mf,mala.go.mf$BLCA), "WebGestalt" = evalRelevance(blca_WebGestalt_mf,mala.go.mf$BLCA),"Enrichr"=evalRelevance(blca_enrichr_mf,mala.go.mf$BLCA),"PANTHER"=evalRelevance(blca_panther_mf,mala.go.mf$BLCA),"DAVID"=evalRelevance(blca_david_mf,mala.go.mf$BLCA))
relevancescores_brca_bp <- c("clusterProfiler"=evalRelevance(brca_clusterProfiler_bp,mala.go.bp$BRCA),"gProfiler"=evalRelevance(brca_gProfiler_bp,mala.go.bp$BRCA), "WebGestalt" = evalRelevance(brca_WebGestalt_bp,mala.go.bp$BRCA),"Enrichr"=evalRelevance(brca_enrichr_bp,mala.go.bp$BRCA),"PANTHER"=evalRelevance(brca_panther_bp,mala.go.bp$BRCA),"DAVID"=evalRelevance(brca_david_bp,mala.go.bp$BRCA))
relevancescores_brca_mf <- c("clusterProfiler"=evalRelevance(brca_clusterProfiler_mf,mala.go.mf$BRCA),"gProfiler"=evalRelevance(brca_gProfiler_mf,mala.go.mf$BRCA), "WebGestalt" = evalRelevance(brca_WebGestalt_mf,mala.go.mf$BRCA),"Enrichr"=evalRelevance(brca_enrichr_mf,mala.go.mf$BRCA),"PANTHER"=evalRelevance(brca_panther_mf,mala.go.mf$BRCA),"DAVID"=evalRelevance(brca_david_mf,mala.go.mf$BRCA))

#obtain theoretical optimum relevance scores
optimumscores_blca_bp <- c("clusterProfiler"=compOpt(mala.go.bp$BLCA,gs.ids = blca_clusterProfiler_bp$GENE.SET),"gProfiler"=compOpt(mala.go.bp$BLCA,gs.ids = blca_gProfiler_bp$GENE.SET), "WebGestalt" = compOpt(mala.go.bp$BLCA,gs.ids = blca_WebGestalt_bp$GENE.SET),"Enrichr"=compOpt(mala.go.bp$BLCA,gs.ids = blca_enrichr_bp$GENE.SET),"PANTHER"=compOpt(mala.go.bp$BLCA,gs.ids = blca_panther_bp$GENE.SET),"DAVID"=compOpt(mala.go.bp$BLCA,gs.ids = blca_david_bp$GENE.SET))
optimumscores_blca_mf <- c("clusterProfiler"=compOpt(mala.go.mf$BLCA,gs.ids = blca_clusterProfiler_mf$GENE.SET),"gProfiler"=compOpt(mala.go.mf$BLCA,gs.ids = blca_gProfiler_mf$GENE.SET), "WebGestalt" = compOpt(mala.go.mf$BLCA,gs.ids = blca_WebGestalt_mf$GENE.SET),"Enrichr"=compOpt(mala.go.mf$BLCA,gs.ids = blca_enrichr_mf$GENE.SET),"PANTHER"=compOpt(mala.go.mf$BLCA,gs.ids = blca_panther_mf$GENE.SET),"DAVID"=compOpt(mala.go.mf$BLCA,gs.ids = blca_david_mf$GENE.SET))
optimumscores_brca_bp <- c("clusterProfiler"=compOpt(mala.go.bp$BRCA,gs.ids = brca_clusterProfiler_bp$GENE.SET),"gProfiler"=compOpt(mala.go.bp$BRCA,gs.ids = brca_gProfiler_bp$GENE.SET), "WebGestalt" = compOpt(mala.go.bp$BRCA,gs.ids = brca_WebGestalt_bp$GENE.SET),"Enrichr"=compOpt(mala.go.bp$BRCA,gs.ids = brca_enrichr_bp$GENE.SET),"PANTHER"=compOpt(mala.go.bp$BRCA,gs.ids = brca_panther_bp$GENE.SET),"DAVID"=compOpt(mala.go.bp$BRCA,gs.ids = brca_david_bp$GENE.SET))
optimumscores_brca_mf <- c("clusterProfiler"=compOpt(mala.go.mf$BRCA,gs.ids = brca_clusterProfiler_mf$GENE.SET),"gProfiler"=compOpt(mala.go.mf$BRCA,gs.ids = brca_gProfiler_mf$GENE.SET), "WebGestalt" = compOpt(mala.go.mf$BRCA,gs.ids = brca_WebGestalt_mf$GENE.SET),"Enrichr"=compOpt(mala.go.mf$BRCA,gs.ids = brca_enrichr_mf$GENE.SET),"PANTHER"=compOpt(mala.go.mf$BRCA,gs.ids = brca_panther_mf$GENE.SET),"DAVID"=compOpt(mala.go.mf$BRCA,gs.ids = brca_david_mf$GENE.SET))

#obtain normalized relevance scores
percent_optimum_relevance_scores_blca_bp <- 100*relevancescores_blca_bp/optimumscores_blca_bp
percent_optimum_relevance_scores_blca_mf <- 100*relevancescores_blca_mf/optimumscores_blca_mf
percent_optimum_relevance_scores_brca_bp <- 100*relevancescores_brca_bp/optimumscores_brca_bp
percent_optimum_relevance_scores_brca_mf <- 100*relevancescores_brca_mf/optimumscores_brca_mf

#plot results
boxplot(c(percent_optimum_relevance_scores_blca_bp[[1]],percent_optimum_relevance_scores_blca_mf[[1]],percent_optimum_relevance_scores_brca_bp[[1]],percent_optimum_relevance_scores_brca_mf[[1]]),c(percent_optimum_relevance_scores_blca_bp[[2]],percent_optimum_relevance_scores_blca_mf[[2]],percent_optimum_relevance_scores_brca_bp[[2]],percent_optimum_relevance_scores_brca_mf[[2]]),c(percent_optimum_relevance_scores_blca_bp[[3]],percent_optimum_relevance_scores_blca_mf[[3]],percent_optimum_relevance_scores_brca_bp[[3]],percent_optimum_relevance_scores_brca_mf[[3]]),c(percent_optimum_relevance_scores_blca_bp[[4]],percent_optimum_relevance_scores_blca_mf[[4]],percent_optimum_relevance_scores_brca_bp[[4]],percent_optimum_relevance_scores_brca_mf[[4]]),c(percent_optimum_relevance_scores_blca_bp[[5]],percent_optimum_relevance_scores_blca_mf[[5]],percent_optimum_relevance_scores_brca_bp[[5]],percent_optimum_relevance_scores_brca_mf[[5]]),c(percent_optimum_relevance_scores_blca_bp[[6]],percent_optimum_relevance_scores_blca_mf[[6]],percent_optimum_relevance_scores_brca_bp[[6]],percent_optimum_relevance_scores_brca_mf[[6]]),names = c("clusterProfiler","gProfiler","WebGestalt","Enrichr","PANTHER","DAVID"),xlab = "enrichment analysis tools",ylab = "%optimum relevance score", main = "Comparison of Gene set enrichment analysis tools")
