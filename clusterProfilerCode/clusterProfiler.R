library(tidyverse)
library(GO.db)
library(clusterProfiler)
library(fgsea)
library(msigdbr)

## BRCA dataset
#load data
model.results<- read_csv("../geneListFilter/filtered_BRCA_DEseq_list.csv")
#Define significant gene
signif<-model.results
signif.entrez<-unique(signif$SYMBOL)

#Analyse using custom_BP
custom_BP <-read.gmt("../data/go-bp_2021_gene_set.gmt")
custom_BP.gene_symbol<- select(custom_BP, term, gene)
enrich.GO<- enricher(gene=signif.entrez, pvalueCutoff=1, TERM2GENE = custom_BP.gene_symbol)
head(enrich.GO@result)
write.csv(enrich.GO@result,"brca_clusterProfiler_output/brca_cluseterProfiler_BP_custom_noThreshold.csv", row.names = FALSE)

#Analyse using custom_MF
custom_MF <-read.gmt("../data/go-mf_2021_gene_set.gmt")
custom_MF.gene_symbol<- select(custom_MF, term, gene)
#note that for our project, we output geneset lists regardless of significance level
enrich.GO<- enricher(gene=signif.entrez, pvalueCutoff=1, TERM2GENE = custom_MF.gene_symbol)
head(enrich.GO@result)
write.csv(enrich.GO@result,"brca_clusterProfiler_output/brca_cluseterProfiler_MF_custom_noThreshold.csv", row.names = FALSE)

## BLCA dataset
#load data
model.results<- read_csv("../geneListFilter/filtered_BLCA_DEseq_list.csv")
#Define significant gene
signif_l<-model.results
signif_l.entrez<-unique(signif_l$SYMBOL)

#Analyse using custom_BP
enrich.GO<- enricher(gene=signif_l.entrez, pvalueCutoff=1, TERM2GENE = custom_BP.gene_symbol)
head(enrich.GO@result)
write.csv(enrich.GO@result,"blca_clusterProfiler_output/blca_cluseterProfiler_BP_custom_noThreshold.csv", row.names = FALSE)

#Analyse using custom_MF
#note that for our project, we output geneset lists regardless of significance level
enrich.GO<- enricher(gene=signif_l.entrez, pvalueCutoff=1, TERM2GENE = custom_MF.gene_symbol)
head(enrich.GO@result)
write.csv(enrich.GO@result,"blca_clusterProfiler_output/blca_cluseterProfiler_MF_custom_noThreshold.csv", row.names = FALSE)
