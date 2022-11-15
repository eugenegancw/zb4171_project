#load the library
library(GSEABenchmarkeR)
#load RNA seq data from package -> these are breast cancer and bladder cancer datasets
tcga <- loadEData("tcga", nr.datasets=2)

#run some DE method on dataset. In this case, we used DESeq2 
tcgadeseq <- runDE(tcga, de.method = "DESeq2",padj.method = "flexible")
#export DE results to csv files, sorted by adjusted Pvalue
write.csv(as.data.frame(rowData(tcgadeseq[[1]])[order(rowData(tcgadeseq[[1]])$ADJ.PVAL),]),file = "BLCAdeseq2sorted.csv")
write.csv(as.data.frame(rowData(tcgadeseq[[2]])[order(rowData(tcgadeseq[[2]])$ADJ.PVAL),]),file = "BRCAdeseq2sorted.csv")