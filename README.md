Codes in this repo has been written by Eugene Gan Cheong Wei, Keith Chew Zikai and Shi Shu Yuan for ZB4171 - Pathway enrichment analysis tool comparison.
## **Obtain tcga data set, specifically BLCA & BRCA**
\- The analysis and its parameters are found under zb4171\_project/obtain\_tcgadata.R
## **Pre-processing for quantitative analysis**
\- The analysis and its parameters are found under zb4171\_project/geneListFilter/filterGL.ipynb

\-  The filtered gene lists are found under zb4171\_project/geneListFilter/filtered folder
## **Tools enrichment analysis**
All the outputs, including clusterProfiler, are found under zb4171\_project/Tools\_output folder

[**WebGestalt**](http://www.webgestalt.org/)

**Basic parameters:**

1. Use the respective custom gene set for “Functional Database” 
1. Copy & paste the filtered gene list for Gene List section
1. Select "genome" for Reference Set

**Advanced parameters:**

1. Select "FDR" and input 1 for significance level 

Gene set annotation used: GO:BP 2021, GO:MF 2021 (custom)

[**g:Profiler**](https://biit.cs.ut.ee/gprofiler/gost)

\- Use “g:GOSt”

\- Copy & paste the filtered gene list for query section

**Options:**

1. Tick "ordered query"

**Advanced options:**

1. Tick "All results"

**Data sources:**

1. Click “clear all”

**Bring your data (Custom GMT):**

1. Use the respective custom gene set for Functional Database

Gene set annotation used: GO:BP 2021, GO:MF 2021 (custom)

[**DAVID**](https://david.ncifcrf.gov/summary.jsp)

\- Use “Functional Annotation Tool”

**Step 1A:**

1. Paste the filtered gene list into box

**Step 2:**

1. Select “OFFICIAL\_GENE\_SYMBOL” 
1. Input “9606” or “Homo Sapiens” for Step 2a

**Step 3:**

1. tick “Gene List”

**Step 4:**

1. Click “Submit List” 

Gene set annotation used:  GO:BP 2016, GO:MF 2016

**ClusterProfiler**

\- The analysis and its parameters is found under zb4171\_project/clusterProfilerCode/clusterProfiler.R

Gene set annotation used: GO:BP 2021, GO:MF 2021 (custom)

[**PANTHER**](http://www.pantherdb.org/)

**Step 1:**

1. paste the filtered gene list in box 
1. select list type as “ID List”

**Step 2:**

1. select organism “Homo sapiens”

**Step 3:**

1. select “Statistical overrepresentation test” for Step 3
1. select “GO biological process / molecular function complete” for Step 3

\- Click “submit”

**Default whole-genome lists:**

1. select “Homo sapiens”

\- Click “Upload list”

\- Tick “Binomial” and “Calculate False Discovery Rate”

\- Click “Launch analysis”

\- After it launches and shows the results for FDR<0.05, change the setting to “show all results”.

Gene set annotation used: GO:BP 2022, GO:MF 2022

[**Enrichr**](https://maayanlab.cloud/Enrichr/)

\- Paste the filtered gene list in the box on the right

\- Click “Submit”

\- Select “Ontologies”

\- Select “GO Biological Process / Molecular Function 2021”

\- Select “Table” and download the data

\- Gene set annotation used: GO:BP 2021, GO:MF 2021
## **Quantitative analysis**
\- The analysis and its parameters are found under zb4171\_project/run\_gseabenchmarkeR.R


