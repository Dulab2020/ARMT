# ARMT
## Introduction:
ARMT is auto RNA-seq data minning tool.<br>
The aim of ARMT is : <br>
    i) integrate and facilitate the downstream analysis of RNA-seq data,<br>
    ii) provide the way to analyze geneset according to GSVA, <br>
    iii) explore the correlation of gene expression and mutation.<br>
Author:<br>
Guanda Huang —— 202010108432@mail.scut.edu.cn<br>
Hongli Du —— hldu@scut.edu.cn<br>
Year: 2021<br>
License: GPL(>=3)<br>
#### The packages from Bioconductor including:
maftools<br>
GSVA<br>
GSVAdata<br>
limma<br>
GSEAbase<br>
org.Hs.eg.db<br>
edgeR<br>
survival<br>
clusterProfiler<br>
DOSE<br>
## Function:
Provide TCGA clinical data<br>
Create .gmt file for arbitary geneset<br>
Normalization of counts matrix (TPM)<br>
Gene set variant analysis<br>
Survival analysis<br>
Differential analysis<br>
Correlation analysis<br>
Analysis for multiple sets of data<br>
Plot the mutation information<br>
## Feature:
Easy to use<br>
Automation<br>
Visualization<br>
Comprehensiveness<br>
Integrating GSVA score analysis<br>
## Structure:
![Figure0](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure0.png)
## Dependencies:
R >=4.0.2, Rstudio, R packages(devtools or remotes)<br>
The other dependencies are automatically installed. By default, outdated dependencies are automatically upgraded. In interactive sessions you can select a subset of the dependencies to upgrade.
## Installation:
To install this package from Github, please, use the code below.
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("Dulab2020/ARMT")
```
## Quick Star:
The following commands should be used in order to start the graphical user interface(GUI).
```ARMT::run_app()```
## GUI of ARMT
The GUI is developed based on 'shiny' package, and has four page including: Data, Normalization&GSVA, Integration&Analysis, Mutant mapping.<br>
Figure1<br>
The GUI consists of two parts, sidebar panel at left and main panel at right. The sidebar panel is used to input data and adjust parameters; the main panel is used to demonstrate and save the results.

## Data
Figure2<br>
This page of ARMT provides TCGA clinical data and builds .gmt file of arbitrary geneset.<br> 
Through the internal datasets of ARMT, you can get TCGA clinical data of 33 cancer types. <br>
To build a .gmt file, a .csv file should be provided, and it must contain two column as shown in Figure3.<br>
Figure3<br>

## Normalization&GSVA
In this page, ARMT automatically normalize the expression counts matrix to TPM matrix, and use log2(TPM) to process gene set variant analysis(GSVA). The input data must be gene expression matrix with Ensembl ID. <br> 
There are 9 geneset from MSigDb in ARMT internal datasets available for direct selection to process GSVA, and it is also supported to input a .gmt file of arbitrary geneset from 'Data' page. If the input data is from TCGA, you can chose normal, tumor or all samples to process normalize and GSVA(default: tumor).
Figure4<br>

## Integration&Analysis
ARMT can integrate clinical, GSVA, gene expression(TPM), gene mutation data together, and use these data to carry out analysis including: survival analysis&cox proportional hazards regression analysis, differential analysis, enrichment analysis and correlation analysis.<br>

### Integration
The integrated data should contain the common samples. You should input the clinical, GSVA, TPM data by .csv file which can be produced by 'Data' and 'Integration&Analysis' page, and the mutation data should be entered by .maf file. The next analysis is carried out according to these integrated data. You can choose a column in integrated data table as the 'Group Column' to separate the data into multiple group to analyze respectively. Also, you can input your own data in the same format.<br>
Figure5<br>
In this page, clinical data and GSVA data will be merged automatically, and you can enter your interested gene name in your TPM and mutation data to integrate their expression and mutation information by clicking 'Integrate TPM' and 'Integrate Maf' button. The integrated data is showed in main panel.
### Survival analysis&Cox proportional hazards regression analysis
In this page, you should choose two columns in integrated data table as the survival time and status of samples. Each column of integrated data table can be seem as a factor of samples to carried out survival analysis and cox proportional hazards regression analysis.<br> 
Figure6<br>
In survival analysis, if the survival factor  is a continuous variable, the samples will be grouped into two parts(high&low). You can choose multiple factors to obtain multiple analysis results(K-M curve).<br>
Figure7<br>
The cox analysis can focus on single-variable or multiple-variable. If the factor in cox analysis is not numeric variable, one value of this variable should be seem as a reference.<br>

### Differential analysis
ARMT can carry out differential analysis to GSVA score, gene mutation and expression. 
The difference factor should be selected in the integrated data, and the samples will be grouped by it. If the difference factor is a numeric variable, the samples will be grouped into high and low, and the low part will be seem as control group. The group cut off value 't' is a adjustable parameter in ARMT(upper and lower t*100%). If the difference factor is not a numeric variable, you should select two values used to group the samples into experiment group and control group.
You can screen the result by p-value, logFC and FDR(adj.p) after analysis and visualize it in main panel.<br>
Figure8<br>
To analyze difference of GSVA score, the  GSVA result of samples in integrated data should be entered through .csv file, and you can get this file in 'Normalization&GSVA' page.
Figure9<br>
To analyze difference of gene mutation, the mutation information of samples in integrated data should be entered through .maf file. The result is demonstrated with forest plot in main panel.<br>
Figure10<br>
To analyze difference of gene expression(DEG), the counts matrix of samples in integraated data should be entered through .csv file. The result is demonstrated with volcano plot in main panel.<br>

### Enrichment analysis
The enrichment analysis ways in ARMT include GO and KEGG. This enrichment is processed for differential expression gene(DEG) from above differential analysis.
Figure10<br>
You should set up the p-value and q-value cut off before enrichment analysis, and you can select your interested ontology(MF, CC, BP) of GO enrichment. The enriched genes are up regulation and down regulation in DEG, and their enrichment results(up regulation, down regulation and all DEG) are showed in different pages of main panel. These results are demonstrated with bar plot and dot plot.

### Correlation analysis
ARMT can calculate the correlation between any continuous variable factors of integration data, such as TPM and GSVA score. The result is demonstrated in main panel with heat map, and it can be descreened by p-value and correlation coefficient(r, Spearman or Pearson).<br>
Figure11<br>
ARMT can calculate correlation coefficient of each pair-factors.<br>
Figure12<br>
ARMT can calculate correlation coefficient between two list of factors.<br>
If the integrated data is grouped by 'Group Column', ARMT can calculate correlation in each group of data by using single factor or all factors in the list 'Factor1'.<br>
Figure13<br>
For each group of data, the correlation coefficient between single factor and 'Factor2' list is put together in one result, and the single factor name is replaced by group name in correlation matrix.<br>
Figure14<br>
When all factors in 'Factor1' is demonstrated, you can choose one result of multiple groups in sidebar panel to show in main panel.

## Mutant mapping
In corporation with 'maftools', ARMT can visualize the gene mutation in .maf file. There are two plot modes: Maftools Summary and Self-defined. The number of top mutant genes plot out can be set in sidebar panel.<br>
Figure15<br>
In Maftools Summary mode, ARMT can plot the maf summary and illustrate the enrichment of known Oncogenic Signaling Pathways in TCGA cohorts. It is also supported to draw the oncoplot of interested pathway completely.<br>
In Self-defined mode, ARMT can plot the map for any genes or genesets.<br>
Figure16<br>
Figure17<br>

## Format of Input File
The input data of ARMT include geneset, gene expression matrix, gene mutation, GSVA score matrix and clinical information.<br>
#### Geneset
In 'Data' page, the geneset file(.csv, Figure3) is used to create .gmt file, and the .gmt file is the main format of geneset input.<br>
####Gene expression matrix
For gene expression, the counts matrix can be normalized to TPM matrix in ARMT.<br>
Figure18<br>
Counts matrix must be in a .csv file. The row name is gene Ensembl ID  and the column name is the sample ID. It is used to normalization and differential analysis.<br>
Figure19<br>
TPM matrix must be in a .csv file. The row names is gene symbol ID and the column names is the sample ID. It can be integrated with other data.<br>
#### Gene mutation
To get the mutation information, ARMT requires standard .maf file like the mutation annotation format file in TCGA.<br>
##### GSVA score matrix
The GSVA score matrix can be obtained in 'Normalization&GSVA' page.<br>
Figure20<br>
The GSVA score is saved in a matrix of .csv file. The row name is the sample ID, and the column name is geneset name.
#### Clinical information
The TCGA clinical information can be obtained in 'Data' page.<br>
Figure21<br>
The clinical data is input by a .csv file, and the row name is sample ID and column name is sample characteristic. Any information about samples in this format file can be input  as clinical data. <br>
In 'example' folder, we have uploaded some example file to illustrate the input format of above data.






