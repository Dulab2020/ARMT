# ARMT
## Introduction:
ARMT is auto RNA-seq data minning tool
The aim of ARMT is : 
    i) integrate and facilitate the downstream analysis of RNA-seq data,
    ii) provide the way to analyze geneset according to GSVA, 
    iii) explore the correlation of gene expression and mutation.
License: GPL(>=3)
##### The packages from Bioconductoe including:
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
The following commands should be used in order to start the graphical user interface.
```ARMT::run_app()```
## Detail:
The main function of this tool is data analysis. The following Figure1 is the input and adjustment interface of the analyzed data, which can input clinical data, TPM matrix, GSVA data and mutation information data.<br>
![Figure1](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure1.png)
*<p align="center">Figure1. Wstation can integrate clinical, GSVA, TPM and mutation data </p>*
*Red area: File entry area<br>*
*Green area: Select the columns of survival time and survival status<br>*
*Blue area: Enter the genes name of interest and insert relevant TPM and mutation information into the sample information table in yellow area<br>*
*Yellow area: Show the integrated sample data of all input<br>*
All data should be input in .csv format with the exception of mutation data input in .maf format. Input data information must have the common samples information.<br>
Recently, RNA-seq data analysis mainly includes survival analysis, differential analysis, and correlation analysis.<br>
Survival time and survival status should be determined in advance in the data input (Figure1 green area). The following Figure2 show the interface of survival analysis.<br>
![Figure2](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure2.PNG)
*<p align="center">Figure2</p>*
*Red area: File entry area & Selective analysis factor<br>*
*Blue area: Show the table for survival analysis<br>*
The survival analysis results are shown in Figure3.
![Figure3](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure3.PNG)
*<p align="center">Figure3. Survival analysis result of GSVA score</p>*
In addition, cox regression model analysis for single-factor & multiple-factor can also be selected in Figure2 red area, and the results are shown in Figure4.
![Figure4](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure4.PNG)
*<p align="center">Figure4. Cox regression model analysis</p>*
The target of differential analysis can be the counts matrix, GSVA matrix or the mutation information of relevant samples, as shown in Figure5 below.
![Figure5](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure5.png)
*<p align="center">Figure5. Differential analysis grouped by GSVA score</p>*
*Red area: File entry area & Selective type of analysis target<br>*
*Blue area: Setting up differential grouping<br>*
*Green area: Filter for the result(Based on p-value, logFC & FDR)<br>*
*Yellow area: Show of result, Figure5 show the result of Counts analysis. The GSVA differential analysis is similar but no volcano plotting.<br>*
The differential result of the mutation .maf file is shown in a forest plotting as Figure6.
![Figure6](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure6.PNG)
*<p align="center">Figure6. The differential result of the mutation</p>*
The results after differential analysis can be filtered according to p-value, FDR and logFC, and the filtered differential genes can be further enrichment analysis as shown in Figure7.
![Figure7](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure7.png)
*<p align="center">Figure7. Enrichment analysis</p>*
*Red area: Filter for the enrichment analysis(pcut & qcut)<br>*
*Blue area: Selective ONTOLOGY for GO enrichment(ALL、MF、CC & BP)<br>*
*Green area: The number of pathways show in the pictures.<br>*
*Yellow area: The result of enrichment analysis. The result of GO is similar to that of KEGG.<br>*
The correlation analysis is directed toward the input data of Figure1. One set of features can be selected for pairwise correlation analysis or two sets of those can be selected for correlation analysis between each other as shown in Figure8.
![Figure8](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure8.png)
*<p align="center">Figure8. Correlation analysis between TPM and GSVA score</p>*
*Red area: Selective features for correlation analysis.<br>*
*Blue area: Filter of result(based on r & p-value).<br>*
*Green area: Correlation matrix.<br>*
*Yellow area: Table to show the detail of pairwise features correlation.<br>*
Besides, there is also the heatmap about correlation as shown in Figure9.
![Figure9](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure9.PNG)
*<p align="center">Figure9. Heatmap from correlation analysis</p>*
For the convenience of use, Wstation also provides TCGA data download function, and has the internal TCGA clinical data for use. It is recommended to directly obtain the internal data, as it is better suited for this tool. This function is on the first page of Wstation as shown in Figure10.
![Figure10](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure10.png)
*<p align="center">Figure10</p>*
*Red area: Selective cancer type & data type(Counts or Clinical).<br>*
*Blue area: File entry area for self-identified .gmt file.<br>*
*Green area: The internal clinical data.<br>*
This page also provides the function to generate the .gmt file, with two columns list table containing the gene sets name and genes name in a .csv file, the corresponding .gmt file can be generated for the following GSVA calculation.
The calculation of GSVA requires the input of the Counts matrix (.csv) and the gene set file (.gmt). The gene set of MSigDB has been built in Wstation and can be selected here. Figure11 is the page to calculate GSVA score.
![Figure11](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure11.png)
*<p align="center">Figure11. The interface of GSVA calculation</p>*
*Red area: Counts file entry area & Sample type selection<br>*
*Blue area: Selective MSigDB gene sets & Gene set file entry area<br>*
*Green area: GSVA matrix & TPM matrix.<br>*
For exploring the relationship between RNA-seq data and gene mutations, Wstation integrates some functions of maftools. It requires the .maf file. The general state of mutation can be plot out in summary as shown in Figure12.
![Figure12](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure12.png)
*<p align="center">Figure12. The summary result of maftools</p>*
*Red area: File entry area.<br>*
*Blue area: Plot mode selection including ‘Summary’ & ‘Self-identified’.<br>*
*Green area: The maximum number of genes to plot out.<br>*
*Yellow area: Graphic of mutation information.<br>*
In summary plot mode, the pathway display section can choose to display the mutation of the selected pathway, as shown in Figure13.
![Figure13](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure13.png)
*<p align="center">Figure13</p>*
*Red area: Mutations in all known TCGA cohort pathways.<br>*
*Blue area: All genes mutation in the selected pathway oncoplot.<br>*
The main function of this part is the self-identified plot mode, which allows to select the interesting genes or pathway to plot out. The pathways should be input in a .gmt file as gene sets. Figure14 demonstrates the self-identified plot mode.
![Figure14](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure14.png)
*<p align="center">Figure14. The self-identified plot of some genes' mutation</p>*
*Black area: Selection of genes or pathways. Figure14 show the genes selection.<br>*
*Red area: Selective types of mutation & Entry of genes name.<br>*
*Blue area: The distribution of genes mutation is showed in the boxplot.<br>*
*Green area: Co-occurrence or mutually exclusive between pairwise genes.<br>*
*Yellow area: Genes data and sample data entry area. These data can be showed with mutation information(Figure15).<br>*
With the input data of samples’ GSVA score and the differential FDR of genes, the RNA-seq data analysis result will be demonstrated with mutation as shown in Figure15.
![Figure15](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure15.png)
*<p align="center">Figure15</p>*
The sample data and gene data are required to have two columns in .csv file. The first column is the samples name or genes name, and the second column is the continuous variable data. These data should have the common samples and genes with .maf file.<br>
Also, if the you are interested in pathways, the pathway can be selected, as shown in Figure16.
![Figure16](https://github.com/Dulab2020/Wstation/blob/master/Figure/Figure16.png)
*<p align="center">Figure16</p>*
*Red area: The .gmt file entry & The selection of pathways in the .gmt file. <br>*
*Blue area: Mutation state of all pathways.<br>*
*Green area: An oncoplot similar to Figure15 to demonstrate all genes mutation in the selected pathways.<br>*













