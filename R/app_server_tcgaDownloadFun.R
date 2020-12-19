#' @import TCGAbiolinks
#' @import dplyr
#' @import DT
#' @import SummarizedExperiment
#' @import shiny

downloadClinical <- function(cancer_type){
  #下载临床数据
clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
write.csv(clinical,file = paste0('TCGAdownload/',paste(cancer_type,"clinical.csv",sep = "-")))
}

downloadCounts <- function(cancer_type){
  #下载rna-seq的counts数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "Gene Expression Quantification", 
                    workflow.type = "HTSeq - Counts")
  GDCdownload(query, method = "api", files.per.chunk = 100)
  expdat <- GDCprepare(query = query)
  count_matrix<-assay(expdat)
  write.csv(count_matrix,file = paste0('TCGAdownload/',paste(cancer_type,"Counts.csv",sep = "-")))
}

downloadmiRNA <- function(cancer_type){
  #下载miRNA数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Transcriptome Profiling", 
                    data.type = "miRNA Expression Quantification", 
                    workflow.type = "BCGSC miRNA Profiling")
  
  GDCdownload(query, method = "api", files.per.chunk = 50)
  expdat <- GDCprepare(query = query)
  count_matrix<-assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"miRNA.csv",sep = "-"))
}

downloadCNV <- function(cancer_type){
  #下载Copy Number Variation数据
  query <- GDCquery(project = cancer_type, 
                    data.category = "Copy Number Variation", 
                    data.type = "Copy Number Segment")
  
  GDCdownload(query, method = "api", files.per.chunk = 50)
  expdat <- GDCprepare(query = query)
  count_matrix<-assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"Copy-Number-Variation.csv",sep = "-"))
}

downloadMethylation <- function(cancer_type){
  #下载甲基化数据
  query.met <- GDCquery(project =cancer_type,
                        legacy = TRUE,
                        data.category = "DNA methylation")
  GDCdownload(query.met, method = "api", files.per.chunk = 300)
  expdat <- GDCprepare(query = query)
  count_matrix<-assay(expdat)
  write.csv(count_matrix,file = paste(cancer_type,"methylation.csv",sep = "-"))
}

#下面填入要下载的癌症种类和要下载的数据
downloadTCGAdata <- function(cancerList, dataList){
  dir.create('TCGAdownload', showWarnings = FALSE)
  progress <- Progress$new(min = 1, max = 5)
  for (i in cancerList) {
    t <- 1
    progress$set(message = paste0(paste0('Download ',i),' data:'), detail = 'This may take a while...')
    cancer <- paste("TCGA",i,sep="-")
    if ('clinical' %in% dataList){downloadClinical(cancer)}
    progress$set(value = t)
    t <- t+1
    if ('counts' %in% dataList){downloadCounts(cancer)}
    progress$set(value = t)
    t <- t+1
    if ('miRNA' %in% dataList){downloadmiRNA(cancer)}
    progress$set(value = t)
    t <- t+1
    if ('CNV' %in% dataList){downloadCNV(cancer)}
    progress$set(value = t)
    t <- t+1
    if ('methylation' %in% dataList){downloadMethylation(cancer)}
    progress$set(value = t)
    t <- t+1
  }
}
