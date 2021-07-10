#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinythemes
#' @import DT
#' @import shinyjs
#' @noRd
app_ui <- function(request) {
  data("interdata", package = 'ARMT')
  geneLength <- interdata$gLen
  geneSet <- interdata$gSet
  fluidPage( theme = shinytheme('flatly'),
             tabsetPanel(useShinyjs(),
                         selected = 'Data',
                         #获取文件的UI
                         tabPanel('Data', titlePanel('Data'),
                                  sidebarLayout(
                                    sidebarPanel(
                                      h4('Clinical Data'),
                                      checkboxGroupInput('cancerType', 'Choose the cancer type:', inline = T ,choices = allCancer),
                                      checkboxInput('all','All'),
                                      downloadButton('getInterClinical', 'Get internal TCGA clinical data'),
                                      br(),
                                      br(),
                                      h4('Geneset Data'),
                                      fileInput('genesetCsv', 'Input geneset .csv file:',
                                                accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                      actionButton('creatGmt', 'Creat gmt file')
                                    ),
                                    mainPanel(
                                      uiOutput('titleGeneset'),
                                      DT::dataTableOutput('genesetMatrixShow'),
                                      br(),
                                      uiOutput('titleInterClinical'),
                                      DT::dataTableOutput('interClinicalShow')
                                    )
                                  )       
                         ),
                         #GSVA的UI
                         tabPanel('Normalization & GSVA',titlePanel("Normalization & GSVA"),
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput('countsMatrix', 'Input counts file with Ensembl ID:',
                                                accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                      selectInput('sampleType', 'Sample Type', choices = c('All' = 'all', 'Tumor' = 'tumor', 'Normal'= 'normal'), selected = 'tumor'),
                                      checkboxGroupInput('genesetlist', 'Choose the gene set in MSigDB:',
                                                         choices = c(
                                                           'H: hallmark gene sets' = 'h',
                                                           'C1: positional gene sets' = 'c1',
                                                           'C2: curated gene sets' = 'c2',
                                                           'C3: regulatory target gene sets' = 'c3',
                                                           'C4: computational gene sets' = 'c4',
                                                           'C5: ontology gene sets' = 'c5',
                                                           'C6: oncogenic signature gene sets' = 'c6',
                                                           'C7: immunologic signature gene sets' = 'c7',
                                                           'C8: cell type signature gene sets' = 'c8'
                                                         )
                                      ),
                                      fileInput('gmtFile', 'Input gene set file(.gmt):', accept = '.gmt'),
                                      actionButton('GSVA', 'Calculate GSVA score')
                                    ),
                                    mainPanel(
                                      uiOutput('titleGSVA'),
                                      DT::dataTableOutput('gsvaShow'),
                                      uiOutput('gsvaDownload'),
                                      br(),
                                      uiOutput('titleTpm'),
                                      DT::dataTableOutput('tpmShow'),
                                      uiOutput('tpmDownload')
                                    )
                                  )
                         ), 
                         #临床信息整理的UI
                         tabPanel('Integration & Analysis',titlePanel("Data Integration"),
                                  fluidRow(
                                    column(width = 3,
                                           fileInput('clinicalData', 'Input clinical file(.csv):',
                                                     accept='.csv'),
                                           fileInput('gsvaData', 'Input GSVA file(.csv):',
                                                     accept='.csv'),
                                           fileInput('tpmData', 'Input TPM file(.csv):',
                                                     accept='.csv'),
                                           fileInput('mafData', 'Input mutation maf file:',
                                                     accept=c('.maf.gz', '.maf'))
                                    ),
                                    column(width = 3,
                                           uiOutput('chooseGroup'),
                                           fluidRow(
                                             textAreaInput('keyGeneList', 'Input Gene names', width = 330, height = 220, placeholder = 'TP53,EEF2,RPS2......'),
                                             column(width = 6,
                                                    actionButton('addTpm','Integrate TPM')),
                                             column(width = 6,
                                                    actionButton('addMaf','Integrate MAF'))
                                           )
                                    ),
                                    column(width = 6,column(width = 12,
                                                            DT::DTOutput('clinDataShow'), style = "overflow-x: scroll;"),
                                           downloadButton('getNewClinical', 'Save')
                                    )
                                  ),
                                  #数据分析的UI
                                  fluidRow(hr()),
                                  fluidRow(titlePanel('Analysis')),
                                  fluidRow(tabsetPanel(
                                    #生存分析
                                    tabPanel('Survival',
                                             sidebarLayout(
                                               sidebarPanel(
                                                 uiOutput('chooseTime'),
                                                 uiOutput('chooseStatus'),
                                                 selectInput('surWay', 'Choose analysis way:', choices = 
                                                               c('Kaplan-Meier Curve' = 'sur', 'Univarivity Cox Analysis' = 'singleCox', 'Multivarivity Cox Analysis' = 'multipleCox')),
                                                 uiOutput('surGroupUI'),
                                                 uiOutput('chooseSurFactor'),
                                                 uiOutput('coxRefUI'),
                                                 uiOutput('chooseEvent'),
                                                 actionButton('calSur', 'Calculate')
                                               ), 
                                               mainPanel(
                                                 sliderInput('surPlotSize', 'Size:', min = 100, max = 800, value = 300),
                                                 sliderInput('forestSize', 'Size:', min = 600, max = 1600, value = 800),
                                                 uiOutput('surShow'),
                                                 hr(),
                                                 DT::DTOutput('showtest')
                                               )
                                             )
                                    ),
                                    #差异分析
                                    tabPanel('Differential Analysis', sidebarLayout(
                                      sidebarPanel(
                                        h3('Differential Analysis'),
                                        fileInput('deaData', 'Input the file to differential analysis:', accept = c('.csv', '.maf.gz', '.maf')),
                                        radioButtons('deaWay', 'Choose the type of differential object: ', choices = c('Gene expression' = 'edg','GSVA score' = 'lm', 'Mutation' = 'maf'), inline = TRUE),
                                        selectInput('deaNorm', 'Normalize way for counts:', choices = c('TMM', 'TMMwsp', 'RLE', 'upperquartile', 'none')),
                                        checkboxInput('deaTCGAFlag', 'TCGA', value = TRUE),
                                        uiOutput('chooseDeaFactor'),
                                        sliderInput('groupCutOff', 'The group cut off value:', min = 0, max = 0.5, value = 0.3, step = 0.01),
                                        fluidRow(column(width = 6, uiOutput('chooseExperience')),
                                                 column(width = 6, uiOutput('chooseControl'))
                                        ),
                                        actionButton('calDiffer', 'Calculate'),
                                        uiOutput('deaGroupUI'),
                                        sliderInput('logFCCutOff', 'logFC cut off(>):', min = 0, max = 10, value = 2, step = 0.5),
                                        sliderInput('pCutOff', 'P-value cut off(<):', min = 0, max = 0.2, value = 0.05, step = 0.01),
                                        sliderInput('fdrCutOff', 'FDR(or adj.P) cut off(<):', min = 0, max = 0.2, value = 0.05, step = 0.01),
                                        hr(),
                                        h3('Enrichment Analysis'), 
                                        fluidRow(radioButtons('enrichGene', 'Genes for analysis:', choices = c('DEG' = 'DEG','Inputed gene list' = 'geneInput'), inline = TRUE),
                                                 textAreaInput('enrichGeneInput', 'Input Gene names for analysis', width = 415, height = 100, placeholder = 'TP53,EEF2,RPS2......'),
                                                 sliderInput('enrichPCut', 'p Cut Off(<):', min = 0, max = 0.2, value = 0.05, step = 0.01),
                                                 sliderInput('enrichQCut', 'q Cut Off(<):', min = 0, max = 0.2, value = 0.05, step = 0.01),
                                                 selectInput('enrichOnto', 'GO:ONTOLOGY', choices = c('ALL', 'MF', 'CC', 'BP')),
                                                 column(width = 6, actionButton('calGO', 'GO Enrichment')),
                                                 column(width = 6, actionButton('calKegg', 'KEGG Enrichment'))),
                                      ),
                                      mainPanel(
                                        uiOutput('showDea'),
                                        uiOutput('volcanoPictureUI'),
                                        downloadButton('volcanoSave', 'Save(.pdf)'), #火山图下载代码放在deaScreen中
                                        hr(),
                                        uiOutput('goHead'), #GO结果设置以及标题
                                        uiOutput('goShow'),
                                        uiOutput('keggHead'), #KEGG结果设置以及标题
                                        uiOutput('keggShow')
                                      )
                                    )),
                                    #相关性分析
                                    tabPanel('Correlation', sidebarLayout(
                                      sidebarPanel(
                                        h2('Correlation Analysis'),
                                        uiOutput('chooseCorFactor1'),
                                        uiOutput('chooseCorFactor2'),
                                        selectInput('corWay', 'Correlation Coefficient:', choices = c('Spearman'= 'spearman', 'Pearson' = 'pearson')),
                                        actionButton('calCor', 'Calculate'),
                                        uiOutput('corActiveUI'),
                                        sliderInput('corrCut', 'r Cut Off(>):', min = 0, max = 1, value = 0, step = 0.01),
                                        sliderInput('corpCut', 'p Cut Off(<):', min = 0, max = 2, value = 0.05, step = 0.01)
                                      ),
                                      mainPanel(
                                        uiOutput('corShow'),
                                        sliderInput('heatSize', 'Heatmap Size:', min = 10, max = 100, value = 20),
                                        uiOutput('corHeatShow'))
                                    ))
                                  ))
                         ),
                         #maf文件可视化
                         tabPanel('Mutant mapping', titlePanel('Mutation Visualization'),
                                  sidebarLayout(sidebarPanel(
                                    fileInput('mafVisual', 'Input mutation profile(maf):',
                                              accept=c('.maf.gz','.maf')),
                                    checkboxInput('mafTCGAFlag', 'TCGA', value = TRUE),
                                    selectInput('mafShowMode', 'Plot Mode:', choices = c('Self-defined' = 'self',
                                                                                         'Maftools Summary' = 'sum')),
                                    numericInput('mafTop', 'Top:', min = 1, max = 50, step = 1, value = 20),
                                    uiOutput('mafIn'),
                                    actionButton('mafPlot', 'Plot Out')
                                  ), 
                                  mainPanel(
                                    uiOutput('mafShow')
                                  ))
                         )
             )
             
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'ARMT'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

