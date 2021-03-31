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
  data("interdata", package = 'Wstation')
  geneLength <- interdata$gLen
  geneSet <- interdata$gSet
  fluidPage( theme = shinytheme('flatly'),
             tabsetPanel(useShinyjs(),
                         #获取文件的UI
                         tabPanel('New', titlePanel('Creater'),
                                  sidebarLayout(
                                    sidebarPanel(
                                      h4('TCGADownload'),
                                      checkboxGroupInput('cancerType', 'Choose the cancer type:', inline = T ,choices = allCancer),
                                      checkboxInput('all','All'),
                                      downloadButton('getInterClinical', 'Get internal clinical data'),
                                      br(),
                                      br(),
                                      h4('Geneset Creater'),
                                      fileInput('genesetCsv', 'Choose geneset .csv file:',
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
                         tabPanel('Counts',titlePanel("GSVA"),
                                  sidebarLayout(
                                    sidebarPanel(
                                      fileInput('countsMatrix', 'Choose counts file:',
                                                accept=c('text/csv', 'text/comma-separated-values,text/plain')),
                                      selectInput('sampleType', 'Sample Type', choices = c('All' = 'all', 'Tumor' = 'tumor', 'Normal'= 'normal'), selected = 'tumor'),
                                      checkboxGroupInput('genesetlist', 'Choose the gene set:',
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
                                      fileInput('gmtFile', 'Choose *.gmt file:', accept = '.gmt'),
                                      actionButton('GSVA', 'Calculate')
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
                         tabPanel('Clinical',titlePanel("Clinical"),
                                  fluidRow(
                                    column(width = 3,
                                           fileInput('clinicalData', 'Choose clinical .csv file:',
                                                     accept='.csv'),
                                           fileInput('gsvaData', 'Choose GSVA .csv file:',
                                                     accept='.csv'),
                                           fileInput('tpmData', 'Choose TPM .csv file:',
                                                     accept='.csv'),
                                           fileInput('mafData', 'Choose maf.gz file:',
                                                     accept=c('.maf.gz', '.maf'))
                                    ),
                                    column(width = 3,
                                           uiOutput('chooseTime'),
                                           uiOutput('chooseStatus'),
                                           uiOutput('chooseGroup'),
                                           fluidRow(
                                             column(width = 9,
                                                    textAreaInput('keyGeneList', 'Gene names', height = 110, placeholder = 'TP53,EEF2,RPS2......')
                                             ),
                                             column(width = 1,
                                                    br(),
                                                    actionButton('addTpm','TPM'),
                                                    br(),
                                                    br(),
                                                    actionButton('addMaf','MAF')
                                             )
                                           )
                                    ),
                                    column(width = 6,column(width = 12,
                                                            DT::DTOutput('clinDataShow'), style = "overflow-x: scroll;"),
                                           downloadButton('getNewClinical', 'Save')
                                    )
                                  ),
                                  #数据分析的UI
                                  fluidRow(hr()),
                                  fluidRow(tabsetPanel(
                                    #生存分析
                                    tabPanel('Survival',
                                             sidebarLayout(
                                               sidebarPanel(
                                                 selectInput('surWay', 'Choose your calculate way:', choices = 
                                                               c('Survival' = 'sur', 'Single Cox' = 'singleCox', 'Multiple Cox' = 'multipleCox')),
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
                                    tabPanel('DEA', sidebarLayout(
                                      sidebarPanel(
                                        h3('Differential Analysis'),
                                        fileInput('deaData', 'Input the file to analysis:', accept = c('.csv', '.maf.gz', '.maf')),
                                        radioButtons('deaWay', 'Choose the tool to DEA: ', choices = c('Counts' = 'edg','GSVA' = 'lm', 'Maf' = 'maf'), inline = TRUE),
                                        checkboxInput('deaTCGAFlag', 'TCGA', value = TRUE),
                                        uiOutput('chooseDeaFactor'),
                                        sliderInput('groupCutOff', 'The group cut off:', min = 0, max = 0.5, value = 0.3, step = 0.01),
                                        fluidRow(column(width = 6, uiOutput('chooseExperience')),
                                                 column(width = 6, uiOutput('chooseControl'))
                                        ),
                                        actionButton('calDiffer', 'Calculate'),
                                        uiOutput('deaGroupUI'),
                                        sliderInput('logFCCutOff', 'logFC cut off:', min = 0, max = 2, value = 1, step = 0.01),
                                        sliderInput('pCutOff', 'P-value cut off:', min = 0, max = 1, value = 0.05, step = 0.01),
                                        sliderInput('fdrCutOff', 'FDR(or adj.P) cut off:', min = 0, max = 1, value = 0.02, step = 0.01),
                                        hr(),
                                        h3('Enrich Analysis'), 
                                        fluidRow(sliderInput('enrichPCut', 'P Cut Off:', min = 0, max = 1, step = 0.01, value = 0.1),
                                                 sliderInput('enrichQCut', 'q Cut Off:', min = 0, max = 1, step = 0.01, value = 0.2),
                                                 selectInput('enrichOnto', 'GO:ONTOLOGY', choices = c('ALL', 'MF', 'CC', 'BP')),
                                                 column(width = 6, actionButton('calGO', 'GO Analysis')),
                                                 column(width = 6, actionButton('calKegg', 'KEGG Analysis'))),
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
                                    tabPanel('COR', sidebarLayout(
                                      sidebarPanel(
                                        h2('Correlation'),
                                        uiOutput('chooseCorFactor1'),
                                        uiOutput('chooseCorFactor2'),
                                        selectInput('corWay', 'Correlation Coefficient:', choices = c('Spearman'= 'spearman', 'Pearson' = 'pearson')),
                                        actionButton('calCor', 'Calculate'),
                                        uiOutput('corActiveUI'),
                                        sliderInput('corrCut', 'r Cut Off:', min = 0, max = 1, value = 0),
                                        sliderInput('corpCut', 'p Cut Off:', min = 0, max = 1, value = 1)
                                      ),
                                      mainPanel(
                                        uiOutput('corShow'),
                                        sliderInput('heatSize', 'Heatmap Size:', min = 10, max = 100, value = 20),
                                        uiOutput('corHeatShow'))
                                    ))
                                  ))
                         ),
                         #maf文件可视化
                         tabPanel('Maftools', titlePanel('Mutation Visualization'),
                                  sidebarLayout(sidebarPanel(
                                    fileInput('mafVisual', 'Choose maf.gz file:',
                                              accept=c('.maf.gz','.maf')),
                                    checkboxInput('mafTCGAFlag', 'TCGA', value = TRUE),
                                    selectInput('mafShowMode', 'Plot Mode:', choices = c('Self-defined' = 'self',
                                                                                         'Summary' = 'sum')),
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
      app_title = 'Wstation'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

