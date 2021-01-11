#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinythemes
#' @import DT
#' @import shinyjs
#' @import data.table
#' @import maftools
#' @import GSVA
#' @import GSVAdata
#' @import maftools
#' @import survival
#' @import survminer
#' @import edgeR
#' @import limma
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import DOSE
#' @import ggplot2
#' @import ggpubr
#' @import psych
#' @import pheatmap
#' @import GSEABase
#' @import ggcorrplot
#' @import patchwork
#' @noRd
app_server <- function( input, output, session ) {
  options(shiny.maxRequestSize=-1) # Remove limit of upload
  options(shiny.deprecation.messages=FALSE)
  options(warn =-1)
  #下载TCGA文件
  observe({
    allFlag <- input$all
    if (allFlag){updateCheckboxGroupInput(session, 'cancerType','Choose the cancer type:',inline = T, 
                                          choices = allCancer, selected =  allCancer)}
    else{updateCheckboxGroupInput(session, 'cancerType','Choose the cancer type:',inline = T, 
                                  choices = allCancer)}
  })
  observeEvent(input$downloadTCGA, {downloadTCGAdata(input$cancerType, input$dataList)})
  #获取内部临床信息文件
  chooseInterClinical <- reactive({
    if(!is.null(input$cancerType)){
      output$titleInterClinical <- renderUI({h3('ClinicalData')})
      shinyjs::show(id = 'titleInterClinical')}
    else{
      shinyjs::hide(id = 'titleInterClinical')}
    choInterClinical <- interClinical[which(interClinical$type == input$cancerType),]
    naFlag <- apply(choInterClinical, 2, function(x) all(x =='#N/A'))
    choInterClinical[,which(!naFlag)]
  })
  output$interClinicalShow<-DT::renderDataTable({chooseInterClinical()}, edit = TRUE)
  #创建gmt文件
  genesetMatrix <- reactive(readMatrix(input$genesetCsv, FALSE))
  output$titleGeneset <- renderUI({
    if(is.null(input$genesetCsv))
      NULL
    else
      h3('GeneSet')
  })
  observeEvent(input$creatGmt,{
    writeGmt(genesetMatrix(), paste(strsplit(input$genesetCsv$name, '[.]')[[1]][1], '.gmt', sep='') )
  })
  output$genesetMatrixShow <- DT::renderDataTable(genesetMatrix())
  #tpm,tpmMatrix改为tpm列表，包含有counts矩阵和分类类型
  tpmMatrix <- reactive({
    countsData <- readMatrix(input$countsMatrix)
    if (!is.null(countsData)){
      output$tpmDownload <- renderUI({
        fluidRow(column(width = 4, downloadButton('downloadTPM', 'Save TPM')),
                 column(width = 4, downloadButton('downloadCounts', 'Save Counts')))
      })
      output$titleTpm <- renderUI({h2('TPM')})
      if(input$sampleType == 'all'){
        list(tpm = countsToTPM(countsData), counts = countsData, type = 'all')}
      else if(input$sampleType == 'tumor'){
        cts <- classifyCounts(countsData)$tm
        list(tpm = countsToTPM(cts), counts = cts, type = 'tumor')
      }
      else if(input$sampleType ==  'normal'){
        cts <- classifyCounts(countsData)$nr
        list(tpm = countsToTPM(cts), counts = cts, type = 'norm')
      }
    }
  })
  output$tpmShow <- DT::renderDataTable({tpmMatrix()$tpm})
  #GSVA
  gsvaMatrix <- reactive({
    if(any(!is.null(input$genesetlist),!is.null(input$gmtFile))){
      gsvaCal(tpmMatrix()$tpm, input$genesetlist, gmtFlag = !is.null(input$gmtFile), gmtPath = input$gmtFile$datapath)
    }    
  })
  observeEvent(input$GSVA,{
    if(any(!is.null(input$genesetlist), !is.null(input$gmtFile))){
      gsvamat <- gsvaMatrix()
      output$gsvaShow <- DT::renderDataTable({gsvamat})
      output$titleGSVA <- renderUI({h2('GSVA')})
      output$gsvaDownload <- renderUI({downloadButton('downloadGSVA', 'Download GSVA')})
      shinyjs::show(id = 'gsvaShow')
      shinyjs::show(id = 'titleGSVA')
      shinyjs::show(id = 'gsvaDownload')
    }
    else{
      shinyjs::hide(id = 'gsvaShow')
      shinyjs::hide(id = 'titleGSVA')
      shinyjs::hide(id = 'gsvaDownload')
    }
  })
  #临床信息整理动态UI
  observe({
    if(all(is.null(input$tpmData),is.null(input$mafData))){
      shinyjs::hide(id = 'keyGeneList')
      shinyjs::hide(id = 'addTpm')
      shinyjs::hide(id = 'addMaf')
    }
    else{
      shinyjs::show(id = 'keyGeneList')
      shinyjs::show(id = 'addTpm')
      shinyjs::show(id = 'addMaf')
    }
  })
  observe({
    #临床信息整理区域
    if(!is.null(input$clinicalData)){
      clinData <<- readMatrix(input$clinicalData)
      output$chooseTime <- renderUI(selectInput('surTime', 'Survival Time:', 
                                                colnames(clinData), selected = 'OS.time'))
      output$chooseStatus <- renderUI(selectInput('surSta', 'Survival Status:', 
                                                  colnames(clinData), selected = 'OS'))
      output$chooseGroup <- renderUI(selectInput('groupInfo', 'Group Column:', 
                                                 c('No group',colnames(clinData)), selected = 'No group'))
      shinyjs::show(id = 'getNewClinical')
    }
    else{shinyjs::hide(id = 'getNewClinical')
      clinData <<- NULL}
    newClin <<- clinData
  })
  gsvaClinalData <- reactive({
    input$clinicalData
    if(!is.null(input$gsvaData)){clinData <- data.frame(merge(readMatrix(input$gsvaData)
                                                              , clinData, by = 0, all = FALSE), row.names = 1, stringsAsFactors = FALSE)}
    clinData   
  })
  observe({newClin <<- gsvaClinalData()})
  #临床信息整理
  newClinalData <- reactive({
    input$clinicalData
    input$gsvaData
    input$addTpm
    input$addMaf
    tempClin <- newClin
    if(input$addTpm == 0){addTpmFlag <<- input$addTpm}
    if(input$addMaf == 0){addMafFlag <<- input$addMaf}
    if(input$addTpm != addTpmFlag){
      tempClin <- data.frame(merge(t(readMatrix(input$tpmData))[,unlist(strsplit(input$keyGeneList, split = ',')), drop = FALSE]
                                   , tempClin, by = 0, all = FALSE), row.names = 1, stringsAsFactors = FALSE)
      addTpmFlag <<- input$addTpm
    }
    if(input$addMaf != addMafFlag){
      tempClin <- data.frame(merge(getMaf(unlist(strsplit(isolate(input$keyGeneList), split = ',')), input$mafData$datapath)
                                   , tempClin, by = 0, all = FALSE), row.names = 1, stringsAsFactors = FALSE)
      addMafFlag <<- input$addMaf
    }
    newClin <<- tempClin
    tempClin  
  })
  output$clinDataShow <- DT::renderDT(newClinalData(), options=list(pageLength = 2), selection = list(target = 'column'))
  #分析因子选择框
  output$chooseSurFactor <- renderUI(selectInput('surFactor', label = 'Factor:', multiple = TRUE, choices = colnames(newClinalData())))
  output$chooseDeaFactor <- renderUI(selectInput('deaFactor', label = 'Group by: ', choices = colnames(newClinalData())))
  output$chooseCorFactor1 <- renderUI(selectInput('corFactor1', 'Factor1:', choices = colnames(newClinalData()), multiple = TRUE))
  output$chooseCorFactor2 <- renderUI(selectInput('corFactor2', 'Factor2:', choices = colnames(newClinalData()), multiple = TRUE))
  #生存分析计算矩阵
  surClin <- reactive({
    if(!is.null(newClinalData())){
      if(!is.numeric(newClinalData()[,input$surSta])){
        factorSta <- tolower(levels(as.factor(newClinalData()[,input$surSta])))
        if(any(all(factorSta == c('alive','dead')), all(factorSta == c('dead','alive')))){
          shinyjs::hide('deadEvent')
          toSurStatus(newClinalData(), input$surSta, 'dead')
        }
        else{
          output$chooseEvent <- renderUI(selectInput('deadEvent', 'Input the dead event: ', choices = factorSta))
          shinyjs::show('deadEvent')
          toSurStatus(newClinalData(), input$surSta, input$deadEvent)
        }
      }
      else{shinyjs::hide('deadEvent')
        newClinalData()}
    }
  })
  output$showtest <- DT::renderDT({
    tryCatch({surClin()[unique(c(input$surTime, input$surSta, input$surFactor))]}, error = function(x){NULL})
  },  options=list(pageLength = 1))
  shinyjs::hide('surPlotSize')
  #生存分析
  surResult <- eventReactive(input$calSur,{
    if(input$surWay == 'sur'){
      allPlot <- surAnalysis(surClin(), input$surTime, input$surSta, input$surFactor[[1]])$plot
      for(fac in input$surFactor){
        if(fac == input$surFactor[[1]]){next}
        allPlot <- allPlot + surAnalysis(surClin(), input$surTime, input$surSta, fac)$plot
        
      }
      allPlot
    }
    else if(input$surWay == 'singleCox'){
      singleCox(surClin(), input$surTime, input$surSta, input$surFactor[[1]])
    }
    else if(input$surWay == 'multipleCox'){
      multipleCox(surClin(), input$surTime, input$surSta, input$surFactor)
    }
  })
  output$surShow <- renderUI({
    if(!is.null(surResult())){
      if(isolate(input$surWay) == 'sur'){
        output$surPlot <- renderPlot(surResult(),
                                     width = ceiling(sqrt(length(isolate(input$surFactor))))*input$surPlotSize,
                                     height = round(sqrt(length(isolate(input$surFactor))))*input$surPlotSize)
        shinyjs::show('surPlotSize')
        plotOutput('surPlot')
      }
      else if(isolate(input$surWay) == 'singleCox'){
        shinyjs::hide('surPlotSize')
        output$singleCoxShow <- renderPrint(surResult())
        verbatimTextOutput('singleCoxShow')
      }
      else if(isolate(input$surWay) == 'multipleCox'){
        shinyjs::hide('surPlotSize')
        output$multipleCoxShow <- renderPrint(summary(surResult()$re))
        output$multipleCoxForest <- renderPlot(ggforest(surResult()$or, data = sbToLine(surClin()), main = 'Hazard ratio', 
                                                        cpositions = c(0.01,input$forestColSite,0.4), 
                                                        fontsize = input$forestFont, refLabel = '1', noDigits = 4), width = 1000)
        fluidRow(verbatimTextOutput('multipleCoxShow'),
                 column(width = 4, sliderInput('forestColSite', 'Column Site:', min = 0.02, max = 0.35, step = 0.01, value = 0.2)),
                 column(width = 4, sliderInput('forestFont', 'Font Size:', min = 0.5, max = 1.5, step = 0.1, value = 1)),
                 plotOutput('multipleCoxForest'))
      }
    }
    else{NULL}
  })
  #差异分析动态UI
  observe({
    temp <- tryCatch(newClinalData()[,input$deaFactor], error = function(x){NULL})
    if(!is.numeric(temp)){
      factorDea <- levels(as.factor(temp))
      output$chooseExperience <- renderUI(selectInput('exGroup','Experience group', choices = factorDea))
      output$chooseControl <- renderUI(selectInput('ctGroup','Control group', choices = factorDea))
      shinyjs::show('exGroup')
      shinyjs::show('ctGroup')
    }
    else{
      shinyjs::hide('exGroup')
      shinyjs::hide('ctGroup')
    }
  })
  deaObj <- reactive({if(input$deaWay == 'maf'){read.maf(input$deaData$datapath, isTCGA = input$deaTCGAFlag)}
    else{readMatrix(input$deaData)}
  })
  #差异分析
  deaResult <- eventReactive(input$calDiffer,{
    if(!is.null(deaObj())){
      groupCondit <- newClinalData()[input$deaFactor]
      if(is.numeric(groupCondit[,1])){
        groupCondit[,1] <- seriesToDiscrete(groupCondit[,1], input$groupCutOff)
        ex <- 'High'
        ct <- 'Low'
      }
      else{
        ex <- input$exGroup
        ct <- input$ctGroup
      }
      if(input$deaWay == 'edg'){
        transferID(deaEdgeR(deaObj(), groupCondit, ex, ct))
      }
      else if(input$deaWay == 'lm'){
        deaLimma(deaObj(), groupCondit, ex, ct)
      }
      else if(input$deaWay == 'maf'){
        mafDiffer(deaObj(), groupCondit, ex, ct)
      }
      else{print('?')}
    }
  })
  deaScreen <- reactive({
    if(!is.null(deaResult())){
      if(input$deaWay == 'edg'){
        differResult <- deaResult()[deaResult()$PValue <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$FDR <= input$fdrCutOff, , drop = FALSE]
      }
      else if(input$deaWay == 'lm'){
        differResult <- deaResult()[deaResult()$P.Value <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$adj.P.Val <= input$fdrCutOff, , drop = FALSE]
      }
      else if(input$deaWay == 'maf'){
        differResult <- deaResult()$results
        differResult <- differResult[differResult$pval <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$adjPval <= input$fdrCutOff, , drop = FALSE]
      }
      #logFC过滤没maf差异的份
      if(input$deaWay != 'maf'){
        if(input$logFCCutOff == 0)
        {differResult$Status <- cut(differResult$logFC, c(-Inf, input$logFCCutOff, Inf), c('Down','Up'))}
        else
        {differResult$Status <- cut(differResult$logFC, c(-Inf, -input$logFCCutOff, input$logFCCutOff, Inf), c('Down','None','Up'))}
      }
      if(input$deaWay == 'edg'){
        output$volcanoPicture <- renderPlot(plotVolcano(differResult))
        shinyjs::show('volcanoPicture')
      }
      #maf的图借用火山图的output通道
      else if(input$deaWay == 'maf'){
        output$volcanoPicture <- renderPlot(forestPlot(deaResult(), 
                                                       pVal = input$pCutOff,
                                                       fdr = input$fdrCutOff))
        shinyjs::show('volcanoPicture')
      }
      else{shinyjs::hide('volcanoPicture')}
      if(input$deaWay != 'maf'){
        differResult[any(differResult$Status == 'Up', differResult$Status == 'Down'), , drop = FALSE]
      }
      else{differResult}
    }
    else{
      shinyjs::hide('volcanoPicture')
      NULL}
  })
  output$showDea <- renderUI({
    if(!is.null(deaScreen())){
      shinyjs::show('showDea')
      output$showScr <- DT::renderDT(deaScreen(), options=list(pageLength = 3))
      fluidRow(DT::DTOutput('showScr'),
               downloadButton('getDea', 'Save'))
    }
    else{shinyjs::hide('showDea')}
  })
  #富集分析
  enrichDeg <- reactive({list(up = symToEnt(deaScreen()[deaScreen()$Status == 'Up', ,drop = FALSE]),
                              down = symToEnt(deaScreen()[deaScreen()$Status == 'Down', ,drop = FALSE]))
  })
  goResult <- eventReactive(input$calGO, {
    goProgress <- Progress$new(min=1, max=4)
    on.exit(goProgress$close())
    goProgress$set(message = 'Calculation in eGO',
                   detail = 'This may take a while...',
                   value = 1)
    goUp <- eGO(enrichDeg()$up, input$enrichPCut, input$enrichQCut, input$enrichOnto)
    goProgress$set(value = 2)
    goDown <- eGO(enrichDeg()$down, input$enrichPCut, input$enrichQCut, input$enrichOnto)
    goProgress$set(value = 3)
    goAll <- eGO(rbind(enrichDeg()$up, enrichDeg()$down), input$enrichPCut, input$enrichQCut, input$enrichOnto)
    goProgress$set(value = 4)
    list(up = goUp, down = goDown, all = goAll)
  })
  keggResult <- eventReactive(input$calKegg, {
    keggProgress <- Progress$new(min=1, max=4)
    on.exit(keggProgress$close())
    keggProgress$set(message = 'Calculation in eKEGG',
                     detail = 'This may take a while...',
                     value = 1)
    keggUp <- eKegg(enrichDeg()$up, input$enrichPCut, input$enrichQCut)
    keggProgress$set(value = 2)
    keggDown <- eKegg(enrichDeg()$down, input$enrichPCut, input$enrichQCut)
    keggProgress$set(value = 3)
    keggAll <- eKegg(rbind(enrichDeg()$up, enrichDeg()$down), input$enrichPCut, input$enrichQCut)
    keggProgress$set(value = 4)
    list(up = keggUp, down = keggDown, all = keggAll)
  })
  output$goShow <- renderUI({
    if(!is.null(goResult())){
      shinyjs::show('goShow')
      output$goShowMatrixUp <- DT::renderDT(goResult()$up[[1]], options=list(pageLength = 2))
      output$goShowMatrixDown <- DT::renderDT(goResult()$down[[1]], options=list(pageLength = 2))
      output$goShowMatrixAll <- DT::renderDT(goResult()$all[[1]], options=list(pageLength = 2))
      output$goShowDotUp <- renderPlot(plotDot(goResult()$up[[2]], 'GO', input$goShowNum))
      output$goShowDotDown <- renderPlot(plotDot(goResult()$down[[2]], 'GO', input$goShowNum))
      output$goShowDotAll <- renderPlot(plotDot(goResult()$all[[2]], 'GO', input$goShowNum))
      output$goShowBarUp <- renderPlot(plotBar(goResult()$up[[2]], 'GO', input$goShowNum))
      output$goShowBarDown <- renderPlot(plotBar(goResult()$down[[2]], 'GO', input$goShowNum))
      output$goShowBarAll <- renderPlot(plotBar(goResult()$all[[2]], 'GO', input$goShowNum))
      fluidRow(h2('GO Result'),
               numericInput('goShowNum', 'Show Number', min = 3, max = 20, value = 5),
               tabsetPanel(
                 tabPanel('Up', DT::DTOutput('goShowMatrixUp'),
                          downloadButton('getGOMatrixUp', 'Save'),
                          plotOutput('goShowDotUp'), 
                          plotOutput('goShowBarUp')),
                 tabPanel('Down', DT::DTOutput('goShowMatrixDown'),
                          downloadButton('getGOMatrixDown', 'Save'),
                          plotOutput('goShowDotDown'),
                          plotOutput('goShowBarDown')),
                 tabPanel('All', DT::DTOutput('goShowMatrixAll'),
                          downloadButton('getGOMatrixAll', 'Save'),
                          plotOutput('goShowDotAll'),
                          plotOutput('goShowBarAll'))
               )
      )    
    }
    else{shinyjs::hide('goShow')}
  })
  output$keggShow <- renderUI({
    if(!is.null(keggResult())){
      shinyjs::show('keggShow')
      output$keggShowMatrixUp <- DT::renderDT(keggResult()$up[[1]], options=list(pageLength = 2))
      output$keggShowMatrixDown <- DT::renderDT(keggResult()$down[[1]], options=list(pageLength = 2))
      output$keggShowMatrixAll <- DT::renderDT(keggResult()$all[[1]], options=list(pageLength = 2))
      output$keggShowDotUp <- renderPlot(plotDot(keggResult()$up[[2]], 'KEGG', input$keggShowNum))
      output$keggShowDotDown <- renderPlot(plotDot(keggResult()$down[[2]], 'KEGG', input$keggShowNum))
      output$keggShowDotAll <- renderPlot(plotDot(keggResult()$all[[2]], 'KEGG', input$keggShowNum))
      output$keggShowBarUp <- renderPlot(plotBar(goResult()$up[[2]], 'KEGG', input$keggShowNum))
      output$keggShowBarDown <- renderPlot(plotBar(goResult()$down[[2]], 'KEGG', input$keggShowNum))
      output$keggShowBarAll <- renderPlot(plotBar(goResult()$all[[2]], 'KEGG', input$keggShowNum))
      fluidRow(h2('KEGG Result'),
               numericInput('keggShowNum', 'Show Number:', min = 3, max = 20, value = 5),
               tabsetPanel(
                 tabPanel('Up', DT::DTOutput('keggShowMatrixUp'),
                          downloadButton('getKEGGMatrixUp', 'Save'),
                          plotOutput('keggShowDotUp'),
                          plotOutput('keggShowBarUp')),
                 tabPanel('Down', DT::DTOutput('keggShowMatrixDown'),
                          downloadButton('getKEGGMatrixDown', 'Save'),
                          plotOutput('keggShowDotDown'),
                          plotOutput('keggShowBarDown')),
                 tabPanel('All', DT::DTOutput('keggShowMatrixAll'),
                          downloadButton('getKEGGMatrixAll', 'Save'),
                          plotOutput('keggShowDotAll'),
                          plotOutput('keggShowBarAll'))
               )
      )    
    }
    else{shinyjs::hide('keggShow')}
  })
  #相关性分析
  corResult <- eventReactive(input$calCor, {
    if(!is.null(input$corFactor1)){
      if(!is.null(input$corFactor2)){corCal(newClinalData(), input$corFactor1, input$corFactor2, input$corWay)}
      else{corCal(newClinalData(), input$corFactor1, input$corFactor1, input$corWay)}           
    }
  })
  corMatResult <- reactive({
    if(!is.null(corResult())){
      corResult()$mat
    }
    else{NULL}
  })
  corLsResult <- reactive({
    if(!is.null(corResult())){
      corRe <- corResult()$ls
      corRe <- corRe[abs(corRe$r) >= input$corrCut,,drop = FALSE]
      corRe <- corRe[corRe$p <= input$corpCut,,drop = FALSE]
      corRe
    }
    else{NULL}
  })
  output$corMatShow <- DT::renderDT(corMatResult())
  output$corLsShow <- DT::renderDT(corLsResult())
  output$corHeat <- renderPlot({
    corHeatMat <- corMatResult()
    corHeatMat <- corMatScreen(corHeatMat, corResult()$ls, input$corrCut, input$corpCut)
    plotHeat(corHeatMat)
  })
  output$corShow <- renderUI({
    if(!is.null(corMatResult())){
      shinyjs::show('corShow')
      fluidPage(h2('Correlation Matrix'),
                DT::DTOutput('corMatShow'),
                downloadButton('getCorMat', 'Save'),
                hr(),
                DT::DTOutput('corLsShow'),
                downloadButton('getCorLs', 'Save'),
                hr(),
                plotOutput('corHeat')
      )
    }
    else{shinyjs::hide('corShow')}
  })
  #maf动态UI
  output$mafIn <- renderUI({
    if(input$mafShowMode == 'sum'){NULL}
    else if(input$mafShowMode == 'self'){
      fluidRow(
        fileInput('topIn', 'Input sample data:', accept = '.csv'),
        fileInput('rightIn', 'Input genes data:', accept = '.csv'),
        selectInput('mafMutationType', 'Select a mutation Type:', 
                    choices = extrcactVariantType(mafObj()), multiple = TRUE),
        uiOutput('gpUI'),
        radioButtons('mafGeneOrPath', 'Select your interested object:', choices = c('Genes', 'Pathway'))
      )
    }
  })
  observe({
    if(!is.null(input$mafGeneOrPath)){
      if(input$mafGeneOrPath == 'Genes'){
        output$gpUI <- renderUI(textAreaInput('mafGene', 'Input the genes name:', placeholder = 'TP53,PRR11,SPP1...'))}
      else if(input$mafGeneOrPath == 'Pathway'){
        output$gpUI <- renderUI(fluidPage(fileInput('mafGmt', 'Input a .gmt file:', accept = '.gmt'),
                                          uiOutput('choosePathway')
        ))}
    } 
  })
  mafObj <- eventReactive(input$mafVisual, {
    readMafProgress <- Progress$new(min = 0, max = 3)
    on.exit(readMafProgress$close())
    readMafProgress$set(message = 'Read maf file.',
                        detail = 'This may take a while...', value = 1)
    originMaf <- read.maf(input$mafVisual$datapath, isTCGA = input$mafTCGAFlag)
    readMafProgress$set(value = 2)
    originMaf@data$Variant_Classification <- originMaf@data$Variant_Type
    readMafProgress$set(value = 3)
    originMaf
  })
  genesetTable <- eventReactive(input$mafGmt, {gmtToDataframe(input$mafGmt$datapath)})
  output$choosePathway <- renderUI(if(!is.null(genesetTable())){
    selectInput('mafKeyPathway', 'Choose pathways', 
                choices = unique(genesetTable()[,2]), multiple = TRUE)})
  sampleData <- reactive({
    if(!is.null(input$topIn)){readMatrix(input$topIn,FALSE)}
    else{NULL}
  })
  genesData <- reactive({
    if(!is.null(input$rightIn)){readMatrix(input$rightIn,FALSE)}
    else{NULL}
  })
  samOrder <- reactive({
    if(!is.null(sampleData())){
      sD <- sampleData()[order(sampleData()[,2], decreasing = TRUE),]
      sD[,1]
    }
    else{NULL}
  })
  observeEvent(input$mafPlot, {
    mafTemp <- mafObj()
    if(input$mafShowMode == 'sum'){
      output$mafSum <- renderPlot(plotmafSummary(maf = mafTemp, top = input$mafTop))
      output$mafSumWaterFall <- renderPlot(oncoplot(maf = mafTemp, top = input$mafTop), width = 1000)
      output$mafSumPathway <- renderCachedPlot({sump <- OncogenicPathways(mafTemp)
      output$selectSumPathway <- renderUI(selectInput('mafSumKeyPath', 'Choose a known pathway:', choices = rev(sump$Pathway)))}, 
      cacheKeyExpr = {sump <- OncogenicPathways(mafTemp)})
      output$mafSumKeyPathPlot <- renderPlot(if(!is.null(input$mafSumKeyPath)){PlotOncogenicPathways(mafTemp, pathways = input$mafSumKeyPath)})
      output$mafShow <- renderUI({fluidRow(plotOutput('mafSum'),
                                           plotOutput('mafSumWaterFall'),
                                           plotOutput('mafSumPathway'),
                                           uiOutput('selectSumPathway'),
                                           plotOutput('mafSumKeyPathPlot'))})
      
    }
    #TP53,CDH1,LPR2,MDN1
    else if(input$mafShowMode == 'self'){
      mafMutType <<- input$mafMutationType
      mafMutGene <- unlist(strsplit(input$mafGene, ','))
      if(input$mafGeneOrPath == 'Genes'){
        mafScreen <- subsetMaf(mafTemp, genes = mafMutGene, query = 'Variant_Classification %in% mafMutType')
        output$mafVaf <- renderPlot(plotVaf(maf = mafScreen))
        output$mafSomatic <- renderCachedPlot(somaticInteractions(mafScreen), cacheKeyExpr = {somaticInteractions(mafScreen)})
        output$mafGeneWaterFall <- renderPlot(oncoplot(mafScreen, top = input$mafTop, 
                                                       topBarData = sampleData(),
                                                       rightBarData = genesData(),
                                                       bgCol = 'white',
                                                       colors = c('SNP' = 'orange','blue','red', 'INS' = 'green', 'DEL' = 'purple'),
                                                       drawColBar = !is.null(input$topIn),
                                                       drawRowBar = !is.null(input$rightIn),
                                                       sampleOrder = samOrder()))
        output$mafShow <- renderUI({fluidRow(plotOutput('mafVaf'),
                                             plotOutput('mafSomatic'),
                                             plotOutput('mafGeneWaterFall'))})
      }
      else if(input$mafGeneOrPath == 'Pathway'){
        mafScreen <- subsetMaf(mafTemp, query = 'Variant_Classification %in% mafMutType')
        output$mafPathSummary <- renderCachedPlot(OncogenicPathways(mafScreen, pathways = genesetTable()), cacheKeyExpr = {OncogenicPathways(mafScreen, pathways = genesetTable())})
        output$mafPathWaterFall <- renderPlot(oncoplot(mafScreen, top = input$mafTop, 
                                                       topBarData = sampleData(),
                                                       bgCol = 'white',
                                                       drawColBar = !is.null(input$topIn),
                                                       drawRowBar = !is.null(input$rightIn),
                                                       sampleOrder = samOrder(),
                                                       pathways = genesetTable()[genesetTable()$Geneset_Name %in% input$mafKeyPathway, ,drop = FALSE]))
        output$mafShow <- renderUI(fluidRow(
          plotOutput('mafPathSummary'),
          plotOutput('mafPathWaterFall')
        ))
      }
    }
  })
  
  #文件获取区域
  {
    output$downloadTPM <- downloadHandler(
      filename = function() {paste(paste(strsplit(input$countsMatrix$name, '[.]')[[1]][1], tpmMatrix()$type, sep = '_'), '_TPM.csv', sep='') },
      content = function(file) {fwrite(tpmMatrix()$tpm, file, row.names = TRUE)})
    output$downloadCounts <- downloadHandler(
      filename = function() {paste(paste(strsplit(input$countsMatrix$name, '[.]')[[1]][1], tpmMatrix()$type, sep = '_'), '.csv', sep='') },
      content = function(file) {fwrite(tpmMatrix()$counts, file, row.names = TRUE)})
    output$downloadGSVA <- downloadHandler(
      filename = function() {paste(strsplit(input$countsMatrix$name, '[.]')[[1]][1], '_GSVA.csv', sep='')},
      content = function(file) {fwrite(gsvaMatrix(),file, row.names = TRUE)})
    output$getInterClinical <- downloadHandler(
      filename = function() { paste(input$cancerType[[1]],'internal_Clinical.csv', sep='_') },
      content = function(file) {fwrite(chooseInterClinical(), file, row.names = TRUE)}) 
    output$getNewClinical <- downloadHandler(
      filename = function() { paste(strsplit(input$clinicalData$name, '[.]')[[1]][1], '_new.csv', sep='') },
      content = function(file) {fwrite(newClinalData(), file, row.names = TRUE)})
    output$getDea <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_dea.csv', sep='') },
      content = function(file) {fwrite(deaScreen(), file, row.names = TRUE)})
    output$getGOMatrixUp <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_Up.csv', sep='') },
      content = function(file) {fwrite(goResult()$up[[1]], file, row.names = TRUE)})
    output$getGOMatrixDown <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_Down.csv', sep='') },
      content = function(file) {fwrite.csv(goResult()$down[[1]], file, row.names = TRUE)})
    output$getGOMatrixAll <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_All.csv', sep='') },
      content = function(file) {fwrite(goResult()$all[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixUp <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Up.csv', sep='') },
      content = function(file) {fwrite(keggResult()$up[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixDown <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Down.csv', sep='') },
      content = function(file) {fwrite(keggResult()$up[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixAll <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_All.csv', sep='') },
      content = function(file) {fwrite(keggResult()$all[[1]], file, row.names = TRUE)})
    output$getCorMat <- downloadHandler(
      filename = function() { paste(paste(input$corFactor1[[1]],strsplit(input$clinicalData$name, '[.]')[[1]][1],sep = '_'), '_cor.csv', sep='') },
      content = function(file) {fwrite(corMatResult(), file, row.names = TRUE)})
    output$getCorLs <- downloadHandler(
      filename = function() { paste(paste(input$corFactor1[[1]],strsplit(input$clinicalData$name, '[.]')[[1]][1],sep = '_'), '_corlist.csv', sep='') },
      content = function(file) {fwrite(corLsResult(), file, row.names = TRUE)})
  }
}
