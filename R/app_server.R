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
  #cox参考组的UI
  observe({
    fL <- input$surFactor
    if(is.null(fL)){quoteExpr<-quote(NULL)}
    else{
      quoteExpr <- quote(fluidPage)
      for(t in input$surFactor){
        if(!is.numeric(newClinalData()[[t]])){
          quoteExpr <- append(quoteExpr,
                              bquote(
                                selectInput(paste0(.(t),'Ref'), paste(.(t), 'reference'), choices = levels(as.factor(newClinalData()[[.(t)]])))
                              )
          )
        }
      }
      quoteExpr <- as.call(quoteExpr)
      if(quoteExpr == quote(fluidPage)){quoteExpr <- quote(NULL)}
    }
    output$coxRefUI <- renderUI(quoteExpr, quoted = TRUE)
  })
  observe({
    if(input$surWay == 'sur'){shinyjs::hide('coxRefUI')}
    else{shinyjs::show('coxRefUI')}
  })
  output$surGroupUI <- renderUI({
    if(!is.null(input$groupInfo)){
      if(input$groupInfo == 'No group'){NULL}
      else{selectInput('surGroup', 'Choose a group:', choices = levels(factor(newClinalData()[[input$groupInfo]])))} 
    }
  })
  #生存分析计算矩阵
  surClin <- reactive({
    factorList <- input$surFactor
    if(!is.null(newClinalData())){
      if(!is.numeric(newClinalData()[,input$surSta])){
        factorSta <- tolower(levels(as.factor(newClinalData()[,input$surSta])))
        if(any(all(factorSta == c('alive','dead')), all(factorSta == c('dead','alive')))){
          shinyjs::hide('deadEvent')
          tempSurClin <- toSurStatus(newClinalData(), input$surSta, 'dead')
        }
        else{
          output$chooseEvent <- renderUI(selectInput('deadEvent', 'Input the dead event: ', choices = factorSta))
          shinyjs::show('deadEvent')
          tempSurClin <- toSurStatus(newClinalData(), input$surSta, input$deadEvent)
        }
      }
      else{shinyjs::hide('deadEvent')
        tempSurClin <- newClinalData()}
      #设置因子水平
      for (f in factorList){
        facRef <- input[[paste0(f,'Ref')]]
        if(!is.null(facRef)){
          tempSurClin[[f]] <- relevel(factor(tempSurClin[[f]]), ref = facRef)
        }
      }
      #选组
      if(!is.null(input$surGroup)){tempSurClin <- tempSurClin[tempSurClin[[input$groupInfo]] == input$surGroup, , drop = FALSE]}
      tempSurClin
    }
  })
  output$showtest <- DT::renderDT({
    tryCatch({surClin()[unique(c(input$surTime, input$surSta, input$surFactor))]}, error = function(x){NULL})
  },  options=list(pageLength = 1))
  shinyjs::hide('surPlotSize')
  shinyjs::hide('forestSize')
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
      sCoxResult <- getCoxTable(singleCox(surClin(), input$surTime, input$surSta, input$surFactor[[1]]))
      for(fac in input$surFactor){
        if(fac == input$surFactor[[1]]){next}
        sCox <- getCoxTable(singleCox(surClin(), input$surTime, input$surSta, fac))
        sCoxResult <- rbind(sCoxResult, sCox)
      }
      sCoxResult
    }
    else if(input$surWay == 'multipleCox'){
      getCoxTable(multipleCox(surClin(), input$surTime, input$surSta, input$surFactor))
    }
  })
  output$surShow <- renderUI({
    if(!is.null(surResult())){
      if(isolate(input$surWay) == 'sur'){
        shinyjs::hide('forestSize')
        output$surPlot <- renderPlot(surResult(),
                                     width = ceiling(sqrt(length(isolate(input$surFactor))))*input$surPlotSize,
                                     height = round(sqrt(length(isolate(input$surFactor))))*input$surPlotSize)
        shinyjs::show('surPlotSize')
        plotOutput('surPlot', height = round(sqrt(length(isolate(input$surFactor))))*input$surPlotSize)
      }
      else{
        shinyjs::hide('surPlotSize')
        output$singleCoxShow <- DT::renderDT(surResult())
        output$singleCoxForest <- renderPlot(diyForest(surResult()),
                                             width = input$forestSize,
                                             height = 100+20*nrow(surResult()))
        shinyjs::show('forestSize')
        fluidRow(plotOutput('singleCoxForest', height = 100+20*nrow(surResult())),
                 br(),
                 DT::DTOutput('singleCoxShow'))
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
      shinyjs::hide('groupCutOff')
      shinyjs::show('exGroup')
      shinyjs::show('ctGroup')
    }
    else{
      shinyjs::show('groupCutOff')
      shinyjs::hide('exGroup')
      shinyjs::hide('ctGroup')
    }
  })
  deaObj <- reactive({if(input$deaWay == 'maf'){read.maf(input$deaData$datapath, isTCGA = input$deaTCGAFlag)}
    else{readMatrix(input$deaData)}
  })
  #差异分析的condit分组
  deaGroupList <- reactive({
    if(!is.null(input$groupInfo)){
      if(input$groupInfo == 'No group'){NULL}
      else{levels(as.factor(newClinalData()[[input$groupInfo]]))}
    }
    else{NULL}
  })
  output$deaGroupUI <- renderUI({
    if(!is.null(input$groupInfo)){
      if(input$groupInfo == 'No group'){NULL}
      else{selectInput('deaOutGroup', 'Choose a group to out put:', choices = names(deaResult()))}
    }
    else{NULL}
  })
  #差异分析
  deaResult <- eventReactive(input$calDiffer,{
    if(!is.null(deaObj())){
      #无分组计算
      if(input$groupInfo == 'No group'){
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
      #分组计算
      else{
        resultDea <- list()
        clin <- newClinalData()
        for(q in deaGroupList()){
          groupCondit <- clin[clin[input$groupInfo] == q, ][input$deaFactor]
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
            resultDea[[q]] <- transferID(deaEdgeR(deaObj(), groupCondit, ex, ct))
          }
          else if(input$deaWay == 'lm'){
            resultDea[[q]] <- deaLimma(deaObj(), groupCondit, ex, ct)
          }
          else if(input$deaWay == 'maf'){
            resultDea[[q]] <- mafDiffer(deaObj(), groupCondit, ex, ct)
          }
          else{print('?')}
        }
        resultDea
      }
    }   
  })
  deaScreen <- reactive({
    wayDea <- isolate(input$deaWay)
    if(!is.null(deaResult())){
      if(!is.data.frame(deaResult())){
        if(all(names(deaResult()) == c('results', 'SampleSummary'))){deaToPrint <- deaResult()}
        else{deaToPrint <- deaResult()[[input$deaOutGroup]]}
      }
      else{deaToPrint <- deaResult()}
      if(wayDea == 'edg'){
        differResult <- deaToPrint[deaToPrint$PValue <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$FDR <= input$fdrCutOff, , drop = FALSE]
      }
      else if(wayDea == 'lm'){
        differResult <- deaToPrint[deaToPrint$P.Value <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$adj.P.Val <= input$fdrCutOff, , drop = FALSE]
      }
      else if(wayDea == 'maf'){
        differResult <- deaToPrint$results
        differResult <- differResult[differResult$pval <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$adjPval <= input$fdrCutOff, , drop = FALSE]
      }
      #logFC过滤没maf差异的份
      if(wayDea != 'maf'){
        if(input$logFCCutOff == 0)
        {differResult$Status <- cut(differResult$logFC, c(-Inf, input$logFCCutOff, Inf), c('Down','Up'))}
        else
        {differResult$Status <- cut(differResult$logFC, c(-Inf, -input$logFCCutOff, input$logFCCutOff, Inf), c('Down','None','Up'))}
      }
      if(wayDea == 'edg'){
        output$volcanoPicture <- renderPlot(plotVolcano(differResult))
        shinyjs::show('volcanoPicture')
      }
      #maf的图借用火山图的output通道
      else if(wayDea == 'maf'){
        output$volcanoPicture <- renderPlot(forestPlot(deaToPrint, 
                                                       pVal = input$pCutOff,
                                                       fdr = input$fdrCutOff))
        shinyjs::show('volcanoPicture')
      }
      else{shinyjs::hide('volcanoPicture')}
      if(wayDea != 'maf'){
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
  #相关性分析动态UI
  output$corActiveUI <- renderUI({
    if(!is.null(input$groupInfo)){
      #分组UI
      if(input$groupInfo != 'No group'){
        allGroup <- table(factor(newClinalData()[[input$groupInfo]]))
        fluidPage(
          selectInput('groupCorFactor', 'Single factor for group correlation:', choices = c('All', input$corFactor1)),
          selectInput('corGroup', 'Show Group:', choices = names(allGroup[allGroup > 3]))
        )
      }
      #不分组UI
      else{NULL}
    }
    else{NULL}
  })
  observeEvent(input$groupCorFactor,{
    if(input$groupCorFactor == 'All'){shinyjs::show('corGroup')}
    else{shinyjs::hide('corGroup')}
  })
  #相关性分析
  corResult <- eventReactive(input$calCor, {
    if(!is.null(input$corFactor1)){
      if(!is.null(input$corFactor2)){corCal(newClinalData(), input$corFactor1, input$corFactor2, input$corWay, input$groupInfo, input$groupCorFactor)}
      else{corCal(newClinalData(), input$corFactor1, input$corFactor1, input$corWay, input$groupInfo, input$groupCorFactor)}            
    }
  })
  corMatResult <- reactive({
    if(!is.null(corResult())){
      if(!is.null(input$groupCorFactor)){
        if('All'%in%names(corResult())){
          #分组全比较
          if(input$groupCorFactor == 'All'){
            corResult()[[input$groupCorFactor]]$mat[[input$corGroup]]
          }
          #分组单因子比较
          else{corResult()[[input$groupCorFactor]]$mat}
        }
        #不分组
        else{corResult()$mat}
      }
      #初始状态不分组
      else{corResult()$mat}
    }
    else{NULL}
  })
  corLsResult <- reactive({
    if(!is.null(corResult())){
      if(isolate(input$groupInfo) =='No group'){corRe <- corResult()$ls}
      else{corRe <- corResult()[[input$groupCorFactor]]$ls}
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
    if('All'%in%names(corResult())){
      corScr <- corResult()[[input$groupCorFactor]]$ls
      if(length(colnames(corScr)) == 5){corScr <- corScr[corScr[[input$groupInfo]] == input$corGroup, , drop =FALSE]}
    }
    else{corScr <- corResult()$ls}
    corHeatMat <- corMatScreen(corHeatMat, corScr, input$corrCut, input$corpCut)
    plotHeat(corHeatMat)
  })
  output$corShow <- renderUI({
    if(!is.null(corMatResult())){
      shinyjs::show('corShow')
      fluidPage(
        h2('Correlation Matrix'),
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
