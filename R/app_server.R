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
#' @import R.utils
#' @noRd
app_server <- function( input, output, session ) {
  options(shiny.maxRequestSize=-1) # Remove limit of upload
  options(shiny.deprecation.messages=FALSE)
  options(warn =-1)
  R.utils::setOption("clusterProfiler.download.method", "auto")
  #下载TCGA文件
  observe({
    allFlag <- input$all
    if (allFlag){updateCheckboxGroupInput(session, 'cancerType','Choose the cancer type:',inline = T, 
                                          choices = allCancer, selected =  allCancer)}
    else{updateCheckboxGroupInput(session, 'cancerType','Choose the cancer type:',inline = T, 
                                  choices = allCancer)}
  })
  #获取内部临床信息文件
  chooseInterClinical <- reactive({
    if(!is.null(input$cancerType)){
      output$titleInterClinical <- renderUI({h3('Clinical Data')})
      shinyjs::show(id = 'titleInterClinical')}
    else{
      shinyjs::hide(id = 'titleInterClinical')}
    choInterClinical <- interClinical[which(interClinical$type %in% input$cancerType),]
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
      h3('Gene Set')
  })
  output$genesetMatrixShow <- DT::renderDataTable(genesetMatrix())
  #tpm,tpmMatrix改为tpm列表，包含有counts矩阵和分类类型
  tpmMatrix <- reactive({
    countsData <- readMatrix(input$countsMatrix)
    tpmIn <- NULL
    if(!is.null(countsData)){
      if(input$tpmInFlag){tpmIn <- countsData}
      else{tpmIn <- NULL}
    }
    if (!is.null(countsData)&is.null(tpmIn)){
      output$tpmDownload <- renderUI({
        fluidRow(column(width = 4, downloadButton('downloadTPM', 'Save TPM')),
                 column(width = 4, downloadButton('downloadCounts', 'Save Counts')))
      })
      output$titleTpm <- renderUI({h2('TPM')})
      if(input$sampleType == 'all'){
        matrixList <- list(tpm = countsToTPM(countsData, transID = input$transidFlag), counts = countsData, type = 'all')}
      else if(input$sampleType == 'tumor'){
        cts <- classifyCounts(countsData)$tm
        matrixList <- list(tpm = countsToTPM(cts, transID = input$transidFlag), counts = cts, type = 'tumor')
      }
      else if(input$sampleType ==  'normal'){
        cts <- classifyCounts(countsData)$nr
        matrixList <- list(tpm = countsToTPM(cts, transID = input$transidFlag), counts = cts, type = 'norm')
      }
    }
    if (!is.null(tpmIn)){
      output$tpmDownload <- renderUI({downloadButton('downloadTPM', 'Save TPM')})
      output$titleTpm <- renderUI({h2('TPM')})
      if(input$sampleType == 'all'){
        matrixList <- list(tpm = tpmIn, counts = NULL, type = 'all')}
      else if(input$sampleType == 'tumor'){
        matrixList <- list(tpm = classifyCounts(tpmIn)$tm, counts = NULL, type = 'tumor')
      }
      else if(input$sampleType ==  'normal'){
        matrixList <- list(tpm = classifyCounts(tpmIn)$nr, counts = NULL, type = 'norm')
      }
    }
    if(!is.null(input$countsMatrix)){matrixList}
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
      if(all(quoteExpr == quote(fluidPage))){quoteExpr <- quote(NULL)}
      else{quoteExpr <- as.call(quoteExpr)}
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
  coxForest <- reactive(diyForest(surResult()))
  output$surShow <- renderUI({
    if(!is.null(surResult())){
      if(isolate(input$surWay) == 'sur'){
        shinyjs::hide('forestSize')
        output$surPlot <- renderPlot(surResult(),
                                     width = ceiling(sqrt(length(isolate(input$surFactor))))*input$surPlotSize,
                                     height = round(sqrt(length(isolate(input$surFactor))))*input$surPlotSize)
        shinyjs::show('surPlotSize')
        fluidRow(
          plotOutput('surPlot', height = round(sqrt(length(isolate(input$surFactor))))*input$surPlotSize),
          downloadButton('surPlotSave', 'Save(.pdf)') 
        )
      }
      else{
        shinyjs::hide('surPlotSize')
        output$singleCoxShow <- DT::renderDT(surResult())
        output$singleCoxForest <- renderPlot(coxForest(),
                                             width = input$forestSize,
                                             height = 100+20*nrow(surResult()))
        shinyjs::show('forestSize')
        fluidRow(plotOutput('singleCoxForest', height = 100+20*nrow(surResult())),
                 br(),
                 downloadButton('coxForestSave', 'Save(.pdf)'),
                 hr(),
                 DT::DTOutput('singleCoxShow'),
                 br(),
                 downloadButton('coxResultSave', 'Save')
        )
      }
    }
    else{NULL}
  })
  #差异分析动态UI
  observe(if(input$calDiffer == 0){shinyjs::hide('volcanoSave')})#初始化隐藏火山图下载按钮
  observe({if(input$deaWay != 'edg'){shinyjs::hide('deaNorm')}
           else {shinyjs::show('deaNorm')}})#基因表达差异分析时才有标准化方法选择
  observe({
    temp <- tryCatch(newClinalData()[,input$deaFactor], error = function(x){NULL})
    if(!is.numeric(temp)){
      factorDea <- levels(as.factor(temp))
      output$chooseExperience <- renderUI(selectInput('exGroup','Experiment group', choices = factorDea))
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
  deaObj <- reactive({if(input$deaWay == 'maf'){
    prReadmaf <- Progress$new(min=0, max=2)
    on.exit(prReadmaf$close())
    prReadmaf$set(message = 'Reading maf file',detail = 'This may take a while...', value = 1)
    mafToDea <- read.maf(input$deaData$datapath, isTCGA = input$deaTCGAFlag)
    prReadmaf$set(value = 2)
    mafToDea
  }
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
          transferID(deaEdgeR(deaObj(), groupCondit, ex, ct, input$deaNorm))
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
            resultDea[[q]] <- transferID(deaEdgeR(deaObj(), groupCondit, ex, ct, input$deaNorm))
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
        #非分组
        if(all(names(deaResult()) == c('results', 'SampleSummary'))){deaToPrint <- deaResult()}
        #分组
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
        shinyjs::hide('volcanoPicture')
        shinyjs::hide('volcanoSave')
      }
      else if(wayDea == 'maf'){
        differResult <- deaToPrint$results
        differResult <- differResult[differResult$pval <= input$pCutOff, , drop = FALSE]
        differResult <- differResult[differResult$adjPval <= input$fdrCutOff, , drop = FALSE]
        deaToPrint$results <- differResult
      }
      #logFC过滤没maf差异的份
      if(wayDea != 'maf'){
        if(input$logFCCutOff == 0)
        {differResult$Status <- cut(differResult$logFC, c(-Inf, input$logFCCutOff, Inf), c('Down','Up'))}
        else
        {differResult$Status <- cut(differResult$logFC, c(-Inf, -input$logFCCutOff, input$logFCCutOff, Inf), c('Down','None','Up'))}
      }
      if(wayDea == 'edg'){
        deaVolcano <- reactive(plotVolcano(differResult))
        output$volcanoPicture <- renderPlot(deaVolcano())
        output$volcanoPictureUI <- renderUI(plotOutput('volcanoPicture'))
        output$volcanoSave <- downloadHandler(filename = function(){paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_volcano.pdf', sep='')},
                                              content = function(file){ggsave(file, plot = deaVolcano(), device = 'pdf', dpi = 600)},
                                              contentType = 'pdf')   #火山图输出在这儿
        shinyjs::show('volcanoPicture')
        shinyjs::show('volcanoSave')
      }
      #maf的图借用火山图的output通道
      else if(wayDea == 'maf'){
        deaToPrint$results <- differResult #改完参数过滤后再画图
        plNum <- nrow(deaToPrint$results)  #输出基因的数量
        output$volcanoPicture <- renderPlot(forestPlot(deaToPrint,pVal = 1.1))
        output$volcanoPictureUI <- renderUI(plotOutput('volcanoPicture', width = '8in', height = paste0(4.65+0.15*plNum, 'in')))
        output$volcanoSave <- downloadHandler(filename = function(){paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_mafdeaforest.pdf', sep='')},
                                              content = function(file){pdf(file, width = 8, height = 4.65+0.15*plNum)
                                                forestPlot(deaToPrint, pVal = 1.1)
                                                dev.off()},
                                              contentType = 'pdf')  #森林图输出在这
        shinyjs::show('volcanoPicture')
        shinyjs::show('volcanoSave')
      }
      else{shinyjs::hide('volcanoPicture')}
      if(wayDea != 'maf'){
        differResult[differResult$Status == 'Up'|differResult$Status == 'Down', , drop = FALSE]
      }
      else{differResult}
    }
    else{
      shinyjs::hide('volcanoPicture')
      shinyjs::hide('volcanoSave')
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
  #富集分析动态UI
  observe({if(input$enrichGene == 'DEG'){shinyjs::hide('enrichGeneInput')}
    else {shinyjs::show('enrichGeneInput')}})#自定义富集基因在计算差异基因时隐藏
  #富集分析
  enrichDeg <- reactive({
    #利用以上所得DEG
    if(input$enrichGene == 'DEG'){
      list(up = symToEnt(deaScreen()[deaScreen()$Status == 'Up', ,drop = FALSE]),
           down = symToEnt(deaScreen()[deaScreen()$Status == 'Down', ,drop = FALSE]))
    }
    #自输入基因
    else{
      list(self = bitr(unlist(strsplit(input$enrichGeneInput, split = ',')), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'))
    }
  })
  goResult <- eventReactive(input$calGO, {
    #利用以上所得DEG
    if(length(enrichDeg()) == 2){
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
    }
    #自输入基因
    else if(length(enrichDeg()) == 1){
      goProgress <- Progress$new(min=1, max=2)
      on.exit(goProgress$close())
      goProgress$set(message = 'Calculation in eGO',
                     detail = 'This may take a while...',
                     value = 1)
      eGOresult <- eGO(enrichDeg()$self, input$enrichPCut, input$enrichQCut, input$enrichOnto)
      goProgress$set(value = 2)
      return(list(self = eGOresult))
    }
    else{return(NULL)}
  })
  keggResult <- eventReactive(input$calKegg, {
    #利用以上所得DEG
    if(length(enrichDeg()) == 2){
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
    }
    #自输入基因
    else if(length(enrichDeg()) == 1){
      keggProgress <- Progress$new(min=1, max=2)
      on.exit(keggProgress$close())
      keggProgress$set(message = 'Calculation in eKEGG',
                       detail = 'This may take a while...',
                       value = 1)
      eKEGGresult <- eKegg(enrichDeg()$self, input$enrichPCut, input$enrichQCut)
      keggProgress$set(value = 2)
      return(list(self = eKEGGresult))
    }
    else{return(NULL)}
  })
  
  goDotUp <- reactive(plotDot(goResult()$up[[2]], 'GO', input$goShowNum))
  goBarUp <- reactive(plotBar(goResult()$up[[2]], 'GO', input$goShowNum))
  
  goDotDown <- reactive(plotDot(goResult()$down[[2]], 'GO', input$goShowNum))
  goBarDown <- reactive(plotBar(goResult()$down[[2]], 'GO', input$goShowNum))
  
  goDotAll <- reactive(plotDot(goResult()$all[[2]], 'GO', input$goShowNum))
  goBarAll <- reactive(plotBar(goResult()$all[[2]], 'GO', input$goShowNum))
  
  keggDotUp <- reactive(plotDot(keggResult()$up[[2]], 'KEGG', input$keggShowNum))
  keggBarUp <- reactive(plotBar(keggResult()$up[[2]], 'KEGG', input$keggShowNum))
  
  keggDotDown <- reactive(plotDot(keggResult()$down[[2]], 'KEGG', input$keggShowNum))
  keggBarDown <- reactive(plotBar(keggResult()$down[[2]], 'KEGG', input$keggShowNum))
  
  keggDotAll <- reactive(plotDot(keggResult()$all[[2]], 'KEGG', input$keggShowNum))
  keggBarAll <- reactive(plotBar(keggResult()$all[[2]], 'KEGG', input$keggShowNum))
  
  #利用自输入基因画的图
  goDotSelf <- reactive(plotDot(goResult()$self[[2]], 'GO', input$goShowNum))
  goBarSelf <- reactive(plotBar(goResult()$self[[2]], 'GO', input$goShowNum))
  
  keggDotSelf <- reactive(plotDot(keggResult()$self[[2]], 'KEGG', input$keggShowNum))
  keggBarSelf <- reactive(plotBar(keggResult()$self[[2]], 'KEGG', input$keggShowNum))
  
  goHeight <- reactive({
    if(is.null(input$goShowNum)){n <- 5}
    else{n <- input$goShowNum}
    if(input$enrichOnto == 'ALL'){120+3*n*20}
    else{if(n<=8){280}
      else{280 + (n-8)*20}}  
  })  #GO图高度
  keggHeight <- reactive({
    if(is.null(input$keggShowNum)){n <- 5}
    else{n <- input$keggShowNum}
    if(n<=8){280}
    else{280 + (n-8)*20}}) #KEGG图高度
  output$goHead <- renderUI({
    if(!is.null(goResult())){
      fluidRow(h2('GO Enrichment Result'),
               numericInput('goShowNum', 'Show Number:', min = 3, max = 20, value = 5))
    }
    else{NULL}
  }) #GO结果设置及标题
  output$keggHead <- renderUI({
    if(!is.null(keggResult())){
      fluidRow(h2('KEGG Enrichment Result'),
               numericInput('keggShowNum', 'Show Number:', min = 3, max = 20, value = 5))
    }
    else{NULL}
  }) #KEGG结果设置及标题
  output$goShow <- renderUI({
    if(!is.null(goResult())){
      shinyjs::show('goShow')
      #DEG
      if(isolate(input$enrichGene) == 'DEG'){
        output$goShowMatrixUp <- DT::renderDT(goResult()$up[[1]], options=list(pageLength = 2))
        output$goShowMatrixDown <- DT::renderDT(goResult()$down[[1]], options=list(pageLength = 2))
        output$goShowMatrixAll <- DT::renderDT(goResult()$all[[1]], options=list(pageLength = 2))
        
        output$goShowDotUp <- renderPlot(goDotUp()) 
        output$goShowDotDown <- renderPlot(goDotDown())
        output$goShowDotAll <- renderPlot(goDotAll())
        
        output$goShowBarUp <- renderPlot(goBarUp()) 
        output$goShowBarDown <- renderPlot(goBarDown())
        output$goShowBarAll <- renderPlot(goBarAll())
        
        fluidRow(
          tabsetPanel(
            tabPanel('Up regulation', DT::DTOutput('goShowMatrixUp'),
                     downloadButton('getGOMatrixUp', 'Save'),
                     plotOutput('goShowDotUp', width = 800, height = goHeight()),
                     downloadButton('goDotUpSave', 'Save(.pdf)'),
                     plotOutput('goShowBarUp', width = 800, height = goHeight()),
                     downloadButton('goBarUpSave', 'Save(.pdf)')),
            tabPanel('Down regulation', DT::DTOutput('goShowMatrixDown'),
                     downloadButton('getGOMatrixDown', 'Save'),
                     plotOutput('goShowDotDown', width = 800, height = goHeight()),
                     downloadButton('goDotDownSave', 'Save(.pdf)'),
                     plotOutput('goShowBarDown', width = 800, height = goHeight()),
                     downloadButton('goBarDownSave', 'Save(.pdf)'),),
            tabPanel('All', DT::DTOutput('goShowMatrixAll'),
                     downloadButton('getGOMatrixAll', 'Save'),
                     plotOutput('goShowDotAll', width = 800, height = goHeight()),
                     downloadButton('goDotAllSave', 'Save(.pdf)'),
                     plotOutput('goShowBarAll', width = 800, height = goHeight()),
                     downloadButton('goBarAllSave', 'Save(.pdf)'),)
          )
        )
      }
      #自输入基因
      else{
        output$goShowMatrixSelf <- DT::renderDT(goResult()$self[[1]], options=list(pageLength = 2))
        output$goShowDotSelf <- renderPlot(goDotSelf())
        output$goShowBarSelf <- renderPlot(goBarSelf())
        fluidRow(
          DT::DTOutput('goShowMatrixSelf'),
          downloadButton('getGOMatrixSelf', 'Save'),
          plotOutput('goShowDotSelf', width = 800, height = goHeight()),
          downloadButton('goDotSelfSave', 'Save(.pdf)'),
          plotOutput('goShowBarSelf', width = 800, height = goHeight()),
          downloadButton('goBarSelfSave', 'Save(.pdf)')
        )
      }
    }
    else{shinyjs::hide('goShow')}
  })
  output$keggShow <- renderUI({
    if(!is.null(keggResult())){
      shinyjs::show('keggShow')
      #DEG
      if(isolate(input$enrichGene) == 'DEG'){
        output$keggShowMatrixUp <- DT::renderDT(keggResult()$up[[1]], options=list(pageLength = 2))
        output$keggShowMatrixDown <- DT::renderDT(keggResult()$down[[1]], options=list(pageLength = 2))
        output$keggShowMatrixAll <- DT::renderDT(keggResult()$all[[1]], options=list(pageLength = 2))
        
        output$keggShowDotUp <- renderPlot(keggDotUp()) 
        output$keggShowDotDown <- renderPlot(keggDotDown())
        output$keggShowDotAll <- renderPlot(keggDotAll())
        
        output$keggShowBarUp <- renderPlot(keggBarUp()) 
        output$keggShowBarDown <- renderPlot(keggBarDown())
        output$keggShowBarAll <- renderPlot(keggBarAll())
        
        fluidRow(tabsetPanel(
          tabPanel('Up', DT::DTOutput('keggShowMatrixUp'),
                   downloadButton('getKEGGMatrixUp', 'Save'),
                   plotOutput('keggShowDotUp', height = keggHeight()),
                   downloadButton('keggDotUpSave', 'Save(.pdf)'),
                   plotOutput('keggShowBarUp', height = keggHeight()),
                   downloadButton('keggBarUpSave', 'Save(.pdf)')),
          tabPanel('Down', DT::DTOutput('keggShowMatrixDown'),
                   downloadButton('getKEGGMatrixDown', 'Save'),
                   plotOutput('keggShowDotDown', height = keggHeight()),
                   downloadButton('keggDotDownSave', 'Save(.pdf)'),
                   plotOutput('keggShowBarDown', height = keggHeight()),
                   downloadButton('keggBarDownSave', 'Save(.pdf)')),
          tabPanel('All', DT::DTOutput('keggShowMatrixAll'),
                   downloadButton('getKEGGMatrixAll', 'Save'),
                   plotOutput('keggShowDotAll', height = keggHeight()),
                   downloadButton('keggDotAllSave', 'Save(.pdf)'),
                   plotOutput('keggShowBarAll', height = keggHeight()),
                   downloadButton('keggBarAllSave', 'Save(.pdf)'))
          )
        )
      }
      #自输入基因
      else{
        output$keggShowMatrixSelf <- DT::renderDT(keggResult()$self[[1]], options=list(pageLength = 2))
        output$keggShowDotSelf <- renderPlot(keggDotSelf())
        output$keggShowBarSelf <- renderPlot(keggBarSelf())
        fluidRow(
          DT::DTOutput('keggShowMatrixSelf'),
          downloadButton('getKEGGMatrixSelf', 'Save'),
          plotOutput('keggShowDotSelf', width = 800, height = keggHeight()),
          downloadButton('keggDotSelfSave', 'Save(.pdf)'),
          plotOutput('keggShowBarSelf', width = 800, height = keggHeight()),
          downloadButton('keggBarSelfSave', 'Save(.pdf)')
        )
      }
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
  corHeatMap <- reactive({
    corHeatMat <- corMatResult()
    if('All'%in%names(corResult())){
      corScr <- corResult()[[input$groupCorFactor]]$ls
      if(length(colnames(corScr)) == 5){corScr <- corScr[corScr[[input$groupInfo]] == input$corGroup, , drop =FALSE]}
    }
    else{corScr <- corResult()$ls}
    corHeatMat <- corMatScreen(corHeatMat, corScr, input$corrCut, input$corpCut)
    plotHeat(corHeatMat)
  })
  heatMapWidth <- reactive({
    if(is.null(input$heatSize)){n <- 20}
    else{n <- input$heatSize}
    150+nrow(corMatResult())*n})
  heatMapHeight <- reactive({
    if(is.null(input$heatSize)){n <- 20}
    else{n <- input$heatSize}
    150+ncol(corMatResult())*n})
  output$corHeat <- renderPlot(corHeatMap())
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
        hr()
      )
    }
    else{shinyjs::hide('corShow')}
  })
  output$corHeatShow <- renderUI({
    if(!is.null(corMatResult())){
      shinyjs::show('corHeatShow')
      fluidPage(
        plotOutput('corHeat', width = heatMapWidth(), height = heatMapHeight()),
        downloadButton('corHeatSave', 'Save(.pdf)')
      )
    }
    else{shinyjs::hide('corHeatShow')}
  })
  #maf动态UI
  output$mafIn <- renderUI({
    if(input$mafShowMode == 'sum'){NULL}
    else if(input$mafShowMode == 'self'){
      fluidRow(
        fileInput('topIn', 'Input sample data(.csv):', accept = '.csv'),
        fileInput('rightIn', 'Input genes data(.csv):', accept = '.csv'),
        selectInput('mafMutationType', 'Select mutation Type:', 
                    choices = extrcactVariantType(mafObj()), multiple = TRUE),
        uiOutput('gpUI'),
        radioButtons('mafGeneOrPath', 'Select your interested object:', choices = c('Genes', 'Pathway(Gene set)' = 'Pathway'))
      )
    }
  })
  observe({
    if(!is.null(input$mafGeneOrPath)){
      if(input$mafGeneOrPath == 'Genes'){
        output$gpUI <- renderUI(textAreaInput('mafGene', 'Input the genes name:', placeholder = 'TP53,PRR11,SPP1...'))}
      else if(input$mafGeneOrPath == 'Pathway'){
        output$gpUI <- renderUI(fluidPage(fileInput('mafGmt', 'Input gene set file(.gmt):', accept = '.gmt'),
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
    originMaf@data$Variant_Classification <- originMaf@data$VARIANT_CLASS
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
                                                output$selectSumPathway <- renderUI(selectInput('mafSumKeyPath', 'Choose a pathway:', choices = rev(sump$Pathway)))}, 
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
      mafMutGene <<- unlist(strsplit(input$mafGene, ','))
      if(input$mafGeneOrPath == 'Genes'){
        mafScreen <<- subsetMaf(mafTemp, genes = mafMutGene, query = 'VARIANT_CLASS %in% mafMutType')
        output$mafVaf <- renderPlot(plotVaf(maf = mafScreen))
        output$mafSomatic <- renderCachedPlot(somaticInteractions(mafScreen, fontSize = 0.6), cacheKeyExpr = {somaticInteractions(mafScreen, fontSize = 0.6)})
        output$mafGeneWaterFall <- renderPlot(oncoplot(mafScreen, top = input$mafTop, 
                                                       topBarData = sampleData(),
                                                       rightBarData = genesData(),
                                                       bgCol = '#EEEEEE',
                                                       colors = c('SNP' = '#FE817D', 'INS' = '#45BC9C', 'DEL' = '#FFCD6E', 'orange', 'purple'),
                                                       drawColBar = !is.null(input$topIn),
                                                       drawRowBar = !is.null(input$rightIn),
                                                       sampleOrder = samOrder()))
        output$mafShow <- renderUI({fluidRow(plotOutput('mafVaf', width = 80+40*length(mafMutGene), height = 500),
                                             downloadButton('saveVaf', 'Save(.pdf)'),
                                             hr(),
                                             column(width = 1, {br()}),column(width = 11, {plotOutput('mafSomatic')}),
                                             downloadButton('saveMafSomatic', 'Save(.pdf)'),
                                             hr(),
                                             plotOutput('mafGeneWaterFall', width = 1000, height = (600+20*(length(mafMutGene)-20)*(length(mafMutGene)>20))),
                                             downloadButton('saveMafGeneWaterFall', 'Save(.pdf)'))})
      }
      else if(input$mafGeneOrPath == 'Pathway'){
        mafScreen <- subsetMaf(mafTemp, query = 'VARIANT_CLASS %in% mafMutType')
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
    output$creatGmt <- downloadHandler(
      filename = function() {paste(strsplit(input$genesetCsv$name, '[.]')[[1]][1], '.gmt', sep='')},
      content = function(file) {writeGmt(genesetMatrix(), file)})
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
    output$getGOMatrixSelf <- downloadHandler(
      filename = function() { 'GO.csv' },
      content = function(file) {fwrite(goResult()$self[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixUp <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Up.csv', sep='') },
      content = function(file) {fwrite(keggResult()$up[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixDown <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Down.csv', sep='') },
      content = function(file) {fwrite(keggResult()$up[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixAll <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_All.csv', sep='') },
      content = function(file) {fwrite(keggResult()$all[[1]], file, row.names = TRUE)})
    output$getKEGGMatrixSelf <- downloadHandler(
      filename = function() { 'KEGG.csv' },
      content = function(file) {fwrite(keggResult()$self[[1]], file, row.names = TRUE)})
    output$getCorMat <- downloadHandler(
      filename = function() { paste(paste(input$corFactor1[[1]],strsplit(input$clinicalData$name, '[.]')[[1]][1],sep = '_'), '_cor.csv', sep='') },
      content = function(file) {fwrite(corMatResult(), file, row.names = TRUE)})
    output$getCorLs <- downloadHandler(
      filename = function() { paste(paste(input$corFactor1[[1]],strsplit(input$clinicalData$name, '[.]')[[1]][1],sep = '_'), '_corlist.csv', sep='') },
      content = function(file) {fwrite(corLsResult(), file, row.names = TRUE)})
    output$corHeatSave <- downloadHandler(
      filename = function() { paste(paste(input$corFactor1[[1]],strsplit(input$clinicalData$name, '[.]')[[1]][1],sep = '_'), '_corHeatmap.pdf', sep='') },
      content = function(file) {ggsave(file, plot = corHeatMap(), device = 'pdf', 
                                       width = heatMapWidth()/2.5,
                                       height = heatMapHeight()/2.5,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$saveVaf <- downloadHandler(
      filename = function(){paste(strsplit(input$mafVisual$name, '[.]')[[1]][1], '_vaf.pdf', sep='')},
      content = function(file){ggsave(file, plot = plotVaf(mafScreen), device = 'pdf',
                                      width = (80+40*length(mafMutGene))/3,
                                      height = 500/3,
                                      unit = 'mm',
                                      dpi = 600,
                                      limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$saveMafSomatic <- downloadHandler(
      filename = function(){paste(strsplit(input$mafVisual$name, '[.]')[[1]][1], '_geneInteraction.pdf', sep='')},
      content = function(file){
        pdf(file = file)
        somaticInteractions(mafScreen, fontSize = 0.6)
        dev.off()},
      contentType = 'pdf'
    )
    output$saveMafGeneWaterFall <- downloadHandler(
      filename = function(){paste(strsplit(input$mafVisual$name, '[.]')[[1]][1], '_oncoplot.pdf', sep='')},
      content = function(file){
        ggsave(file, plot = oncoplot(mafScreen, top = input$mafTop, 
                                     topBarData = sampleData(),
                                     rightBarData = genesData(),
                                     bgCol = '#EEEEEE',
                                     colors = c('SNP' = '#FE817D', 'INS' = '#45BC9C', 'DEL' = '#FFCD6E', 'orange', 'purple'),
                                     drawColBar = !is.null(input$topIn),
                                     drawRowBar = !is.null(input$rightIn),
                                     sampleOrder = samOrder()), device = 'pdf',
               width = 1000/3,
               height = (600+20*(length(mafMutGene)-20)*(length(mafMutGene)>20))/3,
               units = 'mm',
               dpi = 600,
               limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$surPlotSave <- downloadHandler(
      filename = function() {paste0(strsplit(input$clinicalData$name, '[.]')[[1]][1],'_survival.pdf')},
      content = function(file) {ggsave(file, plot = surResult(), device = 'pdf', 
                                       dpi = 600,
                                       width = ceiling(sqrt(length(isolate(input$surFactor))))*input$surPlotSize/2.5,
                                       height = round(sqrt(length(isolate(input$surFactor))))*input$surPlotSize/2.5,
                                       units = 'mm',
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$coxForestSave <- downloadHandler(
      filename = function() {paste0(strsplit(input$clinicalData$name, '[.]')[[1]][1],'_coxforest.pdf')},
      content = function(file) {ggsave(file, plot = coxForest(), device = 'pdf', 
                                       dpi = 600,
                                       width = input$forestSize/2.5,
                                       height = (100+20*nrow(surResult()))/2.5,
                                       units = 'mm',
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$coxResultSave <- downloadHandler(
      filename = function() {paste0(strsplit(input$clinicalData$name, '[.]')[[1]][1],'_cox.csv')},
      content = function(file) {fwrite(surResult(), file, row.names = TRUE)})
    #GO图下载
    output$goDotUpSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_Up_Dot.pdf', sep='') },
      content = function(file) {ggsave(file, plot = goDotUp(), device = 'pdf', 
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goBarUpSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_Up_Bar.pdf', sep='') },
      content = function(file) {ggsave(file, plot = goBarUp(), device = 'pdf',
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goDotDownSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_Down_Dot.pdf', sep='') },
      content = function(file) {ggsave(file, plot = goDotDown(), device = 'pdf', 
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goBarDownSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_Down_Bar.pdf', sep='') },
      content = function(file) {ggsave(file, plot = goBarDown(), device = 'pdf', 
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goDotAllSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_All_Dot.pdf', sep='') },
      content = function(file) {ggsave(file, plot = goDotAll(), device = 'pdf', 
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goBarAllSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_GO_All_Bar.pdf', sep='') },
      content = function(file) {ggsave(file, plot = goBarAll(), device = 'pdf', 
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goDotSelfSave <- downloadHandler(
      filename = function() { 'GO_Dot.pdf' },
      content = function(file) {ggsave(file, plot = goDotSelf(), device = 'pdf', 
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$goBarSelfSave <- downloadHandler(
      filename = function() { 'GO_Bar.pdf' },
      content = function(file) {ggsave(file, plot = goBarSelf(), device = 'pdf',
                                       width = 800/3,
                                       height = goHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    #kegg图下载
    output$keggDotUpSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Up_Dot.pdf', sep='') },
      content = function(file) {ggsave(file, plot = keggDotUp(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggBarUpSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Up_Bar.pdf', sep='') },
      content = function(file) {ggsave(file, plot = keggBarUp(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggDotDownSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Down_Dot.pdf', sep='') },
      content = function(file) {ggsave(file, plot = keggDotDown(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggBarDownSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_Down_Bar.pdf', sep='') },
      content = function(file) {ggsave(file, plot = keggBarDown(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggDotAllSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_All_Dot.pdf', sep='') },
      content = function(file) {ggsave(file, plot = keggDotAll(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggBarAllSave <- downloadHandler(
      filename = function() { paste(paste(input$deaFactor,strsplit(input$deaData$name, '[.]')[[1]][1],sep = '_'), '_KEGG_All_Bar.pdf', sep='') },
      content = function(file) {ggsave(file, plot = keggBarAll(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggDotSelfSave <- downloadHandler(
      filename = function() { 'KEGG_Dot.pdf' },
      content = function(file) {ggsave(file, plot = keggDotSelf(), device = 'pdf', 
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
    output$keggBarSelfSave <- downloadHandler(
      filename = function() { 'KEGG_Bar.pdf' },
      content = function(file) {ggsave(file, plot = keggBarSelf(), device = 'pdf',
                                       width = 800/3,
                                       height = keggHeight()/3,
                                       units = 'mm',
                                       dpi = 600,
                                       limitsize = FALSE)},
      contentType = 'pdf'
    )
  }
}
