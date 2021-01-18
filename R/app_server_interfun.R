#'  @import GSVA
#'  @import GSVAdata
#'  @import maftools
#'  @import survival
#'  @import survminer
#'  @import edgeR
#'  @import limma
#'  @import clusterProfiler
#'  @import org.Hs.eg.db
#'  @import DOSE
#'  @import ggplot2
#'  @import ggpubr
#'  @import psych
#'  @import pheatmap
#'  @import data.table
#'  @import GSEABase
#'  @import ggcorrplot
data("interdata", package = 'Wstation')
geneLength <- interdata$gLen
geneSet <- interdata$gSet

interClinical <- interdata$clinical
allCancer <- levels(factor(interClinical$type)) 

#读csv文件
readMatrix <- function(csvPath, rowFlag = TRUE){
    if (is.null(csvPath))
      return(NULL)
    getread <- fread(csvPath$datapath, na.strings = NULL, check.name = FALSE, data.table = FALSE, stringsAsFactors = FALSE)
    if(rowFlag){
      rownames(getread) <- getread[,1]
      getread <- getread[,-1, drop = FALSE]
    }
    return(getread)
}
#写gmt文件
writeGmt <- function(gsCsv, gmtfile){
  gs_list <- tapply(gsCsv[,2],as.factor(gsCsv[,1]),function(x) x)
  i <- 1
  sink(gmtfile)
  for(gs in gs_list){
    cat(names(gs_list)[i])
    i<-i+1
    cat('\tNA\t')
    cat(paste(gs,collapse = '\t'))
    cat('\n')
  }
  sink()
}
#counts数据样本分离
classifyCounts <- function(allData){
  if(all(nchar(colnames(allData)[1]) >= 15, substr(colnames(allData)[1], 1, 4) == 'TCGA')){
      normData <- allData[,as.numeric(substr(colnames(allData), 14, 15)) > 10, drop = FALSE]
      tumorData <- allData[, !(colnames(allData) %in% colnames(normData)), drop = FALSE]
      colnames(normData) <- sapply(colnames(normData), function(x){substr(x, 1, 12)})
      normData <- normData[,!duplicated(colnames(normData)), drop = FALSE]
      colnames(tumorData) <- sapply(colnames(tumorData), function(x){substr(x, 1, 12)})
      tumorData <- tumorData[, !duplicated(colnames(tumorData)), drop = FALSE]
      return(list(nr = normData, tm = tumorData))
  }
  else{return(list(nr = allData, tm = allData))}

}
#TPM标准化
countsToTPM <- function(countsMatrix){
  p <- Progress$new(min = 1, max = 12)
  on.exit(p$close())
  p$set(message = 'TPM Normalizing', detail = 'This may take a while...', value = 1)
  rownames(countsMatrix) <- unlist(lapply(strsplit(rownames(countsMatrix), "[.]"), function(x) x[1]))
  p$set(value = 2)
  rownames(geneLength) <- unlist(lapply(strsplit(rownames(geneLength), "[.]"), function(x) x[1]))
  p$set(value = 3)
  countsMatrix <- countsMatrix[rownames(countsMatrix) %in% rownames(geneLength), , drop = FALSE]
  p$set(value = 4)
  geneLength <- geneLength[match(rownames(countsMatrix), rownames(geneLength)), , drop = FALSE]
  p$set(value = 5)
  tmp <- countsMatrix/geneLength$Gene_length 
  p$set(value = 6)
  tpm_counts <- as.data.frame(t(t(tmp)/colSums(tmp))*1000000)   #标准化
  p$set(value = 7)
  #id转换及去重
  ids <- data.frame(symbol = geneLength$Gene_name, mean = apply(tpm_counts, 1, mean))
  p$set(value = 8)
  ids <- ids[order(ids$mean,decreasing = T), ] #ids$symbol按照ids$median中位数从大到小排列的顺序排序，将对应的行赋值为一个新的ids
  p$set(value = 9)
  ids <- ids[!duplicated(ids$symbol), ] #将symbol这一列取取出重复项，'!'为否，即取出不重复的项，去除重复的gene ，保留每个基因最大表达量结果s
  p$set(value = 10)
  tpm_counts <- tpm_counts[rownames(ids),] #新的ids取出ensembl名，将tpm_counts按照取出的这一列中的每一行组成一个新的
  p$set(value = 11)
  rownames(tpm_counts) <- ids$symbol
  #去除在所有样本中表达为零的基因
  p$set(value = 12)
  tpm_counts <- tpm_counts[rowSums(tpm_counts) >0,] 
  p$set(value = 13)
  return(tpm_counts)
}

#GSVA计算
gsvaCal <- function(mat, setlist, gmtFlag, gmtPath){
  if(gmtFlag){
    geneSet$new <- getGmt(gmtPath)
    setlist[length(setlist)+1]<-'new'
  }
  prolen <- length(setlist)+3
  withProgress(message = 'Calculation in GSVA',
               detail = 'This may take a while...', value = 0,
               expr = {
                   mat <- log(mat+1, 2)
                   incProgress(1/prolen)
                   data <- mat[-(1:nrow(mat)),]
                   incProgress(1/prolen)
                   mat <- as.matrix(mat)
                   incProgress(1/prolen)
                   for (gmt in setlist){
                     gset <- geneSet[[gmt]]
                     vgmt <- gsva(mat, gset,verbose=F,parallel.sz=1)
                     incProgress(1/prolen)
                     vgmt <-  as.data.frame(vgmt)
                     data <- rbind(data,vgmt)  }
                     data <- as.data.frame(t(data))
               })
  return(data)
}

#读maf文件并生成所需的突变信息矩阵
getMaf <- function(variantGeneList, mafPath){
  mafContent <- read.maf(mafPath, isTCGA = TRUE)
  maf_data <- mafContent@data
  mafResult <- sampleList <- mafContent@clinical.data
  for(variantGene in variantGeneList){
    maf_gene <- maf_data[which(maf_data$Hugo_Symbol == variantGene)]
    maf_abstract<-apply(sampleList, MARGIN = 1, function(x){ 
      if(x %in% maf_gene$Tumor_Sample_Barcode){st <- 'Variant' 
                                              ty <- paste(unique(maf_gene[which(maf_gene$Tumor_Sample_Barcode == x)]$VARIANT_CLASS), collapse = ',')}
      else{st <- 'None'
          ty <- 'None'}
      return(c(st, ty))
    })
    maf_abstract <- data.frame(t(as.data.frame(maf_abstract)))
    names(maf_abstract) <- c(paste(variantGene, 'Variant_Status', sep = '_'), paste(variantGene, 'Variant_Type', sep = '_'))
    mafResult <- cbind(mafResult,maf_abstract)
  }
  mafResult <- data.frame(mafResult, row.names = 1)
  return(mafResult)
}

#生存状态二值替换
toSurStatus <- function(dt, status ,event){
  dt[,status] = ifelse(tolower(dt[,status])==tolower(event), 1, 0)
  return(dt)
}

#列名符号换下划线
sbToLine <- function(data){
  names(data) <- gsub("[^[:alnum:]]", '_', names(data))
  return(data)
}

#生存分析
surAnalysis <- function(data, survivaltime, sta, t){
  if(is.numeric(data[,t])){data[,t] <- ifelse(data[,t] > median(data[,t]), 'High', 'Low')}
  data[,survivaltime]<-as.numeric(data[,survivaltime])
  data[,t] <- factor(as.character(data[,t]), order = TRUE)
  ort <- t
  t <- gsub("[^[:alnum:]]", '_', t)
  names(data)[names(data) == ort] <- t
  #s 和f要是全局变量
  s<<-Surv(data[,survivaltime],data[,sta])
  f<<-as.formula(paste0('s','~',t))
  fit<-survfit(f[0:3], data=data)
  names(fit$strata) <- sub(t, ort, names(fit$strata))
  suppressMessages(surplot <- ggsurvplot(fit = fit, data = data,  legend.title=t, palette=c("#F95006", "#33484C"),pval = TRUE,surv.median.line = "hv", legend = c(0.7, 0.95)))
  remove(f, envir = parent.env(environment()))
  remove(s, envir = parent.env(environment()))
  #suppressMessages(ggsave(paste0(t,'.pdf'),plot=print(surplot),path= outpath,width=8, height=8))
  return(surplot)
}

#单因素COX
singleCox <- function(data, survivaltime, sta, t){
  data[,survivaltime]<-as.numeric(data[,survivaltime])
  ort <- t
  t <- gsub("[^[:alnum:]]", '_', t)
  names(data)[names(data) == ort] <- t
  s<-Surv(data[,survivaltime],data[,sta])
  f<-as.formula(paste0('s','~',t))
  result <- coxph(f,data=data)
  snum <- result$means
  result <- summary(result)
  rownames(result$coefficients) <- sub(t, ort, rownames(result$coefficients))
  rownames(result$conf.int) <- sub(t, ort, rownames(result$conf.int))
  result$sample_n <- round(snum*result$n)
  return(result)
}

#多因素COX
multipleCox <- function(data, survivaltime, sta, factable){
  data[,survivaltime]<-as.numeric(data[,survivaltime])
  factable<-as.vector(factable)
  orfactable <- factable
  for(i in 1:length(factable)){
    if(!is.numeric(data[, factable[[i]]])){
      data[, factable[[i]]] <- as.factor(data[,factable[[i]]])
    }
    factable[[i]] <- gsub("[^[:alnum:]]", '_', factable[[i]])
    names(data)[names(data) == orfactable[[i]]] <- factable[[i]]
  }
  s<-Surv(data[,survivaltime],data[,sta])
  multif<-as.formula(paste0('s','~',paste(factable, collapse ="+")))
  result<-coxph(multif,data=data)
  for(b in 1:length(factable)){
    names(result$coefficients)[b] <- sub(factable[[b]] ,orfactable[[b]], names(result$coefficients)[[b]])
  }
  snum <- result$means
  result <- summary(result)
  result$sample_n <- round(snum*result$n)
  return(result)
}

#Cox结果提炼
getCoxTable <- function(data){
  a<- as.data.frame(data$coefficients)
  b <- as.data.frame(data$conf.int)
  a = a[c(2,1,3,4, 5)]
  b = b[-1]
  result = merge(a, b, by = "row.names")
  row.names(result) <- result[,1]
  result <- result[,-1]
  result$sample_number <- data$sample_n
  names(result)[names(result) == 'exp(coef)'] <- 'HR'
  names(result)[names(result) == 'Pr(>|z|)'] <- 'P-Value'
  names(result)[names(result) == 'lower .95'] <- 'low95'
  names(result)[names(result) == 'upper .95'] <- 'up95'
  return(result)
}

#自制森林图
diyForest <- function(data){
  data$factor <- row.names(data)
  data$up95 <- round(data$up95, 2)
  data$low95 <- round(data$low95, 2)
  data$HR <- round(data$HR, 2)
  data$`P-Value` <- round(data$`P-Value`, 3)
  unitX <- (max(data$up95) - min(data$low95))/5
  X_max <- ceiling(max(data$up95)/unitX)*unitX
  X_min <- floor(min(data$low95)/unitX)*unitX
  #标记定位
  data$X_low <- X_min - unitX/2 - unitX -unitX
  data$X_up <- X_min -unitX/2 - unitX
  data$X_Hr <- X_min - unitX/2
  data$X_P <- X_max + unitX/2
  data$Y <- c(1:nrow(data))
  #提取标记坐标和标记
  extractData <- function(xname, pname){
    extract <- data[c(xname,'Y',pname)]
    names(extract)[names(extract) == xname] = 'X'
    names(extract)[names(extract) == pname] = 'Label'
    newRow <- data.frame(extract$X[[1]], max(extract$Y)+1, pname)
    names(newRow)<-names(extract)
    extract <- rbind(extract, newRow)
  }
  #标记表
  anno <- rbind(extractData('X_low', 'low95'), extractData('X_up', 'up95'))
  anno <- rbind(anno, extractData('X_Hr', 'HR'))
  anno <- rbind(anno, extractData('X_P', 'P-Value'))
  forestResult<- ggplot(data, aes(HR, factor))+
    geom_point(aes(size = sample_number, col = factor)) +
    geom_errorbarh(aes(xmax = up95, xmin = low95), height = 0.4)+
    geom_vline(aes(xintercept = 1))+
    scale_x_continuous(limits = c(X_min-3*unitX, X_max+unitX), breaks = seq(X_min, X_max, (X_max-X_min)/5)) +
    xlab(' ')+ylab(" ") + theme_bw() + theme(legend.position ="top") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.background = element_rect(fill ='skyblue', color ='red'),
          axis.text=element_text(size=10, face = "bold"),legend.text=element_text(size=11))+
    geom_text(data=anno, aes(x=X, y=Y, label=Label))+
    scale_y_discrete(expand = c(0, 2))+
    scale_color_discrete(guide = FALSE)
  return(forestResult)
}

#edgeR
deaEdgeR <- function(df, condition, h, control){
  pr <- Progress$new(min=1, max=15)
  on.exit(pr$close())
  condition <-  condition[condition[,1] %in% c(h, control), , drop = FALSE]
  pr$set(message = 'Calculation in EdgeR',detail = 'This may take a while...', value = 1)
  condition <- condition[rownames(condition) %in% colnames(df), , drop = FALSE]
  pr$set(value = 2)
  df <- df[,rownames(condition)]
  pr$set(value = 3)
  condition <- factor(condition[,1])
  pr$set(value = 4)
  condition <- relevel(condition, ref = control)
  pr$set(value = 5)
  genelist<-DGEList(counts=df, group = condition)
  pr$set(value = 6)
  keep<-rowSums(cpm(genelist)>1)>=2
  pr$set(value = 7)
  genelist.filted<-genelist[keep,keep.lib.sizes=FALSE]
  pr$set(value = 8)
  genelist.norm <- calcNormFactors(genelist.filted)
  pr$set(value = 9)
  design <- model.matrix(~condition)
  pr$set(value = 10)
  colnames(design) <- levels(condition)
  pr$set(value = 11)
  genelist.Disp<-estimateDisp(genelist.norm, design,robust=TRUE)
  pr$set(value = 12)
  fit<-glmQLFit(genelist.Disp, design, robust=TRUE)
  pr$set(value = 13)
  res <- glmQLFTest(fit)
  pr$set(value = 14)
  result<-topTags(res,n=Inf)$table
  pr$set(value = 15)
  return(result)
}

#Limma
deaLimma <- function(df, condition, h, control){
  pr <- Progress$new(min=1, max=11)
  on.exit(pr$close())
  df <- t(df)
  pr$set(message = 'Calculation in Limma',detail = 'This may take a while...', value = 1)
  condition <- condition[condition[,1] %in% c(h, control), , drop = FALSE]
  pr$set(value = 2)
  condition <- condition[rownames(condition) %in% colnames(df), , drop = FALSE]
  pr$set(value = 3)
  df <- df[,rownames(condition)]
  pr$set(value = 4)
  group_list <- condition[,1]
  pr$set(value = 5)
  group_list <- relevel(factor(group_list), ref = control)
  pr$set(value = 6)
  design <- model.matrix(~0+group_list)
  pr$set(value = 7)
  colnames(design) <- levels(group_list)
  pr$set(value = 8)
  fit <- lmFit(df,design)
  pr$set(value = 9)
  fit2 <- eBayes(fit)
  pr$set(value = 10)
  result <- topTable(fit2, coef=1, n=Inf)
  pr$set(value = 11)
  return(result)
}

#maf差异分析
mafDiffer <- function(mafIndata, condition, h, control){
  mafEx <- subsetMaf(mafIndata, tsb = row.names(condition[condition[,1] == h, ,drop = FALSE]))
  mafCt <- subsetMaf(mafIndata, tsb = row.names(condition[condition[,1] == control, ,drop = FALSE]))
  result <- mafCompare(mafEx, mafCt, 
                       m1Name = paste(colnames(condition), h, sep = ':'),
                       m2Name = paste(colnames(condition), control, sep = ':'))
  return(result)
}

#火山图
plotVolcano <- function(degData){
  staNum <- length(levels(as.factor(degData$Status)))
  if(staNum == 3){colo <- c('blue','black','red')}
  else{colo <- c('blue','red')}
  pic <- ggplot(data = degData, aes(x = logFC, y=-log10(FDR), color = Status)) +
                   geom_point(alpha = 0.5, size = 1.8) +
                   xlab("log2FC") + ylab("-log10(FDR)") +
                   scale_colour_manual(values = colo) +
                   theme_set(theme_set(theme_bw(base_size=20)))
  return(pic)
}

#id转换器
transferID <- function(emMatrix) {
  rownames(emMatrix) <- unlist(lapply(strsplit(rownames(emMatrix), "[.]"), function(x) x[1]))
  rownames(geneLength) <- unlist(lapply(strsplit(rownames(geneLength), "[.]"), function(x) x[1]))
  emMatrix <- emMatrix[rownames(emMatrix) %in% rownames(geneLength), , drop = FALSE]
  geneLength <- geneLength[match(rownames(emMatrix), rownames(geneLength)), , drop = FALSE]
  emMatrix[,'Symbol'] <- geneLength$Gene_name
  return(emMatrix)
}

#数据分箱
seriesToDiscrete <- function(series, cutoff){
  qu <- quantile(series, c(cutoff, 1-cutoff))
  result <- sapply(series, function(x){
    if(x <= qu[1]){return('Low')}
    else if(x >= qu[2]){return('High')}
    else{return('Median')}
  })
  return(result)
}

#富集分析换id
symToEnt <- function(data){
  if(!is.null(data)){
    DEG <- data$Symbol                                              #取差异分析文件
    degID <- bitr(DEG, fromType = "SYMBOL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )   #id类型转换
    return(degID)
  }
}

#GO富集分析
eGO <- function(deg, pc = 0.1, qc = 0.2, onto){
  eg <- enrichGO(
    gene = deg[,1],
    OrgDb = 'org.Hs.eg.db',
    keyType = "SYMBOL",
    ont = onto,         
    pvalueCutoff = pc,
    pAdjustMethod = "BH",
    qvalueCutoff = qc,
    readable = FALSE, )
  ego <- as.data.frame(eg)
  return(list(ego, eg))
}

#KEGG富集分析
eKegg <- function(deg, pc = 0.1, qc = 0.2){
  eKEG<-enrichKEGG(
    gene = deg[,2],
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = pc,
    pAdjustMethod = "BH",
    qvalueCutoff = qc,
    use_internal_data = FALSE)
  eKEGG <- as.data.frame(eKEG)
  return(list(eKEGG, eKEG))
}
#Bar图
plotBar <- function(data, eway, showNum = 5){
  if(eway == 'GO'){
    if('ONTOLOGY' %in% colnames(data@result)){
      pout <- barplot(data, split="ONTOLOGY",showCategory=showNum)+
            facet_grid(ONTOLOGY~., scale="free")
    }
    else{pout <- barplot(data, showCategory=showNum)}
  }
  else if(eway == 'KEGG'){
    pout <- barplot(data,showCategory=showNum)
  }
  return(pout)
}

#Dot图
plotDot <- function(data, eway, showNum = 5){
  if(eway == 'GO'){
    if('ONTOLOGY' %in% colnames(data@result)){
      pout <- barplot(data, split="ONTOLOGY",showCategory=showNum)+
        facet_grid(ONTOLOGY~., scale="free")
    }
    else{pout <- dotplot(data,showCategory=showNum)}
  }
  else if(eway == 'KEGG'){
    pout <- dotplot(data,showCategory=showNum)
  }
  return(pout)
}

#相关系数
corCal <- function(data, x, y, way, group = 'No group', singleFac = 'All'){
  #基本求相关函数
  corCalculation <- function(dataIn){
    corResult <- corr.test(dataIn[x], dataIn[y], method = way, adjust = 'fdr')
    rN <- NULL
    for(i in y){
      rN <- c(rN, paste(x,i,sep = '-'))
    }
    re2 <- corResult$ci
    row.names(re2) <- rN
    re1 <- as.data.frame(corResult$r)
    return(list(Mat = re1, Ls = re2))
  }
  #不分组
  if(group == 'No group'){
    result <- list()
    corRe <- corCalculation(dataIn = data)
    result$mat <- corRe$Mat
    result$ls <- corRe$Ls
  }
  #分组
  else{
    groupList <- levels(factor(data[[group]]))
    result <- list()
    #多组单对多
    for(singleFac in  x){
      Mat<-data.frame()
      Ls<-data.frame()
      for(t in groupList){
        if(length(data[data[[group]] == t,1])<=3){next}
        oneResult <- corr.test(data[data[[group]] == t, singleFac, drop = FALSE],
                               data[data[[group]] == t, y, drop = FALSE],
                               method = way, adjust = 'fdr')
        row.names(oneResult$r)<-t
        Mat <- rbind(Mat, as.data.frame(oneResult$r))
        rN <- NULL
        for(i in y){
          rN <- c(rN, paste(t,i,sep = '-'))
        }
        oneLs <- oneResult$ci
        row.names(oneLs) <- rN
        Ls <- rbind(Ls, oneLs)
      }
      result[[singleFac]] <- list(mat = Mat, ls = Ls)
    }
    #多组比较
    multiMat <- list()
    multiLs <- data.frame()
    for(k in groupList){
      if(length(data[data[[group]] == k,1])<=3){next}
      oneResult <- corCalculation(data[data[[group]] == k,])
      multiMat[[k]] <- oneResult$Mat
      oneResult$Ls[[group]] <- k
      multiLs <- rbind(multiLs, oneResult$Ls)
    }
    result$All <- list(mat = multiMat, ls = multiLs)
  }
  return(result)
}

#相关系数矩阵筛选
corMatScreen <- function(data, scrData, rcut, pcut){
  scrData <- scrData[(scrData$p > pcut)|(abs(scrData$r) < rcut),,drop = FALSE]
  for(i in row.names(scrData)){
    a <- strsplit(i, split = '-')
    data[a[[1]][1],a[[1]][2]] <- NA 
  }
  return(data)
}

#热图
plotHeat <- function(data){
  if(is.null(data)){return(NULL)}
  pHeat <- ggcorrplot(data)
  return(pHeat)
}

#提取maf文件突变类型
extrcactVariantType <- function(mafObject){
  if(is.null(mafObject)){return(NULL)}
  typeList <- levels(as.factor(mafObject@data$Variant_Type))
  return(typeList)
}

#读gmt文件为dataframe
gmtToDataframe <- function(gmtPath){
  if(is.null(gmtPath)){return(NULL)}
  gmtObj <- GSEABase::getGmt(gmtPath)
  result <- data.frame()
  for (var in names(gmtObj)) {
      varData <- data.frame('Gene_Name' = gmtObj[[var]]@geneIds, 
                            'Geneset_Name' = var,
                            check.names = FALSE,
                            stringsAsFactors = FALSE)
      result <- rbind(result, varData)
  }
  return(result)
}


