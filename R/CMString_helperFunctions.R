##' list available expression data
##'
##' list available expression data
##' @return vector listing the expression data sets that are available
##' @export
##' @author Robert Piskol
##'
showAvailableExpressionDataSets <- function(){

  datasets <- data(package="CMString")$results
  datasets.expression <- subset(datasets, grepl("^exprs",datasets[,"Item"]))
  datasets.expression.names <- gsub("exprs.","",datasets.expression[,"Item"])

  return(datasets.expression.names)

}

##' load expression data
##'
##' load expression data
##' @param dataset string specifying which data set should be retrieved. Use showAvailableExpressionDataSets() to list available expression data  [default: "CIT"]
##' @return ExpressionSet with the requested data
##' @export
##' @import Biobase
##' @author Robert Piskol
##'
getExpressionDataSet <- function(dataset = showAvailableExpressionDataSets()){

  dataset <- match.arg(dataset)

  message(paste0("loading ",dataset, " dataset"))
  dataset.exprs <- get(data(list=paste0("exprs.",dataset),package="CMString", envir=environment()))

  return(dataset.exprs)

}

##' subset samples from the results of the CMS classification based on a threshold value for the prection
##'
##' subset samples from the results of the CMS classification based on a threshold value for the prection
##' @param x result of the CMS classification [default: none]
##' @param thresh numeric value specifying the threshold that will be used to select samples [default: none]
##' @return subsetted result of the CMS classification
##' @export
##' @author Robert Piskol
##'
subsetSamplesBySignificance <- function(x, thresh){

  if(missing(x)){
    message("result of CMS classification has not been specified")
    stop();
  }
  if(missing(thresh)){
    message("threshold has not been specified")
    stop();
  }


  sdat.sig <- x$sdat.sig; pred = x$pred; clu.pred = x$clu.pred;  nam.ord = x$nam.ord

  if(!is.null(thresh)){
    is.sig <- apply(pred,1, function(x){any(x>thresh)})
    nam.sig <- names(is.sig[which(is.sig)])
    clu.pred[which(!names(clu.pred) %in% nam.sig)] <- NA
  }

  nam.sig <- names(clu.pred[which(!is.na(clu.pred))])

  sdat.sig <- sdat.sig[,colnames(sdat.sig) %in% nam.sig]
  pred <- pred[rownames(pred) %in% nam.sig, ]
  clu.pred <- clu.pred[names(clu.pred)%in% nam.sig]

  idx <- match(nam.sig, nam.ord)
  nam.ord <- nam.ord[idx]

  return(list(sdat.sig = sdat.sig, pred = pred, clu.pred = clu.pred,
              nam.ord = nam.ord, gclu.f = x$gclu.f, missing.genes=x$missing.genes))
}

##' subset samples from the results of the CMS classification using sample names
##'
##' subset samples from the results of the CMS classification using sample names
##' @param x result of the CMS classification [default: none]
##' @param samples vector of sample ids that will be selected from the CMS classification result [default: NULL]
##' @return subsetted result of the CMS classification
##' @export
##' @author Robert Piskol
##'
subsetSamplesByName <- function(x, samples){

  if(missing(x)){
    message("result of CMS classification has not been specified")
    stop();
  }
  if(missing(samples)){
    message("vector of samples that will be selected has not been specified")
    stop();
  }

  sdat.sig <- x$sdat.sig; pred = x$pred; clu.pred = x$clu.pred;  nam.ord = x$nam.ord;
  nam.sig <- samples
  nam.sig <- intersect(samples, names(clu.pred))
  # if(!is.null(thresh)){
  #   is.sig <- apply(pred,1, function(x){any(x>thresh)})
  #   nam.sig <- names(is.sig[which(is.sig)])
  #
  # }
  #

  #clu.pred[which(!names(clu.pred) %in% nam.sig)] <- NA
  #nam.sig <- names(clu.pred[which(!is.na(clu.pred))])

  sdat.sig <- sdat.sig[,colnames(sdat.sig) %in% nam.sig]
  pred <- pred[rownames(pred) %in% nam.sig, ]
  clu.pred <- clu.pred[names(clu.pred)%in% nam.sig]

  idx <- match(names(clu.pred), nam.ord)
  nam.ord <- nam.ord[idx]
  #pred <- pred[nam.ord,]

  return(list(sdat.sig = sdat.sig, pred = pred, clu.pred = clu.pred,
              nam.ord = nam.ord, gclu.f = x$gclu.f, missing.genes=x$missing.genes))
}

##' extract feature (gene) names/ids for a specific CMS subtype from a cvfit file
##'
##' extract feature (gene) names/ids for a specific CMS subtype from a cvfit file
##' @param cms CMS subtype for which features should be reported. This can be a string "CMS1", "CMS2", "CMS3", "CMS4" or numeric value eg. 1,2,3,4 [default: none]
##' @param idtype one of "Entrez" or "Symbol" [default: Entrez]
##' @param s.stat Value of the penalty parameter lambda for which features will be reported [default: lambda.min]
##' @return vector of feature names
##' @export
##' @import glmnet
##' @author Robert Piskol
##'
getSignatureGenesPerCMS <- function(cms, idtype=c("Entrez","Symbol"), s.stat="lambda.min"){

  fit <- getFit(idtype = idtype)
  cvfitlist <- getCvfitList(idtype = idtype)
  alpha <- getAlpha()
  lambda <- getLambda(cvfitlist, alpha, s.stat)

  signatureGenes <- getSignatureGenesPerCMSFromFit(fit, cms, lambda )
  return(signatureGenes)
}

##' extract feature (gene) names/ids for a specific CMS subtype from a cvfit object
##'
##' extract feature (gene) names/ids for a specific CMS subtype from a cvfit object
##' @param cvfit cvfit object generated using the cv.glmnet function [default: none]
##' @param cms CMS subtype for which features should be reported. This can be a string "CMS1", "CMS2", "CMS3", "CMS4" or numeric value eg. 1,2,3,4 [default: none]
##' @param s.stat Value of the penalty parameter lambda for which features will be reported [default: lambda.min]
##' @return vector of feature names
##' @import glmnet
##' @author Robert Piskol
##'
getSignatureGenesPerCMSFromFit <- function(cvfit, cms, s.stat="lambda.min"){

  if(missing(cvfit)){
    message("cvfit object has not been specified")
    stop();
  }
  if(missing(cms)){
    message("CMS subtype has not been specified")
    stop();
  }
  #library(glmnet)
  #cvfit <- readRDS(classfile)

  cms <- gsub("CCS","",cms)
  cms <- gsub("ccs","",cms)
  cms <- gsub("CMS","",cms)
  cms <- as.numeric(cms)

  #classGenes <- rownames(coef(cvfit, s=s.stat)[[cms]])
  #classGenes <- classGenes[-1]
  signatureGenes <- coef(cvfit, s=s.stat)[[cms]]
  signatureGenes <- names(signatureGenes[signatureGenes[,1]!=0,1])
  signatureGenes <- signatureGenes[-1]
  return(signatureGenes)
}

##' get colors
##'
##' get colors
##' @param cat the category for which to return colors
##' @return vector of colors or color breaks depending on what category was requested
##' @export
##' @author Robert Piskol
##'
getCols <- function(cat=c("COL_MUT", "COL_MUT_BK", "COL_SUB", "HEATMAP_COL3", "HEATMAP_BKS3", "CBBPALETTE")){

  if(missing(cat)){
    message("no color category specified")
    message(paste("please specify one of:",paste(eval(formals("getCols")$cat), collapse=", ")))
    stop();
  }

  #############
  # check args
  cat <- match.arg(cat)

  #############
  # CMS colors
  COL_SUB <- c("#E08C27","#0B5E9B","#C46193","#138F63", gray(.5))

  ###########
  # mutation colors
  COL_MUT <- colors()[99]
  COL_MUT_BK <- "white"
  COL_MUT_NA <- gray(0.5)

  ############
  # heatmap colors
  COL_GE_HIGH <- "darkorange"
  COL_GE_LOW <- "#1483E4"
  COL_GE_ZERO <- "white"

  HEATMAP_COL3 <- c(colorRampPalette(c(COL_GE_LOW, "white"))(120)[1:99],
                    "white",
                    colorRampPalette(c("white", COL_GE_HIGH))(120)[22:120])
  HEATMAP_BKS3 <- c(-20, -99:(-1)*0.02, 1:99*0.02, 20)

  ###############
  # CBBPALETTE
  CBBPALETTE <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


  return(get(cat))


}

##' add missing genes to expression matrix. Gene expression for these genes in a particular sample will be filled with the average of gene expression across all genes in that sample
##'
##' add missing genes to expression matrix. Gene expression for these genes in a particular sample will be filled with the average of gene expression across all genes in that sample
##' @param sdat gene expression matrix to which missing gene expression will be added (default: none)
##' @param missingGenes vector of gene IDs/symbols that will be added
##' @param addJitter boolean value specifying whether a small amount of jitter should be added to the expression of the supplemented genes (default: FALSE)
##' @param verbose specify whether messages should be printed (default: FALSE)
##' @return matrix with added values
##' @author Robert Piskol
##'
fillVals <- function(sdat, missingGenes, addJitter=FALSE, verbose = FALSE){
  if(missing(sdat)){
    message("gene expression matrix has not been specified")
    stop();
  }
  if(missing(missingGenes)){
    message("vector of missing genes has not been specified")
    stop();
  }

  if(length(missingGenes) == 0){
    return(sdat)
  }

  if(verbose){
    message(paste0("filling in expression values for ", length(missingGenes), " missing genes"))
    message(paste(missingGenes, collapse=","))
  }

  m <- rowMeans(sdat, na.rm = T)
  mm <- mean(m, na.rm = T)
  tmp <- array(mm, dim = c(length(missingGenes), ncol(sdat)))

  if(addJitter){
    set.seed(0815)
    tmp <- jitter(tmp, 0.01)
  }
  tmp <- as.data.frame(tmp)

  rownames(tmp) <- missingGenes
  colnames(tmp) <- colnames(sdat)
  sdat.out <- as.matrix(rbind(sdat, tmp))

  return(sdat.out)
}

##' generate random string
##'
##' generate random string
##' @param length integer specifying the length of the string to be generated (default: 10)
##' @return string
##' @author Robert Piskol
##'
randomString <- function(length = 10){
  paste(sample(c(0:9, letters, LETTERS), length, replace = TRUE),
        collapse = "")
}

##' apply log2 scale to matrix if needed
##'
##' apply log2 scale to matrix if needed
##' @param ex a matrix of values that will be log2 transformed (default: none)
##' @param verbose specify whether messages should be printed (default: FALSE)
##' @return matrix
##' @author Robert Piskol
##'
logIfNeeded <- function(ex, verbose = FALSE){

  if(missing(ex)){
    message("matrix 'ex' has not been specified")
    stop();
  }

  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    if(verbose){message("applying log transform...")}
    ex[which(ex <= 0)] <- NaN
    sdat <- log2(ex)
    return(sdat)
  }
  else{
    if(verbose){message("log transform not necessary")}
    return(ex)
  }
}

##' z-score expression matrix
##'
##' z-score expression matrix
##' @param x a matrix of values that will be z-scored (default: none)
##' @param center integer value to which to center the matrix to (default: NA - matrix will be centered to mean of all values)
##' @return matrix
##' @author Robert Piskol
##'
zscoreVals <- function (x, center = NA){

  if(missing(x)){
    message("matrix 'x' has not been specified")
    stop();
  }

  if (is.na(center))
    center <- mean(x, na.rm=TRUE)
  sd.x <- sd(x, na.rm = TRUE)
  z <- (x - center)/sd.x
  return(z)
}



##' reproduce ggplot's color hues
##'
##' reproduce ggplot's color hues
##' @param n the number of colors to generate
##' @return vector of colors
##' @export
##' @author Robert Piskol
##'
getGgColorHue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

##' get alpha parameter
##'
##' get alpha parameter
##' @return alpha parameter that was used for classifier training
##' @author Robert Piskol
##'
getAlpha <- function(){
  return(0.2)
}

##' get glmnet object with the trained classifier
##'
##' get glmnet object with the trained classifier
##' @param idtype one of "Entrez" or "Symbol" (default: Entrez)
##' @return glmnet fit object

##' @author Robert Piskol
##'
getFit <- function(idtype=c("Entrez", "Symbol")){
  idtype <- match.arg(idtype)

  if(idtype == "Entrez"){
    data(fit.entrez, package="CMString", envir=environment())
    return(fit.entrez)
  }else if(idtype == "Symbol"){
    data(fit.symbol, package="CMString", envir=environment())
    return(fit.symbol)
  }

}

##' get list of objects generated using the cv.glmnet function for multiple values of alpha
##'
##' get list of objects generated using the cv.glmnet function for multiple values of alpha
##' @param idtype one of "Entrez" or "Symbol" (default: "Entrez")
##' @return list of cv.glmnet objects at different alpha levels
##' @author Robert Piskol
##'
getCvfitList <- function(idtype=c("Entrez", "Symbol")){
  idtype <- match.arg(idtype)

  if(idtype == "Entrez"){
    data(cvfitlist.entrez, package="CMString", envir=environment())
    return(cvfitlist.entrez)
  }else if(idtype == "Symbol"){
    data(cvfitlist.symbol, package="CMString", envir=environment())
    return(cvfitlist.symbol)
  }

}

##' get value of the penalty parameter lambda from crossvalidation for a specific alpha value
##'
##' get list of objects generated using the cv.glmnet function for multiple values of alpha
##' @param cvfitList list of cv.glmnet objects at different alpha levels
##' @param alpha value of the alpha parameter for which a specific lambda should be retrieved
##' @param s.stat either "lambda.min" or "lambda.1se" (default "lambda.min")
##' @return list of cv.glmnet objects at different alpha levels
##' @author Robert Piskol
##'
getLambda <- function(cvfitList, alpha, s.stat = c("lambda.min", "lambda.1se")){

  s.stat <- match.arg(s.stat)
  lambda <- cvfitList[[which(names(cvfitList) == alpha)]][[s.stat]]
  return(lambda)

}

##' detect whether gene expression matrix is labeled with Entrez gene ids or Symbols
##'
##' detect whether gene expression matrix is labeled with Entrez gene ids or Symbols
##' @param sdat gene expression matrix that will be used for CMS assignment (default: none)
##' @param minoverlapFrac minimum overlap fraction required between rownames of expression matrix and either Entrez IDs or gene Symbols used for classification
##' @return String specifying whether Entrez IDs or Symbols should be used for classification
##' @author Robert Piskol
##'
checkFeatureIds <- function(sdat, minoverlapFrac = .25){


  ids.current <- rownames(sdat)
  #message(ids.current)

  alpha <- getAlpha()

  ###############
  # get entrez IDs used in the classifier
  fit.entrez <- getFit(idtype = "Entrez")
  cvfitlist.entrez <- getCvfitList(idtype = "Entrez")
  lambda.entrez <- getLambda(cvfitlist.entrez, alpha)
  classGenes.entrez <- rownames(coef(fit.entrez, s=lambda.entrez)[[1]])
  classGenes.entrez <- classGenes.entrez[-1]

  ##############
  # get symbols used in the classifier
  fit.symbol <- getFit(idtype="Symbol")
  cvfitlist.symbol <- getCvfitList(idtype = "Symbol")
  lambda.symbol <- getLambda(cvfitlist.symbol, alpha)
  classGenes.symbol <- rownames(coef(fit.symbol, s=lambda.symbol)[[1]])
  classGenes.symbol <- classGenes.symbol[-1]

  #############
  # get overlaps
  ids.shared.entrez <- intersect(ids.current, classGenes.entrez)
  ids.shared.symbol <- intersect(ids.current, classGenes.symbol)

  frac.shared.entrez <- length(ids.shared.entrez)/length(classGenes.entrez)
  frac.shared.symbol <- length(ids.shared.symbol)/length(classGenes.symbol)

  if(frac.shared.entrez > frac.shared.symbol & frac.shared.entrez > minoverlapFrac){
    message("assuming Entrez IDs as gene names")
    return("Entrez")
  }else if(frac.shared.symbol > frac.shared.entrez & frac.shared.symbol > minoverlapFrac){
    message("assuming Symbols as gene names")
    return("Symbol")
  }else{
    message("the supplied gene expression matrix overlaps to less than 25% with both EntrezGene ids or Gene Symbols")
    stop()
  }

}



