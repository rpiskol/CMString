##' classification of sample gene expression profiles into consensus molecular subtypes (CMS)
##'
##' classification of sample gene expression profiles into consensus molecular subtypes (CMS)
##' @param x ExpressionSet object or gene expression matrix for samples (rows) and genes (columns) [default: none]
##' @param s.stat Value of the penalty parameter lambda at which prediction will be made [default: lambda.min]
##' @param applylog boolean specifying whether log2 transformation should be applied to the data [default: TRUE]
##' @param mediancenter boolean specifying whether expression of single genes should be median centered [default: TRUE]
##' @param zScore boolean specifying whether gene expression should be z-scored [default: TRUE]
##' @param thresh minimum probability necessary for CMS call [default: NULL - the CMS with highest probability will be reported regardless of its probability]
##' @return list with results
##' @export
##' @import glmnet
##' @import Biobase
##' @author Robert Piskol
##'
classifyCMString <- function(x, applylog = TRUE, mediancenter=TRUE, zScore=TRUE,s.stat = "lambda.min", thresh=NULL){

  if(missing(x)){
    message("please provide an Expression set or gene expression matrix for classification")
    stop();
  }
  if (class(x) == "ExpressionSet") {
    sdat <- suppressPackageStartupMessages(Biobase::exprs(x))
  } else if (class(x) == "data.frame") {
    sdat <- as.matrix(x)
  } else  if (class(x) == "matrix") {
    sdat <- x
  } else{
    message("please make sure that the provided object is of class 'ExpressionSet', 'data.frame' or 'matrix'")
    stop()
  }


  ###################
  # determine whether expression matrix rownames are gene symbols or entrez IDs
  idtype <- checkFeatureIds(sdat)

  ###################
  # load classifier
  fit <- getFit(idtype = idtype)
  cvfitlist <- getCvfitList(idtype = idtype)
  alpha <- getAlpha()
  lambda <- getLambda(cvfitlist, alpha)

  ##############
  # get matrix of classifier genes with non-zero coefficients
  classGenesMat <- lapply(coef(fit, s=lambda), as.matrix)
  classGenesMat.bind <- do.call(cbind, classGenesMat)
  classGenesMat.bind <- classGenesMat.bind[-1,]
  classGenesMat.bind.clean <- classGenesMat.bind[-which(rowSums(classGenesMat.bind==0)==ncol(classGenesMat.bind)),]
  #classGenes <- rownames(classGenesMat.bind.clean)

  #############
  # get all genes used in training
  classGenes <- rownames(coef(fit, s=lambda)[[1]])
  classGenes <- classGenes[-1]

  ##############
  # get missing genes and fill them in - needed for glmnet predictor
  missingGenes <- setdiff(classGenes, rownames(sdat))

  sdat <- fillVals(sdat, missingGenes, addJitter = TRUE)
  sdat.genes <- rownames(sdat)

  #############
  # scale data
  sdat <- sdat[classGenes,]
  if(applylog){
    sdat <- logIfNeeded(sdat)
  }
  if(mediancenter){
    message("mediancenter")
    sdat.med = sweep(sdat,1, apply(sdat,1,median,na.rm=T))
    sdat <- sdat.med
  }
  if(zScore){
    message("zScore")
    sdat <- zscoreVals(sdat)
  }

  ################
  # predict classes and their probabilities
  clu.pred <- c(predict(fit, newx = t(sdat), s = lambda, type = "class"))
  prob <- predict(fit, newx = t(sdat), s = lambda, type = "response")

  ##################
  # re-arrange samples based on probabilities
  prob.colnames <- colnames(prob)
  prob.rownames <- rownames(prob)
  prob <- matrix(prob,ncol=ncol(classGenesMat.bind))
  colnames(prob) <- prob.colnames; rownames(prob) <- prob.rownames
  maxr <- apply(prob, 1, max)
  postR <- maxr/(1 - maxr)
  clu.pred <- apply(prob, 1, which.max)
  clu.pred <- sort(clu.pred)
  clu.pred.reord <- NULL
  for (cl in sort(unique(clu.pred))) {
    temp <- names(sort(postR[names(clu.pred[clu.pred == cl])]))
    clu.pred.reord <- c(clu.pred.reord, temp)
  }
  clu.pred <- clu.pred[clu.pred.reord]
  nam.ord <- names(clu.pred)

  sdat.sig <- sdat[intersect(rownames(classGenesMat.bind.clean), sdat.genes), nam.ord]

  ####################
  # cluster genes
  gclu.f <- hclust(as.dist(1 - cor(t(sdat.sig))), method = "ward.D2")

  ####################
  # if necessary, mark samples with low classification probability
  if(!is.null(thresh)){
    is.sig <- apply(prob,1, function(x){any(x>thresh)})
    nam.sig <- names(is.sig[which(is.sig)])
    clu.pred[which(!names(clu.pred) %in% nam.sig)] <- NA
  }

  return(list(sdat.sig = sdat.sig, pred = prob, clu.pred = clu.pred,
              nam.ord = nam.ord, gclu.f = gclu.f, missing.genes=missingGenes))


}
