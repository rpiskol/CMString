##' plot pie chart for CMS classification
##'
##' plot pie chart for CMS classification
##' @param df data frame with column that will be used to plot the pie chart [default: none]
##' @param var string specifying which column from the supplied data frame to plot [default: "CMS.glmnet"]
##' @param digits integer specifying to how many digits pie labels should be rounded [default: 1]
##' @param main string specifying the title of the plot [default: ""]
##' @param labTotalN boolean specifying whether the total number of samples in each category should be added in addition to the percentage [default: FALSE]
##' @return vector with percentages of samples in each category
##' @export
##' @author Robert Piskol
##'
pieChart <- function(df, var="CMS.glmnet", digits=1, main="", labTotalN=FALSE){

  if(missing(df)){
    message("data frame has not been specified")
    stop();
  }

  COL_SUB =  getCols("COL_SUB")
  ##############
  # from DF to numbers
  if(!var %in%  colnames(df)){
    message(paste0("variable '", var, "' is not present in the supplied data frame"))
    return();
  }

  tt <- table(df[,var])
  t.perc <- round(100*tt/sum(tt), digits)
  t.perc.text <- paste0(t.perc, "%")
  if(labTotalN){
    t.perc.text <- paste0(t.perc.text, " (n=",tt,")")
  }

  pie(tt, labels = t.perc.text, main = main, col = COL_SUB, clockwise=TRUE)
  legend("topright", names(tt), cex = 0.8,
         fill = COL_SUB)

  return(t.perc)

}

##' plot pie chart for CMS classification
##'
##' plot pie chart for CMS classification
##' @param x vector of values that will be plotted on the x-axis of the contingency table [default: none]
##' @param y vector of values that will be plotted on the y-axis of the contingency table. Needs to be sorted in the same way as x [default: none]
##' @param xlab x-axis label [default: NULL]
##' @param ylab y-axis label [default: NULL]
##' @param main string specifying the title of the plot [default: ""]
##' @return vector with percentages of samples in each category
##' @import dplyr tidyr ggplot2 grid
##' @export
##' @author Robert Piskol
##'
plotContingencyGeneral <- function(x,y,xlab=NULL,ylab=NULL, main=NULL){

  if(missing(x) || missing(y)){
    message("please specify both x and y")
    stop();
  }

  #################
  # create contingency table
  nam <- names(x)
  m<-NULL
  xt <- xtabs(~y + x)

  m <- matrix(xt,nrow=nrow(xt))
  rownames(m) <- rownames(xt)
  colnames(m) <- as.character(colnames(xt))

  ###############
  # get on-/off-diagonal element stats
  lab.share <- intersect(rownames(m), colnames(m))
  n.total <- sum(m[lab.share, lab.share])
  n.diag <- sum(diag(m[lab.share, lab.share]))
  n.offdiag <- n.total - n.diag
  concordance <- n.diag / n.total

  ##################
  # calculate z-scores by which the matrix will be colored
  m.zscore <- t(apply(m,1,zscoreVals))
  m.dat.orig <- data.frame(rnames=rownames(m),m) %>% gather(type,value,-rnames)
  m.dat.zscore <- data.frame(rnames=rownames(m.zscore),m.zscore) %>% gather(type,value,-rnames)

  m.dat <- data.frame(m.dat.zscore, m.dat.orig[,3])
  colnames(m.dat)[c(2,4)] <- c("type","value.orig")

  m.dat$type <- gsub("X","",m.dat$type)

  ##############
  # plot matrix
  p <- ggplot(m.dat, aes(y=rnames,x=type,label=formatC(value.orig,digits=2))) +
    geom_tile(aes(fill=value)) +   geom_text(size=11) +
    theme_bw() + ylab(ylab) + xlab(xlab) +
    scale_fill_gradientn(colours=getCols(cat="HEATMAP_COL3"), guide=FALSE)+
    theme(text = element_text(size=20),
          axis.text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1)#,
    )
  if(!is.null(main)){
    p <- p + ggtitle(main)
  }

  print(p)

  ##############
  # return stats
  stats <- c(n.total,n.diag, n.offdiag, concordance)
  names(stats) <- c("total elements with shared labels", "elements on the diagnal", "off-diagonal elements", "concordance")

  return(stats)
}

##' plot heatmap
##'
##' plot heatmap with gene expression of classifier genes as well as classifications, their probabilities and phenotypic information (mutations)
##' @param x list with classification results (pred, clu.pred, sdat.sig, gclu.f, nam.ord, missing.genes) [default: none]
##' @param pred matrix of class prediction probabilities with samples in rows and classes in the columns.  [default: none]
##' @param clu.pred vector with predicted classes [default: none]
##' @param sdat.sig matrix with gene expression values for the classifier genes with samples in genes in rows and samples in columns [default: none]
##' @param gclu.f hierarchical clustering of classifier genes that specifies row ordering [default: none]
##' @param nam.ord vector of sample IDs that specify column ordering of heatmap [default: none]
##' @param missing.genes vector of genes that will be removed from the heatmap. These are usually genes that were not present in the expression data, but required by the classifier. However any other list can of IDs can be specified [default: NULL]
##' @param cex.axis.labs numeric scaling factor for axis labels (genes, samples) [default: 1]
##' @param mutmat matrix of numeric values that will be translated into colors and attached to the bottom of the heatmap. This can contain mutational information or any other phenotypic information [default: NULL]
##' @param gene.labs boolean specifying whether gene labels should be plotted [default: TRUE]
##' @param samp.labs boolean specifying whether sample labels should be plotted [default: FALSE]
##' @import dplyr tidyr ggplot2 grid
##' @export
##' @author Robert Piskol
##'
figClassifyCMS <- function (x, pred, clu.pred, sdat.sig, gclu.f,
                            nam.ord, missing.genes=NULL,
                            cex.axis.labs=1, mutmat=NULL, gene.labs=TRUE, samp.labs=FALSE)
{

  if(!missing(x)){
    pred.elements <- c("pred", "clu.pred", "sdat.sig", "gclu.f", "nam.ord")

    if(!all(pred.elements %in% names(x))){
      message("The supplied list does not contain all the elements needed to plot the heatmap")
      message(paste("List contains:", names(x)))
      message(paste("Needed:", pred.elements))
      message(paste("Missing:", setdiff(pred.elements, names(x))))

    }

    pred <- x$pred
    clu.pred <- x$clu.pred
    sdat.sig <- x$sdat.sig
    gclu.f <- x$gclu.f
    nam.ord <- x$nam.ord
    if("missing.genes" %in% names(x)){missing.genes <- x$missing.genes}

  } else {

    if(missing(pred)){
      message("missing 'pred'")
      stop()
    }
    if(missing(clu.pred)){
      message("missing 'clu.pred'")
      stop()
    }
    if(missing(sdat.sig)){
      message("missing 'sdat.sig'")
      stop()
    }
    if(missing(gclu.f)){
      message("missing 'gclu.f'")
      stop()
    }
    if(missing(nam.ord)){
      message("missing 'nam.ord'")
      stop()
    }
  }


  ##############
  # set up plot layout and margins
  layout(matrix(c(1:3, 4, 4, 5), 6, 1), heights = c(1.8,20, rep(0.9, 2),3,0.5))
  par(mar = c(0.1, 8, 2.1, 0.5))


  ##################
  # plot cluster predictions
  clu.pred <- clu.pred[nam.ord]
  clu.pred.classes <- na.omit(unique(clu.pred))
  clu.pred <- ifelse(is.na(clu.pred), 5,clu.pred)

  at.breaks <- sapply(seq_along(clu.pred.classes), function(i){max(which(clu.pred == clu.pred.classes[i]))})

  at.marks <- c(at.breaks[1]/2,
                at.breaks[1]+(at.breaks[2]-at.breaks[1])/2,
                at.breaks[2]+(at.breaks[3]-at.breaks[2])/2,
                at.breaks[3]+(at.breaks[4]-at.breaks[3])/2) / length(clu.pred)
  at.lab <- paste0("CMS",clu.pred.classes)

  image(matrix(clu.pred, length(clu.pred), 1), col = getCols(cat="COL_SUB"), axes = FALSE,zlim=c(0,length(getCols(cat="COL_SUB"))))
  axis(side = 3, at = at.marks, labels = at.lab,
       line = -1, tick = FALSE, cex.axis = 1.8)
  box(lwd = 0.5, col = gray(0.5))

  ###################
  # plot expression matrix
  par(mar = c(0.5, 8, 0.1, 0.5))

  sdat4plot <- t(sdat.sig[rev(gclu.f$order), ])
  sdat4plot <- sdat4plot[nam.ord,]

  if(length(missing.genes)>0){
    idx <- na.omit(match(missing.genes, colnames(sdat4plot)))
    sdat4plot <- sdat4plot[,-idx]
  }

  image(sdat4plot, col = getCols(cat="HEATMAP_COL3"),
        breaks = getCols(cat="HEATMAP_BKS3"),
        axes = FALSE)
  gene.table <- NULL

  gene.table <- data.frame(Symbol=colnames(sdat4plot), stringsAsFactors = FALSE)

  if(gene.labs){
    axis(2,at=0:(ncol(sdat4plot)-1)/(ncol(sdat4plot)-1), lab=gene.table$Symbol,las=2,cex.axis=cex.axis.labs)
  }
  if(samp.labs){
    axis(1,at=0:(ncol(sdat.sig)-1)/(ncol(sdat.sig)-1), lab=names(clu.pred),las=2,cex.axis=cex.axis.labs, line=3)
  }
  #axis(2,at=1:ncol(t(sdat.sig)))
  mtext("Genes", side = 2, cex = 1.2, line = 6)
  box(lwd = 0.5, col = gray(0.5))
  par(mar = c(0.1, 8, 0.1, 0.5))

  ##################
  # plot classification probabilities
  barplotmat <- t(pred[nam.ord, ])
  colnames(barplotmat) <- NULL
  barplot(barplotmat, beside = FALSE, axes = FALSE, xaxs = "i",
          yaxs = "i", xlab = "", col = getCols(cat="COL_SUB"), border = NA)
  mtext("Prob.", side = 2, line = 0.2, las = 2, cex = 0.8)
  box(lwd = 0.5, col = gray(0.5))

  ###############
  # plot mutation matrix
  sitecols <- getGgColorHue(3)
  stagecols <- getGgColorHue(4)
  stagecols.discrete <- getCols(cat = "HEATMAP_COL3")[c(1,190)]

  mutmat <- mutmat[nam.ord,]

  if(!is.null(mutmat)){
    image(mutmat, col = c(getCols(cat = "COL_MUT_BK"), getCols(cat = "CBBPALETTE")[c(3,7)], sitecols, stagecols, stagecols.discrete), axes = FALSE, zlim=c(0,11))
    nc <- ncol(mutmat)-1
    box(lwd = 0.5, col = gray(0.5))
    axis(2,(0:nc)/nc,colnames(mutmat), las=1,cex.axis=cex.axis.labs*2)
  }

}
