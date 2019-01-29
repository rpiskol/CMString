## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(CMString)
library(Biobase)

## ------------------------------------------------------------------------
showAvailableExpressionDataSets()

## ------------------------------------------------------------------------
CIT.eset <- getExpressionDataSet(dataset = "CIT")

CIT.eset

head(pData(CIT.eset), n=2)

## ------------------------------------------------------------------------
classification <- classifyCMString(CIT.eset, thresh = .4)

## ------------------------------------------------------------------------
table(classification$clu.pred)

head(classification$pred)

## ------------------------------------------------------------------------
CMS.glmnet <- classification$clu.pred[sampleNames(CIT.eset)]
CIT.eset$CMS.glmnet <- ifelse(is.na(CMS.glmnet), "UNK", 
                                                 paste0("CMS",CMS.glmnet))

## ---- fig.height=5, fig.width=5, fig.align = "center"--------------------
pieChart(pData(CIT.eset), digits = 2)

## ---- fig.height=5, fig.width=5, fig.align = "center"--------------------
plotContingencyGeneral(CIT.eset$CMS.glmnet,
                       CIT.eset$CMS_network,
                       xlab="CMString Classification",
                       ylab="Gold Standard Classification",
                       main = "CIT")

## ---- fig.height=6, fig.width=5, fig.align = "center"--------------------
figClassifyCMS(x=classification)

## ---- fig.height=6, fig.width=5, fig.align = "center"--------------------
figClassifyCMS(x=classification,gene.labs = FALSE, samp.labs = FALSE)

## ---- fig.height=6, fig.width=5, fig.align = "center"--------------------
figClassifyCMS(pred = classification$pred, 
               clu.pred = classification$clu.pred, 
               sdat.sig = classification$sdat.sig,
               gclu.f = rev(classification$gclu.f), 
               nam.ord = rev(classification$nam.ord), 
               missing.genes = classification$missing.genes,
               cex.axis.labs = 0.5, gene.labs = FALSE,samp.labs = FALSE)

## ------------------------------------------------------------------------
Procured.eset <- getExpressionDataSet(dataset = "Procured")

table(Procured.eset$Type)

## ------------------------------------------------------------------------
# classification using Entrez IDs
Procured.exprs.entrez <- exprs(Procured.eset)
rownames(Procured.exprs.entrez) <- fData(Procured.eset)$Entrez
Procured.classification.entrez <- classifyCMString(Procured.exprs.entrez, thresh=.5)

# classification using Symbols
Procured.exprs.symbol <- exprs(Procured.eset)
rownames(Procured.exprs.symbol) <- fData(Procured.eset)$Symbol
Procured.classification.symbol <- classifyCMString(Procured.exprs.symbol, thresh=.5)

xtabs(~Procured.classification.symbol$clu.pred + Procured.classification.entrez$clu.pred)

## ---- fig.height=5, fig.width=4, fig.align = "center", fig.show='hold'----
CMS.glmnet <- Procured.classification.entrez$clu.pred[sampleNames(Procured.eset)]
Procured.eset$CMS.glmnet <- ifelse(is.na(CMS.glmnet), "UNK", paste0("CMS",CMS.glmnet))

# plot pie chart
pieChart(pData(Procured.eset), digits = 0)
# plot heatmap for all samples
figClassifyCMS(x=Procured.classification.entrez,gene.labs = FALSE, samp.labs = FALSE)



## ------------------------------------------------------------------------
Procured.classification.entrez.highConf <- subsetSamplesBySignificance(Procured.classification.entrez, thresh=.5)

## ---- fig.height=5, fig.width=3, fig.align = "center", fig.show='hold'----
sampids.primary <- Procured.eset$Sample[which(Procured.eset$Type == "Primary")]
sampids.met <- Procured.eset$Sample[which(Procured.eset$Type == "Met")]

Procured.classification.entrez.primary <- subsetSamplesByName(Procured.classification.entrez, 
                                                              sampids.primary)
Procured.classification.entrez.met <- subsetSamplesByName(Procured.classification.entrez, 
                                                          sampids.met)

# plot heatmap for primary samples
figClassifyCMS(x=Procured.classification.entrez.primary, gene.labs = FALSE)
# plot heatmap for metastases
figClassifyCMS(x=Procured.classification.entrez.met, gene.labs = FALSE)


## ------------------------------------------------------------------------

# as Entrez IDs
ids.entrez <- lapply(1:4,function(i){getSignatureGenesPerCMS(i, "Entrez" ,"lambda.min")})

# as Symbols
ids.symbol <- lapply(1:4,function(i){getSignatureGenesPerCMS(i, "Symbol" ,"lambda.min")})

sapply(ids.entrez, length)
sapply(ids.symbol, length)

length(unique(unlist(ids.symbol)))
sort(unique(unlist(ids.symbol)))



## ------------------------------------------------------------------------
sessionInfo()

