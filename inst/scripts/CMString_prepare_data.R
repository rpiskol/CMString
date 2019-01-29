

############################
# get fit and cvfitlist objects
vers <- "v16"; alpha <- 0.2 # like v7 but all folds for final model

fit.entrez.file <- paste0("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/classifiers/Traindata_nanostring-classifier.entrez.lasso.uniqgene.log2.scale.alpha",alpha,".final.",vers,".rda")
fit.symbol.file <- paste0("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/classifiers/Traindata_nanostring-classifier.lasso.uniqgene.log2.scale.alpha",alpha,".final.",vers,".rda")

fit.entrez <- readRDS(fit.entrez.file)
fit.symbol <- readRDS(fit.symbol.file)

cvfitlist.entrez.file <- paste0("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/Traindata_nanostring-classifier.entrez.lasso.uniqgene.log2.scale.cvfitList.",vers,".rda")
cvfitlist.symbol.file <- paste0("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/Traindata_nanostring-classifier.lasso.uniqgene.log2.scale.cvfitList.",vers,".rda")

cvfitlist.entrez <- readRDS(cvfitlist.entrez.file)
cvfitlist.symbol <- readRDS(cvfitlist.symbol.file)


########################
# save fit and cvfitlist objects
save(fit.entrez, file="./data/fit.entrez.RData")
save(fit.symbol, file="./data/fit.symbol.RData")
save(cvfitlist.entrez, file="./data/cvfitlist.entrez.RData")
save(cvfitlist.symbol, file="./data/cvfitlist.symbol.RData")


#######################
# save CRC cell line set
CRC.celllines.eset <- readRDS(file="/gne/home/piskolr/Projects/Collab/deSeM_lab/2018_12_NGS2531_CMS_celllines_invivo/results/NGS2531.eset.entrez.annot.Igis3.0.rds")
pheno.fields <- c("SAMID", "INVENTORY_SAMPLE_NAME", "PRIMARY_TISSUE", "cellline", "CMS.known","CMS.glmnet")
pData(CRC.celllines.eset) <- pData(CRC.celllines.eset)[,pheno.fields]
assayDataElement(CRC.celllines.eset, "voomed.norm") <- NULL
exprs.crc.celllines <- CRC.celllines.eset
save(exprs.crc.celllines, file="./data/exprs.crc.celllines.RData")

#######################
# CIT and MVRM data sets
CIT.eset <- readRDS("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/datasets/CIT.curated.eset.rds")
MVRM.eset <- readRDS("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/datasets/MVRM.curated.eset.rds")

sampleNames(CIT.eset) <- CIT.eset$Sample
sampleNames(MVRM.eset) <- MVRM.eset$Sample

fData(CIT.eset)$Probe.ID <- rownames(CIT.eset)
rownames(CIT.eset) <- fData(CIT.eset)$Gene.EG

fData(MVRM.eset)$Probe.ID <- rownames(MVRM.eset)
rownames(MVRM.eset) <- fData(MVRM.eset)$Gene.EG

assayDataElement(CIT.eset, "se.exprs") <- NULL
assayDataElement(MVRM.eset, "se.exprs") <- NULL


exprs.CIT <- CIT.eset
exprs.MVRM <- MVRM.eset

save(exprs.CIT, file="./data/exprs.CIT.RData")
save(exprs.MVRM, file="./data/exprs.MVRM.RData")

######################
# Procured samples
Procured.eset <- readRDS("/gne/home/piskolr/Projects/Collab/deSauvage_lab/2018_11_CRC_and_met/results/Procured_eset_curated.rds")
exprs.Procured <- Procured.eset

save(exprs.Procured, file="./data/exprs.Procured.RData")

