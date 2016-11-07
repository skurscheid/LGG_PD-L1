

i1 <- intersect(glioma_annotations$Case, colnames(tcga_gliomas_rna_seq))
lgg_annotations <- glioma_annotations[which(rownames(glioma_annotations) %in% i1 & glioma_annotations$Study == "Glioblastoma multiforme"),]
table(gbm_annotations$IDH.status, useNA = "ifany")
table(gbm_annotations$Supervised.DNA.Methylation.Cluster, useNA = "ifany")
gbm_expression <- tcga_gliomas_rna_seq[, intersect(colnames(tcga_gliomas_rna_seq), rownames(gbm_annotations))]
gbm_expression <- log2(gbm_expression + 1)

plotData <- as.data.frame(glioma_annotations[i1, c("Study", "IDH.status", "IDH.codel.subtype")])
plotData <- cbind(plotData, as.numeric(log2(tcga_gliomas_rna_seq[grep("CD274", rownames(tcga_gliomas_rna_seq)), i1] +1)))
colnames(plotData)[4] <- "PDL1"
levels(plotData$Study) <- c("LGG", "GBM")
plotData$combined <- NA 
plotData[which(plotData$Study == "GBM" & plotData$IDH.codel.subtype == "IDHwt"),]$combined <- "GBM_IDHwt"
plotData[which(plotData$Study == "GBM" & plotData$IDH.codel.subtype == "IDHmut-non-codel"),]$combined <- "GBM_IDHmut"
plotData[which(plotData$Study == "LGG" & plotData$IDH.codel.subtype == "IDHwt"),]$combined <- "LGG_IDHwt"
plotData[which(plotData$Study == "LGG" & plotData$IDH.codel.subtype == "IDHmut-non-codel"),]$combined <- "LGG_IDHmut_nonCoDel"
plotData[which(plotData$Study == "LGG" & plotData$IDH.codel.subtype == "IDHmut-codel"),]$combined <- "LGG_IDHmut_coDel"
plotData$combined <- as.factor(plotData$combined)
table(plotData$combined, useNA = "ifany")
boxplot(plotData$PDL1 ~ plotData$combined)
