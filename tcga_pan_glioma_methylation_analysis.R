# "tcga_pan_glioma_analysis.R" needs to be run first
require(ggplot2)
# local functions
load_tcga_450k_level3_data <- function(filename = NULL){
  stopifnot(!is.null(filename))
  tab <- readr::read_tsv(filename, skip = 1)
  id <- unlist(lapply(strsplit(filename, "\\."), function(x) x[6]))
  id <- unlist(lapply(strsplit(id, "-"), function(x) paste(x[1:3], collapse = "-")))
  colnames(tab)[2] <- id
  tab <- as.data.frame(tab[,c(1:2)])
  rownames(tab) <- tab[,1]
  return(tab)
}

# load Infinium450k annotation data
load("~/Data/References/Annotations/Platforms/gr.annot450k.rda")

setwd("~/mount/gduserv/Data/TCGA/GBM/Infinium450k/")

gdc_manifest_file <- "gdc_manifest.2016-10-15T02_24_40.986744.tsv"
gdc_manifest <- readr::read_tsv(gdc_manifest_file)

tcga.gbm.450k <- apply(gdc_manifest, 1, function(x) {
  t1 <- load_tcga_450k_level3_data(paste(x["id"], 
                                         x["filename"],
                                         sep = "/"))
  return(t1)
})

tcga.gbm.450k <- do.call("cbind", tcga.gbm.450k)
tcga.gbm.450k <- tcga.gbm.450k[,colnames(tcga.gbm.450k)[grep("Composite", colnames(tcga.gbm.450k), invert=T)]]
save(tcga.gbm.450k, file = "tcga.gbm.450k.rda")

tcga.lgg.450k.path <- "~/mount/gduserv/Data/TCGA/LGG/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3"
tcga.lgg.450k.files <- list.files("~/mount/gduserv/Data/TCGA/LGG/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3")
tcga.lgg.450k <- lapply(tcga.lgg.450k.files, function(x){
  t1 <- load_tcga_450k_level3_data(paste(tcga.lgg.450k.path,
                                         x, 
                                         sep = "/"))
  return(t1)
})
tcga.lgg.450k <- do.call("cbind", tcga.lgg.450k)
tcga.lgg.450k <- tcga.lgg.450k[,colnames(tcga.lgg.450k)[grep("Composite", colnames(tcga.lgg.450k), invert=T)]]

#
load("/Users/u1001407/Data/Collaborations/LGG_PDL1/tcga_gliomas_rna_seq.rda")

ids <- intersect(rownames(glioma_annotations), colnames(tcga.gbm.450k))

# GBM data
# get CD274 probes
idss <- intersect(ids, colnames(tcga_gliomas_rna_seq))
GBM.PDL1.methylation_expression <- as.data.frame(cbind(IDH_status = as.character(glioma_annotations[idss, "IDH.status"]),
                                                       t(tcga.gbm.450k[names(gr.annot450k[grep("CD274", gr.annot450k$UCSC_RefGene_Name)]), idss]),
                                                       CD274 = as.numeric(tcga_gliomas_rna_seq[grep("CD274", rownames(tcga_gliomas_rna_seq)), idss])))
for (i in colnames(GBM.PDL1.methylation_expression[grep("IDH", colnames(GBM.PDL1.methylation_expression), invert = T)])){
  GBM.PDL1.methylation_expression[,i] <- as.numeric(as.character(GBM.PDL1.methylation_expression[,i]))
}
GBM.PDL1.methylation_expression$IDH_status <- as.factor(GBM.PDL1.methylation_expression$IDH_status)
GBM.PDL1.methylation_expression$IDH_status <- relevel(GBM.PDL1.methylation_expression$IDH_status, ref = "WT")

GBM.PDL1.methylation.plots <- lapply(colnames(GBM.PDL1.methylation_expression[grep("cg", colnames(GBM.PDL1.methylation_expression))]), function(z){
  dat <- GBM.PDL1.methylation_expression[!is.na(GBM.PDL1.methylation_expression$IDH_status), c("IDH_status", z, "CD274")]
  colnames(dat)[2] <- "beta_value"
  cor1 <- round(cor(dat[,2], dat[,3]), 2)
  N_WT <- length(which(dat$IDH_status == "WT"))
  N_MUT <- length(which(dat$IDH_status == "Mutant"))
  ttest <- t.test(beta_value ~ IDH_status, data = dat)
  p1 <- ggplot2::ggplot(data = dat, mapping = aes(x = IDH_status, y = beta_value, fill = IDH_status)) + geom_boxplot() + ggtitle(paste("TCGA GBM PD-L1 CpG ", z, "\nCorrelation with PD-L1 expression ", cor1, sep = ""))
  p1 <- p1 + scale_x_discrete(paste("IDH Status [pval=", round(ttest$p.value, 4), "]", sep = "")) + scale_fill_discrete("Number of samples", labels = c(paste("WT [", N_WT, "]", sep = ""),
                                                                                                                                                        paste("Mutant [", N_MUT, "]", sep = "")))
  ttest <- t.test(beta_value ~ IDH_status, data = dat)
  pdf(paste("~/Data/Collaborations/LGG_PDL1/GBM/TCGA_GBM_PD-L1_methylation_", z, ".pdf", sep = ""))
  print(p1)
  dev.off()
  return(list(plot = p1,
              ttest = ttest))
})
names(GBM.PDL1.methylation.plots) <- colnames(GBM.PDL1.methylation_expression[grep("cg", colnames(GBM.PDL1.methylation_expression))])

# LGG data

lgg.ids <- intersect(colnames(tcga.lgg.450k), rownames(glioma_annotations))
lgg.ids <- intersect(lgg.ids, colnames(tcga_gliomas_rna_seq))
LGG.PDL1.methylation_expression <- as.data.frame(cbind(IDH_status = as.character(glioma_annotations[lgg.ids, "IDH.status"]),
                                                       chr1p19q_status = as.character(glioma_annotations[lgg.ids, "X1p.19q.codeletion"]),
                                                       t(tcga.lgg.450k[names(gr.annot450k[grep("CD274", gr.annot450k$UCSC_RefGene_Name)]), lgg.ids]),
                                                       CD274 = as.numeric(tcga_gliomas_rna_seq[grep("CD274", rownames(tcga_gliomas_rna_seq)), lgg.ids])))
for (i in colnames(LGG.PDL1.methylation_expression[grep("status", colnames(LGG.PDL1.methylation_expression), invert = T)])){
  LGG.PDL1.methylation_expression[,i] <- as.numeric(as.character(LGG.PDL1.methylation_expression[,i]))
}
table(LGG.PDL1.methylation_expression$IDH_status, LGG.PDL1.methylation_expression$chr1p19q_status)
lgg.idh_mut.codel <- which(LGG.PDL1.methylation_expression$IDH_status == "Mutant" & LGG.PDL1.methylation_expression$chr1p19q_status == "codel")
lgg.idh_mut.non_codel <- which(LGG.PDL1.methylation_expression$IDH_status == "Mutant" & LGG.PDL1.methylation_expression$chr1p19q_status == "non-codel")
lgg.idh_wt <- which(LGG.PDL1.methylation_expression$IDH_status == "WT")
LGG.PDL1.methylation_expression$molecular_subtype <- NA
LGG.PDL1.methylation_expression[lgg.idh_mut.codel,]$molecular_subtype <- "IDHmut_1p19qcodel"
LGG.PDL1.methylation_expression[lgg.idh_mut.non_codel,]$molecular_subtype <- "IDHmut_1p19qnorm"
LGG.PDL1.methylation_expression[lgg.idh_wt,]$molecular_subtype <- "IDHwt"
LGG.PDL1.methylation_expression$molecular_subtype <- as.factor(LGG.PDL1.methylation_expression$molecular_subtype)
LGG.PDL1.methylation_expression$molecular_subtype <- relevel(LGG.PDL1.methylation_expression$molecular_subtype, ref = "IDHwt")

LGG.PDL1.methylation.plots <- lapply(colnames(LGG.PDL1.methylation_expression[grep("cg", colnames(LGG.PDL1.methylation_expression))]), function(z){
  dat <- LGG.PDL1.methylation_expression[!is.na(LGG.PDL1.methylation_expression$molecular_subtype), c("molecular_subtype", z, "CD274")]
  mst <- table(dat$molecular_subtype)
  colnames(dat)[2] <- "beta_value"
  aov1 <- aov(formula = beta_value ~ molecular_subtype, data = dat)
  sum_aov1 <- summary(aov1)
  cor1 <- round(cor(dat[, 2], dat[,3]), 2)
  p1 <- ggplot2::ggplot(data = dat, mapping = aes(x = molecular_subtype, y = beta_value, fill = molecular_subtype)) + geom_boxplot() + ggtitle(paste("TCGA LGG PD-L1 CpG ",z, "\nCorrelation with PD-L1 expression ", cor1, sep = ""))
  p1 <- p1 + scale_x_discrete(paste("Molecular subtype [pval=", round(sum_aov1[[1]]$`Pr(>F)`, 4), " (AOV)]", sep = "")) + scale_fill_discrete("Number of samples", labels = sapply(1:length(mst), function(x) paste(names(mst[x]), " [", mst[x], "]", sep = "")))
  pdf(paste("~/Data/Collaborations/LGG_PDL1/LGG/TCGA_LGG_PD-L1_methylation_", z, ".pdf", sep = ""))
  print(p1)
  dev.off()
  return(list(plot = p1,
              anova = aov1))
})
names(LGG.PDL1.methylation.plots) <- colnames(LGG.PDL1.methylation_expression[grep("cg", colnames(LGG.PDL1.methylation_expression))])
