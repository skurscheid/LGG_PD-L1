# "tcga_pan_glioma_analysis.R" needs to be run first
require(ggplot2)
require(snowfall)


# local functions ---------------------------------------------------------
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


# global variables --------------------------------------------------------
cpus = 16

# Data pre-processing -----------------------------------------------------
# load Infinium450k annotation data
load("~/Data/TCGA/GBM/Ceccarelli et al. 2016/glioma_annotations.rda")
load("~/Data/References/Annotations/Platforms/gr.annot450k.rda")
setwd("~/Data/Collaborations/LGG_PDL1/")

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

tcga.lgg.450k.path <- "~/Data/TCGA/LGG/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3"
tcga.lgg.450k.files <- list.files("~/Data/TCGA/LGG/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3")
snowfall::sfInit(parallel = T,
       cpus = cpus)
sfExport("tcga.lgg.450k.files")
sfExport("tcga.lgg.450k.path")
sfExport("load_tcga_450k_level3_data")
tcga.lgg.450k <- sfLapply(tcga.lgg.450k.files, function(x){
  t1 <- load_tcga_450k_level3_data(paste(tcga.lgg.450k.path,
                                         x, 
                                         sep = "/"))
  return(t1)
})
tcga.lgg.450k <- do.call("cbind", tcga.lgg.450k)
tcga.lgg.450k <- tcga.lgg.450k[,colnames(tcga.lgg.450k)[grep("Composite", colnames(tcga.lgg.450k), invert=T)]]

#


# prepare data tables -----------------------------------------------------
# load pre-processed data
load("~/Data/Collaborations/LGG_PDL1/tcga_gliomas_rna_seq.rda")
load("~/Data/TCGA/GBM/Infinium450k/tcga.gbm.450k.rda")
load("~/Data/Collaborations/LGG_PDL1/tcga.lgg.450k.rda")

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
LGG.PDL1.methylation_expression <- as.data.frame(cbind(IDH_status = as.character(glioma_annotations[lgg.ids,]$IDH.status),
                                                       grade = as.character(glioma_annotations[lgg.ids,]$Grade),
                                                       histology = as.character(glioma_annotations[lgg.ids,]$Histology),
                                                       chr1p19q_status = as.character(glioma_annotations[lgg.ids,]$X1p.19q.codeletion),
                                                       t(tcga.lgg.450k[names(gr.annot450k[grep("CD274", gr.annot450k$UCSC_RefGene_Name)]), lgg.ids]),
                                                       CD274 = as.numeric(tcga_gliomas_rna_seq[grep("CD274", rownames(tcga_gliomas_rna_seq)), lgg.ids])))
for (i in colnames(LGG.PDL1.methylation_expression[grep("cg|CD", colnames(LGG.PDL1.methylation_expression), invert = F)])){
  LGG.PDL1.methylation_expression[,i] <- as.numeric(as.character(LGG.PDL1.methylation_expression[,i]))
}
LGG.PDL1.methylation_expression$IDH_status <- relevel(LGG.PDL1.methylation_expression$IDH_status, ref = "WT")

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


# PD-L1 visualisation of data using Gviz ----------------------------------------
require(Gviz)
require(annotatr)
gr.annot450k.uscs <- gr.annot450k
seqlevels(gr.annot450k.uscs) <- paste("chr", seqlevels(gr.annot450k.uscs), sep = "")

geneSymbol <- "CD274"
genome <- "hg19"
annots <- c("hg19_cpgs", "hg19_cpg_shores", "hg19_cpg_shelves")
annotations <- annotatr::build_annotations(genome = genome, annotations = annots)
gr.data450k <- gr.annot450k.uscs[grep("CD274", gr.annot450k.uscs$UCSC_RefGene_Name)]
strand(gr.data450k) <- "*"

# preparing data
LGGData <- LGG.PDL1.methylation_expression[, grep("cg", colnames(LGG.PDL1.methylation_expression))]
LGGData <- as.matrix(LGGData)
class(LGGData) <- "numeric"

methExpCorLGG <- sapply(colnames(LGG.PDL1.methylation_expression)[grep("cg", colnames(LGG.PDL1.methylation_expression))], function(x){
  IDHwt <- cor(as.numeric(LGG.PDL1.methylation_expression[, x]), 
           log2(as.numeric(LGG.PDL1.methylation_expression[, "CD274"])+1), 
           method = "spearman")
})

# as per Monika's feedback, only plot top 2 negatively correlated probes
pdf("TCGA_LGG_PDL1_methylation_expression_XY_plots_top2.pdf")
par(mfrow = c(2,1))
methExpXYPlotLGG <- sapply(names(sort(methExpCorLGG)[c(1,2)]), function(x){
  c1 <- cor(as.numeric(LGG.PDL1.methylation_expression[, x]), 
            log2(as.numeric(LGG.PDL1.methylation_expression[, "CD274"])+1), 
            method = "spearman")
  c1 <- round(c1, 2)
  p1 <-plot(as.numeric(LGG.PDL1.methylation_expression[, x]), 
            log2(as.numeric(LGG.PDL1.methylation_expression[, "CD274"])+1),
            xlab = paste("Probe ID ", x, " [beta value]", sep = ""),
            ylab = "PD-L1 [log2(FPKM+1)]",
            xlim = c(0,1),
            main = paste("TCGA LGG\nCorrelation ", c1, " [Spearman]", sep = ""),
            pch = c(16, 17)[as.numeric(glioma_annotations[rownames(LGGData),]$IDH.status)],
            col = c("red", "blue")[as.numeric(glioma_annotations[rownames(LGGData),]$IDH.status)])
  legend("topright", legend = c("IDHmut", "WT"), col = c("red", "blue"), pch = c(16,17))
  abline(lm(log2(as.numeric(LGG.PDL1.methylation_expression[, "CD274"])+1) ~ as.numeric(LGG.PDL1.methylation_expression[, x])), lwd = 2)
  lines(lowess(log2(as.numeric(LGG.PDL1.methylation_expression[, "CD274"])+1) ~ as.numeric(LGG.PDL1.methylation_expression[, x])), col = "green", lwd = 2)
})
dev.off()

GBMData <- GBM.PDL1.methylation_expression[, grep("cg", colnames(GBM.PDL1.methylation_expression))]
GBMData <- as.matrix(GBMData)
class(GBMData) <- "numeric"


# calculating correlations for GBM data
methExpCorGBM <- sapply(colnames(GBM.PDL1.methylation_expression)[grep("cg", colnames(GBM.PDL1.methylation_expression))], function(x){
  IDHwt <- cor(as.numeric(GBM.PDL1.methylation_expression[, x]), 
               log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1), 
               method = "spearman")
})

# XY plots for GBM data
pdf("TCGA_GBM_PDL1_methylation_expression_XY_plots.pdf")
par(mfrow = c(3,2))
methExpXYPlotLGG <- sapply(colnames(GBMData), function(x){
  c1 <- cor(as.numeric(GBM.PDL1.methylation_expression[, x]), 
            log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1),
            method = "spearman")
  c1 <- round(c1, 2)
  p1 <-plot(as.numeric(GBM.PDL1.methylation_expression[, x]), 
            log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1),
            xlab = paste("Probe ID ", x, " [beta value]", sep = ""),
            ylab = "PD-L1 [log2(FPKM+1)]",
            xlim = c(0,1),
            main = paste("Correlation ", c1, " [Spearman]", sep = ""),
            pch = c(16, 17)[as.numeric(glioma_annotations[rownames(GBMData),]$IDH.status)],
            col = c("red", "blue")[as.numeric(glioma_annotations[rownames(GBMData),]$IDH.status)])
  legend("topright", legend = c("IDHmut", "WT"), col = c("red", "blue"), pch = c(16,17))
  abline(lm(log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1) ~ as.numeric(GBM.PDL1.methylation_expression[, x])), lwd = 2)
  lines(lowess(log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1) ~ as.numeric(GBM.PDL1.methylation_expression[, x])), col = "green", lwd = 2)
})
dev.off()

# plot the data for the top 2 CpGs from LGG data
cpgs <- names(sort(methExpCorLGG)[c(1,2)])
pdf("TCGA_GBM_PDL1_methylation_expression_XY_plots_top2.pdf")
par(mfrow = c(2,1))
methExpXYPlotLGG <- sapply(cpgs, function(x){
  c1 <- cor(as.numeric(GBM.PDL1.methylation_expression[, x]), 
            log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1), 
            method = "spearman")
  c1 <- round(c1, 2)
  p1 <-plot(as.numeric(GBM.PDL1.methylation_expression[, x]), 
            log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1),
            xlab = paste("Probe ID ", x, " [beta value]", sep = ""),
            ylab = "PD-L1 [log2(FPKM+1)]",
            xlim = c(0,1),
            main = paste("TCGA GBM\nCorrelation ", c1, " [Spearman]", sep = ""),
            pch = c(16, 17)[as.numeric(glioma_annotations[rownames(GBMData),]$IDH.status)],
            col = c("red", "blue")[as.numeric(glioma_annotations[rownames(GBMData),]$IDH.status)])
  legend("topright", legend = c("IDHmut", "WT"), col = c("red", "blue"), pch = c(16,17))
  abline(lm(log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1) ~ as.numeric(GBM.PDL1.methylation_expression[, x])), lwd = 2)
  lines(lowess(log2(as.numeric(GBM.PDL1.methylation_expression[, "CD274"])+1) ~ as.numeric(GBM.PDL1.methylation_expression[, x])), col = "green", lwd = 2)
})
dev.off()


# preparing Gviz tracks
gtrack <- GenomeAxisTrack()

itrack <- IdeogramTrack(genome = genome, chromosome = "chr9")

biomTrack <- BiomartGeneRegionTrack(genome = genome,
                                    symbol = geneSymbol,
                                    name = "PD-L1\n[Ensembl]",
                                    showId = F)

gr.pdl1 <- GRanges(seqnames = seqnames(biomTrack@range)[1],
                   IRanges(start = start(biomTrack@range)[1] - 5000,
                           end = end(biomTrack@range)[length(end(biomTrack@range))] + 5000),
                   strand = strand(biomTrack@range)[1])

cpgTrack <- AnnotationTrack(subsetByOverlaps(annotations[grep("islands", annotations$type)], gr.pdl1), 
                            col = "lightgreen",
                            fill = "lightgreen",
                            name = "CpG Islands")
displayPars(cpgTrack) <- list(rotation.title = 0, cex.title = 0.6)

pos450kTrack <- AnnotationTrack(gr.annot450k.uscs[grep("CD274", gr.annot450k.uscs$UCSC_RefGene_Name)],
                             col = "black",
                             fill = "black",
                             name = "450k probes")
displayPars(pos450kTrack) <- list(rotation.title = 0, cex.title = 0.6)

LGGdata450kTrackByIDH <- DataTrack(gr.data450k, 
                                data = LGGData, 
                                type = c("a"), 
                                groups = LGG.PDL1.methylation_expression$IDH_status,
                                name = "LGG\n[average beta values/group]")

gr.data450k.boxplots <- gr.data450k
width(gr.data450k.boxplots) <- 100
LGGdata450kTrackByIDHboxplot<- DataTrack(gr.data450k.boxplots, 
                                   data = LGGData, 
                                   type = c("boxplot"), 
                                   groups = LGG.PDL1.methylation_expression$IDH_status,
                                   name = "LGG\n[beta values]")

LGGdataMethExpCor <- DataTrack(gr.data450k, 
                            data = methExpCorLGG,
                            type = c("a"),
                            name = "Spearman\nCorrelation",
                            col = "darkblue")
# plotting
pdf("~/Data/Collaborations/LGG_PDL1/Gviz_plot_PDL1_LGG.pdf", paper = "a4r")
plotTracks(list(gtrack, 
                itrack,
                biomTrack, 
                cpgTrack, 
                pos450kTrack, 
                LGGdata450kTrackByIDH, 
                LGGdataMethExpCor), 
           extend.left = 500, 
           extend.right = 500,
           main = "PD-L1 Methylation - TCGA LGG Data")
dev.off()

pdf("~/Data/Collaborations/LGG_PDL1/Gviz_plot_PDL1_LGG_TSS_detail.pdf", paper = "a4r")
plotTracks(list(gtrack, 
                itrack,
                biomTrack, 
                cpgTrack, 
                pos450kTrack, 
                LGGdata450kTrackByIDH,
                LGGdata450kTrackByIDHboxplot, 
                LGGdataMethExpCor),
           from = start(biomTrack@range)[1] - 800,
           to = end(biomTrack@range)[1] + 500,
           main = "PD-L1 Methylation\nTCGA LGG Data - TSS Detail")
dev.off()


# GBM
# preparing data

GBMdata450kTrackByIDH <- DataTrack(gr.data450k, 
                                data = GBMData, 
                                type = c("a"), 
                                groups = GBM.PDL1.methylation_expression$IDH_status,
                                name = "GBM\n[beta values]")

GBMdata450kTrackByIDHboxplot <- DataTrack(gr.data450k.boxplots, 
                                   data = GBMData, 
                                   type = c("boxplot"), 
                                   groups = GBM.PDL1.methylation_expression$IDH_status,
                                   name = "GBM\n[beta values]")

GBMdataMethExpCor <- DataTrack(gr.data450k, 
                            data = methExpCorGBM,
                            type = c("a"),
                            name = "Spearman\nCorrelation",
                            groups = "CD274",
                            col = "darkblue")
# plotting
pdf("~/Data/Collaborations/LGG_PDL1/Gviz_plot_PDL1_GBM.pdf", paper = "a4r")
plotTracks(list(gtrack, 
                itrack,
                biomTrack, 
                cpgTrack, 
                pos450kTrack, 
                GBMdata450kTrackByIDH, 
                GBMdataMethExpCor), 
           extend.left = 2000, 
           extend.right = 2000,
           main = "PD-L1 Methylation - TCGA GBM Data")
dev.off()

pdf("~/Data/Collaborations/LGG_PDL1/Gviz_plot_PDL1_GBM_TSS_details.pdf", paper = "a4r")
plotTracks(list(gtrack, 
                itrack,
                biomTrack, 
                cpgTrack, 
                pos450kTrack, 
                GBMdata450kTrackByIDH,
                GBMdata450kTrackByIDHboxplot,
                GBMdataMethExpCor), 
           from = start(biomTrack@range)[1] - 800,
           to = end(biomTrack@range)[1] + 500,
           main = "PD-L1 Methylation\nTCGA GBM Data - TSS Details")
dev.off()
levels(LGG.PDL1.methylation_expression$IDH_status)


