require(gdata)
require(readr)
require(ggplot2)

# local functions ---------------------------------------------------------

load_gdc_tcga_rna_seq_data <- function(id = NULL, filename = NULL) {
  stopifnot(! is.null(id))
  stopifnot(! is.null(filename))
  fn <- paste(id, "/", filename, sep = "")
  stopifnot(file.exists(fn))
  tab <- readr::read_tsv(fn)
  return(tab)}

load_tcga_agilent_level3_data <- function(filename = NULL) {
  stopifnot(! is.null(filename))
  stopifnot(file.exists(filename))
  tab <- readr::read_tsv(filename, comment = "Composite")
  tcga_id <- unlist(lapply(strsplit(colnames(tab)[2], "-"), function(x) paste(x[1:3], collapse = "-")))
  colnames(tab) <- c("gene_name", tcga_id)
  return(tab)}

iterate_gdc_manifest <- function(gdc_manifest = NULL, type = "normalized_results"){
  stopifnot(length(grep("data.frame", class(gdc_manifest))) == 1)
  sub <- grep(type, gdc_manifest$filename)
  l1 <- lapply(sub, function(x){
    t1 <- load_gdc_tcga_rna_seq_data(gdc_manifest[x,1], gdc_manifest[x, 2])
    return(t1)
  })
  names(l1) <- unlist(lapply(strsplit(unlist(gdc_manifest[sub, "filename"]), "\\."), function(x) x[3]))
  return(l1)
}

# local variables ---------------------------------------------------------
wd <- "~/Data/NCI/GDC/TCGA_Gliomas/R_Analysis"
gdc_manifest_file <- "gdc_manifest.2016-09-25T02_42_46.049303.tsv"
tcga_gbm_metadata_file <- "unc.edu_GBM.IlluminaHiSeq_RNASeqV2.mage-tab.1.4.0/unc.edu_GBM.IlluminaHiSeq_RNASeqV2.1.4.0.sdrf.txt"
tcga_lgg_metadata_file <- "unc.edu_LGG.IlluminaHiSeq_RNASeqV2.mage-tab.1.15.0/unc.edu_LGG.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt"
setwd(wd)

# load the Ceccarelli et al annotation ------------------------------------
glioma_annotations <- read.xls("~/Data/TCGA/GBM/Ceccarelli et al. 2016/mmc2-3.xlsx", 
                               sheet = 1,
                               skip = 1)
glioma_annotations$Case <- as.character(glioma_annotations$Case)
rownames(glioma_annotations) <- glioma_annotations$Case
save(glioma_annotations, file = "glioma_annotations.rda")

# Data pre-processing ---------------------------------------------------
# load the RNA-Seq data
# first GDC manifest
setwd("~/Data/NCI/GDC/TCGA_Gliomas")
gdc_manifest <- readr::read_tsv(gdc_manifest_file)
# data -returns a list
tcga_gliomas_rna_seq <- iterate_gdc_manifest(gdc_manifest)
# sample annotations
tcga_gbm_metadata <- read.table(tcga_gbm_metadata_file, sep = "\t", header = T, as.is = T, stringsAsFactors = F)
tcga_lgg_metadata <- read.table(tcga_lgg_metadata_file, sep = "\t", header = T, as.is = T, stringsAsFactors = F)
tcga_glioma_metadata <- rbind(tcga_gbm_metadata[which(tcga_gbm_metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes_normalized"),],
                              tcga_lgg_metadata[which(tcga_lgg_metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes_normalized"),])

tcga_barcodes <- lapply(names(tcga_gliomas_rna_seq), function(x) {
  print(x)
  w1 <- which(tcga_glioma_metadata$"Extract.Name" == x)
  bc <- tcga_glioma_metadata[w1, "Comment..TCGA.Barcode."]
  bc <- unlist(lapply(strsplit(bc, "-"), function(x) paste(x[1:3], collapse = "-")))
  #df <- data.frame(extract_name = x, tcga_id = bc)
  return(bc)
})

names(tcga_gliomas_rna_seq) <- tcga_barcodes
tcga_gliomas_rna_seq <- do.call("cbind", tcga_gliomas_rna_seq)
rownames(tcga_gliomas_rna_seq) <- tcga_gliomas_rna_seq[,1] 
tcga_gliomas_rna_seq <- tcga_gliomas_rna_seq[, grep("normalized_count", colnames(tcga_gliomas_rna_seq))]
colnames(tcga_gliomas_rna_seq) <- unlist(lapply(strsplit(colnames(tcga_gliomas_rna_seq), "\\."), function(x) x[1]))
#save(tcga_gliomas_rna_seq, file = "R_Analysis/tcga_gliomas_rna_seq.rda")
load("~/Data/Collaborations/LGG_PDL1/tcga_gliomas_rna_seq.rda")

# actual expression analysis ----------------------------------------------
i1 <- intersect(glioma_annotations$Case, colnames(tcga_gliomas_rna_seq))
dim(glioma_annotations[i1,])
dim(tcga_gliomas_rna_seq[,i1])
colnames(glioma_annotations)
table(glioma_annotations[i1,]$IDH.status)
table(glioma_annotations[i1,]$Study)
gbm_annotations <- glioma_annotations[which(rownames(glioma_annotations) %in% i1 & glioma_annotations$Study == "Glioblastoma multiforme"),]
table(gbm_annotations$IDH.status, useNA = "ifany")
table(gbm_annotations$Supervised.DNA.Methylation.Cluster, useNA = "ifany")
gbm_expression <- tcga_gliomas_rna_seq[, intersect(colnames(tcga_gliomas_rna_seq), rownames(gbm_annotations))]
gbm_expression <- log2(gbm_expression + 1)
# use IDH status
pdf("TCGA_GBM_PD-L1.pdf")
boxplot(unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$IDH.status == "WT"),])]),
        unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$IDH.status == "Mutant"),])]),
        names = c("WT [N = 144]", "Mutant [N = 11]"),
        main = "TCGA - GBM [N = 155] - PD-L1", ylab = "log2(TPM + 1)"
        )
dev.off()
t.test(unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$IDH.status == "WT"),])]),
        unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$IDH.status == "Mutant"),])]))

# use G-CIMP status
N_GCIMP_pos <- length(which(gbm_annotations$Original.Subtype == "G-CIMP"))
N_GCIMP_neg <- length(which(gbm_annotations$Original.Subtype != "G-CIMP"))
N_GCIMP_tot <- N_GCIMP_pos + N_GCIMP_neg
pdf("TCGA_GBM_PD-L1_by_GCIMP.pdf")
boxplot(unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$Original.Subtype != "G-CIMP"),])]),
        unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$Original.Subtype == "G-CIMP"),])]),
        names = c(paste("G-CIMP- [N =", N_GCIMP_neg,"]", sep = ""), 
                  paste("G-CIMP+ [N =", N_GCIMP_pos, "]", sep = "")),
        main = paste("TCGA - GBM [N =", N_GCIMP_tot, "] - PD-L1", sep = ""),
        ylab = "log2(TPM + 1)"
)
dev.off()
t.test(unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$Original.Subtype != "G-CIMP"),])]),
       unlist(gbm_expression[grep("CD274", rownames(gbm_expression)), rownames(gbm_annotations[which(gbm_annotations$Original.Subtype == "G-CIMP"),])]))

# load TCGA GBM Agilent Level 3 data --------------------------------------
setwd("~/Data/TCGA/GBM/Agilent")
files <- c(list.files("Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/", full.name = T),
           list.files("Expression-Genes/UNC__AgilentG4502A_07_2/Level_3/", full.name = T))

tcga_gbm_agilent_level_3_data <- lapply(files, function(x) {
  load_tcga_agilent_level3_data(x)})

tcga_gbm_agilent_level_3_data <- do.call(cbind, tcga_gbm_agilent_level_3_data)
rownames(tcga_gbm_agilent_level_3_data) <- tcga_gbm_agilent_level_3_data[,1]
tcga_gbm_agilent_level_3_data <- tcga_gbm_agilent_level_3_data[, grep("TCGA", colnames(tcga_gbm_agilent_level_3_data))]
save(tcga_gbm_agilent_level_3_data, file = "tcga_gbm_agilent_level_3_data.rda")
load("~/Data/Collaborations/LGG_PDL1/tcga_gbm_agilent_level_3_data.rda")

i1 <- intersect(glioma_annotations$Case, colnames(tcga_gbm_agilent_level_3_data))
tcga_gbm_agilent_level_3_data <- tcga_gbm_agilent_level_3_data[, i1]
tcga_gbm_agilent_level_3_data  <- as.matrix(tcga_gbm_agilent_level_3_data )
tcga_gbm_agilent_annotations <- glioma_annotations[i1,]
dim(tcga_gbm_agilent_annotations)
dim(tcga_gbm_agilent_level_3_data)
table(tcga_gbm_agilent_annotations$IDH.status, useNA = "ifany")
table(tcga_gbm_agilent_annotations$IDH.codel.subtype, useNA = "ifany")
table(tcga_gbm_agilent_annotations$Original.Subtype, useNA = "ifany")

N_WT <- length(which(tcga_gbm_agilent_annotations$IDH.status == "WT"))
N_Mutant <- length(which(tcga_gbm_agilent_annotations$IDH.status == "Mutant"))
N_Total <- N_WT + N_Mutant
exp_mt <- as.numeric(tcga_gbm_agilent_level_3_data[grep("CD274", rownames(tcga_gbm_agilent_level_3_data)), rownames(tcga_gbm_agilent_annotations[which(tcga_gbm_agilent_annotations$IDH.status == "Mutant"),])])
exp_wt <- as.numeric(tcga_gbm_agilent_level_3_data[grep("CD274", rownames(tcga_gbm_agilent_level_3_data)), rownames(tcga_gbm_agilent_annotations[which(tcga_gbm_agilent_annotations$IDH.status == "WT"),])])
boxplot(exp_wt,
        exp_mt,
        names = c(paste("WT [N = ", N_WT, "]", sep = ""), 
                  paste("Mutant [N = ", N_Mutant, "]", sep = "")),
        main = paste("TCGA - GBM [N = ", N_Total, "] - PD-L1 (CD274)", sep = ""), 
        ylab = "log2"
)
t.test(exp_wt, exp_mt)


gtable <- data.frame(CD274 = as.numeric(tcga_gbm_agilent_level_3_data[grep("CD274", rownames(tcga_gbm_agilent_level_3_data)), i1]),
                     Original_subtype = tcga_gbm_agilent_annotations[i1, ]$Original.Subtype,
                     IDH_Status = tcga_gbm_agilent_annotations[i1, ]$IDH.status)
gtable$IDH_Status <- as.character(gtable$IDH_Status)
gtable[which(is.na(gtable$IDH_Status)), ]$IDH_Status <- "undetermined"
gtable$CIMP_status <- NA
gtable[which(gtable$Original_subtype == "G-CIMP"),]$CIMP_status <- "G-CIMP"
gtable[which(is.na(gtable$CIMP_status)),]$CIMP_status <- "other"
table(gtable$CIMP_status)
gtable$CIMP_status <- as.factor(gtable$CIMP_status)
gtable$CIMP_status <- relevel(gtable$CIMP_status, ref = "other")
gtable$IDH_Status <- as.factor(gtable$IDH_Status)
gtable$IDH_Status <- relevel(gtable$IDH_Status, ref = "WT")
table(gtable$IDH_Status)

N_GCIMP <- length(which(gtable$CIMP_status == "G-CIMP"))
N_other <- length(which(gtable$CIMP_status == "other"))
p1 <- ggplot(data = gtable, mapping = aes(x = CIMP_status, y = CD274, fill = CIMP_status)) + geom_boxplot() + ggtitle(paste("TCGA - Agilent data - Level 3 - N = [", nrow(gtable), "]", sep = ""))
pdf("TCGA_GBM_PD-L1_by_GCIMP_Agilent.pdf")
p1 + scale_x_discrete("CIMP Status") + scale_fill_discrete("Number of samples", labels = c(paste("Other [", N_other, "]", sep = ""),
                                                                                          paste("Mutant [", N_MUT, "]", sep = "")))

dev.off()
t.test(gtable$CD274 ~ gtable$CIMP_status)

gtable[which(is.na(gtable$IDH_Status)), ] <- "NA"
N_WT <- length(which(gtable$IDH_Status == "WT"))
N_MUT <- length(which(gtable$IDH_Status == "Mutant"))
N_NA <- length(which(gtable$IDH_Status == "undetermined"))
table(gtable$IDH_Status)
p2 <- ggplot(data = gtable, mapping = aes(x = IDH_Status, y = CD274, fill = IDH_Status)) + geom_boxplot() + ggtitle(paste("TCGA - Agilent data - Level 3 - N = [", nrow(gtable), "]", sep = ""))
pdf("TCGA_GBM_PD-L1_by_IDH_Agilent.pdf")
p2 + scale_x_discrete("IDH Status") + scale_fill_discrete("Number of samples", labels = c(paste("WT [", N_WT, "]", sep = ""),
                                                                                          paste("Mutant [", N_MUT, "]", sep = ""),
                                                                                          paste("Und. [", N_NA, "]", sep = "")))
dev.off()
summary(aov(gtable$CD274 ~ gtable$IDH_Status))
dev.off()
