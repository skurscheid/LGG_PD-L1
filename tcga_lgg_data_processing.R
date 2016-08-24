#---------- Expression & DNA methylation data on GDU cluster
path <- "/Volumes/gduserv/Data/TCGA/LGG"

#---------- Import list of LGG RNA-Seq samples --------------
lgg.RNASeq.samples <- read.table(paste(path, "METADATA", "UNC__IlluminaHiSeq_RNASeqV2", "unc.edu_LGG.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt", sep = "/"), header = T, as.is = T, sep = "\t")
lgg.RNASeq.samples$sample_id <- unlist(lapply(strsplit(lgg.RNASeq.samples$Comment..TCGA.Barcode., "-"), function(x) paste(x[1:3], collapse = "-")))

# only retain rows with RSEM gene-level normalized data files & filter out duplicated IDs
lgg.RNASeq.samples <- lgg.RNASeq.samples[grep("RSEM_genes_normalized", lgg.RNASeq.samples$Protocol.REF.4),]
lgg.RNASeq.samples <- lgg.RNASeq.samples[!duplicated(lgg.RNASeq.samples$sample_id),]

# get the unique files
lgg.sampleIDs <-lgg.RNASeq.samples$sample_id
lgg.rna_seq_files <- lgg.RNASeq.samples$Derived.Data.File
names(lgg.rna_seq_files) <- lgg.sampleIDs

#---------- TCGA LGG RNA-Seq data processing --------------
path.rna_seq <- paste(path, "/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", sep = "")

# using normalized data
# PD-L1 (CD274) Entrez gene ID: 29126
PDL1 <- grep(29126, rownames(lgg.rna_seq.norm.log2))
for (i in 1:length(lgg.rna_seq_files)){
  if (i == 1){
    temp <- read.table(paste(path.rna_seq, lgg.rna_seq_files[i], sep = ""), header = T, as.is = T, sep = "\t")
  } else {
    temp2 <- read.table(paste(path.rna_seq, lgg.rna_seq_files[i], sep = ""), header = T, as.is = T, sep = "\t")
    temp <- cbind(temp, temp2$normalized_count)
  }
}
lgg.rna_seq.norm <- temp
rm(temp)
rownames(lgg.rna_seq.norm) <- lgg.rna_seq.norm$gene_id
lgg.rna_seq.norm <- lgg.rna_seq.norm[,-1]
colnames(lgg.rna_seq.norm) <- names(lgg.rna_seq_files)
# save(lgg.rna_seq.norm, file = "~/Data/Collaborations/LGG_PDL1/lgg.rna_seq.rda")
# we only retain data for patients with clinical data
lgg.rna_seq.norm <- lgg.rna_seq.norm[,which(colnames(lgg.rna_seq.norm) %in% rownames(clinpatient))]
# log2 transformation of data, adding offset of 1 to avoid -Inf
lgg.rna_seq.norm.log2 <- log2(lgg.rna_seq.norm + 1)
write.csv(lgg.rna_seq.norm.log2, file = "~/Data/Collaborations/LGG_PDL1/lgg.rna_seq.norm.log2.csv")

# save(lgg.rna_seq.norm.log2, file = "~/Data/Collaborations/LGG_PDL1/lgg.rna_seq.norm.log2.rda")
# PD-L1 expression levels:
# summary(unlist(lgg.rna_seq.norm.log2[grep(29126, rownames(lgg.rna_seq.norm.log2)),]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.830   3.674   3.660   4.358   8.285
# human PD-L1 (CD274, EntrezID: 29126) is located on chr 9p24

#---------- TCGA LGG somatic mutations data processing --------------
somaticMutations <- read.table("~/Data/Collaborations/LGG_PDL1/Somatic_Mutations/BI__IlluminaGA_DNASeq_curated/Level_2/broad.mit.edu__IlluminaGA_curated_DNA_sequencing_level2.maf", header = T, sep = "\t", as.is = T)
somaticMutations$sample_id <- unlist(lapply(strsplit(somaticMutations$Tumor_Sample_Barcode, "-"), function(x) paste(x[1:3], collapse = "-")))
somaticMutations[grep("IDH2", somaticMutations$Hugo_Symbol), c("Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification")]
somaticMutations[grep("IDH1", somaticMutations$Hugo_Symbol), c("Hugo_Symbol", "Entrez_Gene_Id", "Variant_Classification")]
length(somaticMutations[grep("IDH2", somaticMutations$Hugo_Symbol), "sample_id"])
length(somaticMutations[grep("IDH1", somaticMutations$Hugo_Symbol), "sample_id"])

#---------- Import list of GBM RNA-Seq samples --------------
gbmPath <- "/home/skurscheid/Data/TCGA/GBM"
gbm.RNASeq.samples <- read.table(paste(gbmPath, "METADATA", "UNC__IlluminaHiSeq_RNASeqV2", "unc.edu_GBM.IlluminaHiSeq_RNASeqV2.1.4.0.sdrf.txt", sep = "/"), header = T, as.is = T, sep = "\t")
gbm.RNASeq.samples$sample_id <- unlist(lapply(strsplit(gbm.RNASeq.samples$Comment..TCGA.Barcode., "-"), function(x) paste(x[1:3], collapse = "-")))

# only retain rows with RSEM gene-level normalized data files & filter out duplicated IDs
gbm.RNASeq.samples <- gbm.RNASeq.samples[grep("RSEM_genes_normalized", gbm.RNASeq.samples$Protocol.REF.4),]
gbm.RNASeq.samples <- gbm.RNASeq.samples[!duplicated(gbm.RNASeq.samples$sample_id),]

# get the unique files
gbm.sampleIDs <-gbm.RNASeq.samples$sample_id
gbm.rna_seq_files <- gbm.RNASeq.samples$Derived.Data.File
names(gbm.rna_seq_files) <- gbm.sampleIDs

#---------- TCGA LGG RNA-Seq data processing --------------
path.rna_seq <- paste(gbmPath, "/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/", sep = "")

# using normalized data
# PD-L1 (CD274) Entrez gene ID: 29126
for (i in 1:length(gbm.rna_seq_files)){
  if (i == 1){
    temp <- read.table(paste(path.rna_seq, gbm.rna_seq_files[i], sep = ""), header = T, as.is = T, sep = "\t")
  } else {
    temp2 <- read.table(paste(path.rna_seq, gbm.rna_seq_files[i], sep = ""), header = T, as.is = T, sep = "\t")
    temp <- cbind(temp, temp2$normalized_count)
  }
}
gbm.rna_seq.norm <- temp
rm(temp)
rownames(gbm.rna_seq.norm) <- gbm.rna_seq.norm$gene_id
gbm.rna_seq.norm <- gbm.rna_seq.norm[,-1]
colnames(gbm.rna_seq.norm) <- names(gbm.rna_seq_files)
# save(lgg.rna_seq.norm, file = "~/Data/Collaborations/LGG_PDL1/lgg.rna_seq.rda")
# we only retain data for patients with clinical data
# log2 transformation of data, adding offset of 1 to avoid -Inf
gbm.rna_seq.norm.log2 <- log2(gbm.rna_seq.norm + 1)
write.csv(gbm.rna_seq.norm.log2, file = "~/Data/Collaborations/LGG_PDL1/gbm.rna_seq.norm.log2.csv")

# save(lgg.rna_seq.norm.log2, file = "~/Data/Collaborations/LGG_PDL1/lgg.rna_seq.norm.log2.rda")
# PD-L1 expression levels:
# summary(unlist(lgg.rna_seq.norm.log2[grep(29126, rownames(lgg.rna_seq.norm.log2)),]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.830   3.674   3.660   4.358   8.285
# human PD-L1 (CD274, EntrezID: 29126) is located on chr 9p24


#---------- load LGG clinical, mutation, MGMT and CIMP status data provided by Pierre --------------
load("~/Data/Collaborations/LGG_PDL1/LGGcimp450K.rda")
load("~/Data/Collaborations/LGG_PDL1/LGGmutation.rda")
load("~/Data/Collaborations/LGG_PDL1/clinpatient.rda")
load("~/Data/Collaborations/LGG_PDL1/MgmtTcgaLGG.rda")
# some data re-formatting
clinpatient$bcr_patient_barcode <- as(clinpatient$bcr_patient_barcode, "character")
write.csv(clinpatient, file = "~/Data/Collaborations/LGG_PDL1/clinpatient.csv")

#---------- load LGG RNA-Seq data from pre-processing step--------------
load("~/Data/Collaborations/LGG_PDL1/lgg.rna_seq.norm.log2.rda")


# load GBM RNA-Seq data from pre-processing step --------------------------
load("~/Data/Collaborations/LGG_PDL1/gbm.rna_seq.norm.log2.rda")

#---------- Infinium 450k annotation data --------------
load("~/Data/Collaborations/LGG_PDL1/gr.annot450k.rda")
# Check if there are any 450k probes associated with PD-L1/CD274
PDL1synonyms <- c("CD274", "B7-H", "B7H1", "PDL1", "PD-L1", "PDCD1L1", "PDCD1LG1")
h <- unlist(sapply(PDL1synonyms, function(x) grep(x,gr.annot450k$UCSC_RefGene_Name)))
gr.annot450k[h]

#-------------- Boxplot of PD-L1 expression levels, comparing LGG Grade II & III --------------
gradeII <- clinpatient[which(clinpatient$grade == 2), "bcr_patient_barcode"]
gradeIII <- clinpatient[which(clinpatient$grade == 3), "bcr_patient_barcode"]
boxplot(unlist(lgg.rna_seq.norm.log2[PDL1, gradeII]),
        unlist(lgg.rna_seq.norm.log2[PDL1, gradeIII]),
        names = c(paste("Grade II\n[N =", length(gradeII), "]", sep = ""),
                  paste("Grade III\n[N =", length(gradeIII), "]", sep = "")))
# T-test to compare the two population means
t.test(lgg.rna_seq.norm.log2[PDL1, gradeII], 
       lgg.rna_seq.norm.log2[PDL1, gradeIII])

#-------------- collate PD-L1 expression levels and survival data for patients with known IDH1 status --------------
IDHmut <- clinpatient[which(clinpatient$mIDHtot == 1), "bcr_patient_barcode"]
IDHwt <- clinpatient[which(clinpatient$mIDHtot == 0), "bcr_patient_barcode"]
IDHmut1p19qcodel <- clinpatient[which(clinpatient$mIDHtot == 1 & clinpatient$co1p19q == "cd"), "bcr_patient_barcode"]
IDHmut1p19qnorm <- clinpatient[which(clinpatient$mIDHtot == 1 & clinpatient$co1p19q == "n"), "bcr_patient_barcode"]
CIMPpos <- clinpatient[which(clinpatient$hCIMP == 1), "bcr_patient_barcode"]
CIMPneg <- clinpatient[which(clinpatient$hCIMP == 0), "bcr_patient_barcode"]

PDL1.IDHstatusKnownPatients  <- c(unlist(lgg.rna_seq.norm.log2[PDL1, IDHmut1p19qcodel]), 
                                  unlist(lgg.rna_seq.norm.log2[PDL1, IDHmut1p19qnorm]),
                                  unlist(lgg.rna_seq.norm.log2[PDL1, IDHwt]))
PDL1.IDHstatusKnownPatients <- data.frame("PDL1exp" = PDL1.IDHstatusKnownPatients, 
                                          "status" = c(rep("IDHmut1p19qcodel", length(IDHmut1p19qcodel)), 
                                                       rep("IDHmut1p19qnorm", length(IDHmut1p19qnorm)),
                                                       rep("IDHwt", length(IDHwt))),
                                          "months_to_death" = clinpatient[c(IDHmut1p19qcodel, IDHmut1p19qnorm, IDHwt), "days_to_death"]/ 30)
PDL1.IDHstatusKnownPatients$death <- 0
PDL1.IDHstatusKnownPatients[!is.na(PDL1.IDHstatusKnownPatients$months_to_death), "death"] <- 1
write.csv(PDL1.IDHstatusKnownPatients, file = "~/Data/Collaborations/LGG_PDL1/PDL1.IDHstatusKnownPatients.csv")

#-------------- PD-L1 expression levels, comparing IDH mut to wt, irrespective of grade --------------
boxplot(PDL1.IDHstatusKnownPatients[IDHmut, "PDL1exp"],
        PDL1.IDHstatusKnownPatients[IDHwt, "PDL1exp"],
        names = c(paste("IDH mut\n[N =", length(IDHmut), "]", sep = ""),
                  paste("IDH wt\n[N =", length(IDHwt), "]", sep = "")))
# T-test to compare the two population means
t.test(PDL1.IDHstatusKnownPatients[IDHmut, "PDL1exp"],
       PDL1.IDHstatusKnownPatients[IDHwt, "PDL1exp"])

#-------------- PD-L1 expression levels, comparing CIMP pos to neg, irrespective of grade --------------
boxplot(unlist(lgg.rna_seq.norm.log2[PDL1, CIMPpos]),
        unlist(lgg.rna_seq.norm.log2[PDL1, CIMPneg]),
        names = c("CIMP+", "CIMP-"))
# T-test to compare the two population means
t.test(lgg.rna_seq.norm.log2[PDL1, CIMPpos], 
       lgg.rna_seq.norm.log2[PDL1, CIMPneg])


#-------------- PD-L1 expression levels, comparing IDHmut/1p19qcodel to IDmut/1p19qn to IDHwt  --------------
boxplot(PDL1.IDHstatusKnownPatients$PDL1exp ~ PDL1.IDHstatusKnownPatients$status,
        names = c("IDH mut,\n1p19q codel", "IDH mut,\n1p19q norm", "IDH wt"))
# performing ANOVA to check if differences between groups are significant
aov1 <- aov(formula = PDL1exp ~ status, data = PDL1.IDHstatusKnownPatients)
# post-hoc analysis shows that they are
TukeyHSD(aov1)
bartlett.test(PDL1.IDHstatusKnownPatients$PDL1exp ~ PDL1.IDHstatusKnownPatients$status)
# but Bartlett test for homogeneity of variances shows that the assumption is violated
# so the ANOVA has to be taken with a grain of salt
# [possibly due to uneven sized groups?]


#-------------- PD-L1 expression levels, comparing IDHmut/1p19qcodel to IDmut/1p19qn to IDHwt, including GBM  --------------
tab1 <- PDL1.IDHstatusKnownPatients[, c("PDL1exp", "status")]
tab2 <- data.frame("PDL1exp" = unlist(gbm.rna_seq.norm.log2[PDL1,]), "status" = rep("GBM", ncol(gbm.rna_seq.norm.log2)))
tab1 <- rbind(tab1, tab2)
tab1$status <- factor(as.character(tab1$status), levels = c("IDHmut1p19qcodel", "IDHmut1p19qnorm", "IDHwt", "GBM"))


boxplot(tab1$PDL1exp ~ tab1$status,
        axes = FALSE,
        frame = FALSE,
        lwd = 2.75,
        ylim = c(0,10))

axis(1, lwd = 2.75,
     at = c(1,2,3,4), 
     labels = c(paste("IDH mut,\n1p19q codel\n[N = ", length(IDHmut1p19qcodel),"]", sep = ""),
                paste("IDH mut,\n1p19q norm\n[N = ", length(IDHmut1p19qnorm),"]", sep = ""),
                paste("IDH wt\n[N = ", length(IDHwt), "]", sep = ""),
                paste("GBM\n[N = ", ncol(gbm.rna_seq.norm.log2), "]", sep = "")),
     padj = 1)
axis(2, lwd = 2.75)
mtext(side = 2, "PD-L1 expression [log2 FPKM]", line = 2.5, cex = 1.25)

# performing ANOVA to check if differences between groups are significant
aov2 <- aov(formula = PDL1exp ~ status, data = tab1)
# post-hoc analysis shows that they are

TukeyHSD(aov2)
bartlett.test(tab1$PDL1exp ~ tab1$status)
# but Bartlett test for homogeneity of variances shows that the assumption is violated
# so the ANOVA has to be taken with a grain of salt
# [possibly due to uneven sized groups? or because values are derived from count data?]

#-------------- collate PD-L1 expression levels and survival data for patients with known CIMP status [more CIMP cases then known IDH1/2 status] --------------
PDL1 <- grep(29126, rownames(lgg.rna_seq.norm.log2))
CIMPpos1p19qcodel <- clinpatient[which(clinpatient$hCIMP == 1 & clinpatient$co1p19q == "cd"), "bcr_patient_barcode"]
CIMPpos1p19qnorm <- clinpatient[which(clinpatient$hCIMP == 1 & clinpatient$co1p19q == "n"), "bcr_patient_barcode"]
CIMPpos <- clinpatient[which(clinpatient$hCIMP == 1), "bcr_patient_barcode"]
CIMPneg <- clinpatient[which(clinpatient$hCIMP == 0), "bcr_patient_barcode"]

PDL1.CIMPstatusKnownPatients  <- c(unlist(lgg.rna_seq.norm.log2[PDL1, CIMPpos1p19qcodel]), 
                                  unlist(lgg.rna_seq.norm.log2[PDL1, CIMPpos1p19qnorm]),
                                  unlist(lgg.rna_seq.norm.log2[PDL1, CIMPneg]))
PDL1.CIMPstatusKnownPatients <- data.frame("PDL1exp" = PDL1.CIMPstatusKnownPatients, 
                                          "status" = c(rep("CIMPpos1p19qcodel", length(CIMPpos1p19qcodel)), 
                                                       rep("CIMPpos1p19qnorm", length(CIMPpos1p19qnorm)),
                                                       rep("CIMPneg", length(CIMPneg))),
                                          "months_to_death" = clinpatient[c(CIMPpos1p19qcodel, CIMPpos1p19qnorm, CIMPneg), "days_to_death"]/ 30)
PDL1.CIMPstatusKnownPatients$death <- 0
PDL1.CIMPstatusKnownPatients[!is.na(PDL1.CIMPstatusKnownPatients$months_to_death), "death"] <- 1
write.csv(PDL1.CIMPstatusKnownPatients, file = "~/Data/Collaborations/LGG_PDL1/PDL1.CIMPstatusKnownPatients.csv")

#-------------- PD-L1 expression levels, comparing CIMPpos/1p19qcodel to CIMPpos/1p19qn to CIMPneg, including GBM  --------------
tab1 <- PDL1.CIMPstatusKnownPatients[, c("PDL1exp", "status")]
tab2 <- data.frame("PDL1exp" = unlist(gbm.rna_seq.norm.log2[PDL1,]), "status" = rep("GBM", ncol(gbm.rna_seq.norm.log2)))
tab1 <- rbind(tab1, tab2)
tab1$status <- factor(as.character(tab1$status), levels = c("CIMPpos1p19qcodel", "CIMPpos1p19qnorm", "CIMPneg", "GBM"))

boxplot(tab1$PDL1exp ~ tab1$status,
        axes = FALSE,
        frame = FALSE,
        lwd = 2.75,
        ylim = c(0,10))

axis(1, lwd = 2.75,
     at = c(1,2,3,4), 
     labels = c(paste("CIMP+,\n1p19q codel\n[N = ", length(CIMPpos1p19qcodel),"]", sep = ""),
               paste("CIMP+,\n1p19q norm\n[N = ", length(CIMPpos1p19qnorm),"]", sep = ""),
               paste("CIMP-\n[N = ", length(CIMPneg), "]", sep = ""),
               paste("GBM\n[N = ", ncol(gbm.rna_seq.norm.log2), "]", sep = "")),
     padj = 1)
axis(2, 
     lwd = 2.75)
mtext(side = 2, "PD-L1 expression [log2 FPKM]", line = 2.5, cex = 1.25)

# performing ANOVA to check if differences between groups are significant
aov2 <- aov(formula = PDL1exp ~ status, data = tab1)
# post-hoc analysis shows that they are

TukeyHSD(aov2)
bartlett.test(tab1$PDL1exp ~ tab1$status)
# but Bartlett test for homogeneity of variances shows that the assumption is violated
# so the ANOVA has to be taken with a grain of salt
# [possibly due to uneven sized groups? or because values are derived from count data?]

