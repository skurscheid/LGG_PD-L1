# survival analysis in GBM patients with RNA-Seq data, stratified by PD-L1 status
require(exactRankTests)
require(survival)
require(coin)

# TCGA data pre-processing ------------------------------------------------
# data downloaded from TCGA on 2016-01-25
setwd("~/Data/TCGA/GBM")
gbm.clinpatient <- read.table("Clinical/Biotab/nationwidechildrens.org_clinical_patient_gbm.txt", header = T, skip = 1, as.is = T, sep = "\t")
gbm.clinpatient <- gbm.clinpatient[-1, ]
save(gbm.clinpatient, file = "gbm.clinpatient.rda")

# had to convert TAB delimited into CSV due to some escape characters
gbm.clinical_radiation <- read.csv("Clinical/Biotab/nationwidechildrens.org_clinical_radiation_gbm.csv", header = T, skip = 1, as.is = T)
gbm.clinical_radiation <- gbm.clinical_radiation[-1, ]

gbm.clinical_drug <- read.table("Clinical/Biotab/nationwidechildrens.org_clinical_drug_gbm.txt", header = T, skip = 1, as.is = T, sep = "\t")
gbm.clinical_drug <- gbm.clinical_drug[-1, ]

# actual analysis ---------------------------------------------------------
# load pre-processed data
load("~/Data/Collaborations/LGG_PDL1/gbm.clinpatient.rda")

i1 <- intersect(gbm.rna_seq_ids, rownames(gbm.clinpatient))

# first let's have a look at the differences in vital status
vitalStats <- rbind(table(gbm.clinpatient[i1,]$vital_status),
                    table(gbm.clinpatient[-which(rownames(gbm.clinpatient) %in% i1), ]$vital_status))
rownames(vitalStats) <- c("GBM RNA-Seq", "GBM pre-RNA-Seq")

vitalStats[, c("Alive","Dead")]
# looks like the two cohorts are different

fisher.test(vitalStats[, c("Alive","Dead")])
# using Fisher's Exact Test confirms it

# now let's look at the survival data
tab1 <- data.frame(cbind(gbm.clinpatient[, c("days_to_death", "vital_status")],
                         "sample_type" = NA,
                         "PDL1" = NA,
                         "PDL1_class" = NA))
tab1$months_to_death <- as(tab1$days_to_death / 30, "numeric")
tab1[i1,]$sample_type <- "GBM RNA-Seq"
tab1[-which(rownames(tab1) %in% i1),]$sample_type <- "GBM pre-RNA-Seq"
tab1$status <- c(0, 1)[match(tab1$vital_status, c("Alive", "Dead"))]
tab1[i1,]$PDL1 <- unlist(gbm.rna_seq.norm.log2[PDL1, i1])
tab1[i1, ][which(tab1[i1, "PDL1"] > 4.908),]$PDL1_class <- "PDL1 high"
tab1[i1, ][which(tab1[i1, "PDL1"] < 4.908),]$PDL1_class <- "PDL1 low"

# stratification based on PD-L1 < mean(PD-L1) > PD-L1
tab2 <- tab1[i1, ][which(tab1[i1, "vital_status"] %in% c("Dead", "Alive")), ]
tab2$PDL1_class <- factor(tab2$PDL1_class, levels = c("PDL1 low", "PDL1 high"))
surv1 <- survfit(Surv(months_to_death, status) ~ PDL1_class, data = tab2)
lsc1 <- cscores(Surv(tab2$months_to_death, tab2$status), int = TRUE)
pt1 <- perm.test(lsc1 ~ PDL1_class, data = tab2)
plot(surv1, lty = 2:3, main = paste("TCGA Survival, GBM RNA-Seq samples, (n = ", nrow(tab2), ", p-value = ", round(pt1$p.value, digits = 4), ")", sep = ""))
legend("topright", c("PDL1 low", "PD-L1 high"), lty = 2:3)

# trying to re-produce the Nduom et al 2015 analysis

q1 <- quantile(tab2$PDL1, c(seq(0, 1, 0.1)))
surv2 <- lapply(q1[seq(2,10,1)], function(x){
  tab2[which(tab2[, "PDL1"] > x), ]$PDL1_class <- "PDL1 high"
  tab2[which(tab2[, "PDL1"] < x), ]$PDL1_class <- "PDL1 low"
  tab2$PDL1_class <- factor(tab2$PDL1_class, levels = c("PDL1 low", "PDL1 high"))
  sv1 <- survfit(Surv(months_to_death, status) ~ PDL1_class, data = tab2)
  pt1 <- perm.test(lsc1 ~ PDL1_class, data = tab2)
  l1 <- list(sv1, pt1)
  return(l1)
})

q1
