---
title: "TCGA LGG [IDHmut/MGMGunmeth] Survival data"
author: "Sebastian Kurscheid"
date: "10/25/2016"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r overview}
table(lgg_annotations[which(lgg_annotations$MGMT.promoter.status == "Unmethylated" & lgg_annotations$IDH.status == "WT"),]$Vital.status..1.dead., useNA = "always")
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
