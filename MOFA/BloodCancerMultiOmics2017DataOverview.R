library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

library("BloodCancerMultiOmics2017")
# additional
library("Biobase")
library("SummarizedExperiment")
library("DESeq2")
library("reshape2")
library("ggplot2")
library("dplyr")
library("BiocStyle")

library("ExperimentHub")


data("conctab", "drpar", "lpdAll", "patmeta", "day23rep", "drugs",
     "methData", "validateExp", "dds", "exprTreat", "mutCOM",
     "cytokineViab")

samplesPerData = list(
  drpar = colnames(drpar),
  lpdAll = colnames(lpdAll),
  day23rep = colnames(day23rep),
  methData = colnames(methData),
  patmeta = rownames(patmeta),
  validateExp = unique(validateExp$patientID),
  dds = colData(dds)$PatID,
  exprTreat = unique(pData(exprTreat)$PatientID),
  mutCOM = rownames(mutCOM),
  cytokineViab = unique(cytokineViab$Patient)
)

(samples = sort(unique(unlist(samplesPerData))))
length(samples)

# Number of patients per disease
sort(table(patmeta$Diagnosis), decreasing=TRUE)

# Number of samples from pretreated patients
table(!patmeta$IC50beforeTreatment)

# IGHV status of CLL patients
table(patmeta[patmeta$Diagnosis=="CLL", "IGHV"])

channelNames(drpar)
# show viability data for the first 5 patients and 7 drugs in their lowest conc.
assayData(drpar)[["viaraw.1"]][1:7,1:5]

# number of drugs
nrow(drugs)

# type of information included in the object
colnames(drugs)

head(conctab)

channelNames(day23rep)

# show viability data for 48 h time point for all patients marked as
# replicate 1 and 3 first drugs in all their conc.
drugs2Show = unique(fData(day23rep)$DrugID)[1:3]
assayData(day23rep)[["day2rep1"]][fData(day23rep)$DrugID %in% drugs2Show,]

head(validateExp)

head(cytokineViab)

# there is only one channel with the binary type of data for each gene
channelNames(mutCOM)

# the feature data includes detailed information about mutations in
# TP53 and BRAF genes, as well as clone size of 
#del17p13, KRAS, UMODL1, CREBBP, PRPF8, trisomy12 mutations
colnames(fData(mutCOM))

# show count data for the first 5 patients and 7 genes
assay(dds)[1:7,1:5]

# show the above with patient sample ids
assay(dds)[1:7,1:5] %>% `colnames<-` (colData(dds)$PatID[1:5])

# number of genes and patient samples
nrow(dds); ncol(dds)

# patient samples included in the data set
(p = unique(pData(exprTreat)$PatientID))

# type of metadata included for each gene
colnames(fData(exprTreat))

# show expression level for the first patient and 3 first probes
Biobase::exprs(exprTreat)[1:3, pData(exprTreat)$PatientID==p[1]]

# show the methylation for the first 7 CpGs and the first 5 patient samples
assay(methData)[1:7,1:5]

# type of metadata included for CpGs
colnames(rowData(methData))

# number of patient samples screened with the given platform type
table(colData(methData)$platform)

table(fData(lpdAll)$type)

# show viability data for drug ibrutinib, idelalisib and dasatinib
# (in the mean of the two lowest concentration steps) and
# the first 5 patient samples
Biobase::exprs(lpdAll)[which(
  with(fData(lpdAll),
       name %in% c("ibrutinib", "idelalisib", "dasatinib") &
         subtype=="4:5")), 1:5]

#poor connection, try again
eh = ExperimentHub()
obj = query(eh, "CLLmethylation")
meth = obj[["EH1071"]] # extract the methylation data
