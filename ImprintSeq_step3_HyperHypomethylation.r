
rm(list = ls(all = T))

#####################################################################


#########
# Setup #
#########

Calibration_list = c("L4o", "L4n")
SD = 3   # ***
eps = 2.2204e-16   # ***

#####################################################################

#####################################################################


#####################
# Loading libraries #
#####################

library("ggfortify")
library("gplots")
library("ggplot2")
library("grid")
library("pheatmap")
library("RColorBrewer")

#####################################################################

#####################################################################


################
# Loading data #
################

dir = "path to your data"
setwd(dir)

ImprintSeq_Cal_methCpG = read.table(file = paste("ImprintSeq_Cal_methCpG", "txt", sep = "."), header = T, sep = "\t", quote = "", stringsAsFactors = F, check.names = F)
ImprintSeq_methCpG = read.table(file = paste("ImprintSeq_methCpG", "txt", sep = "."), header = T, sep = "\t", quote = "", stringsAsFactors = F, check.names = F)

dir = "path to your data"
setwd(dir)

ImprintSeq_Samples = read.table(file = paste("ImprintSeq_Samples", "txt", sep = "."), header = T, sep = "\t", quote = "", stringsAsFactors = F, check.names = T)
ImprintSeq_Targets = read.table(file = paste("ImprintSeq_Targets", "txt", sep = "."), header = T, sep = "\t", quote = "", stringsAsFactors = F, check.names = T)
ImprintSeq_Targets_Info = read.table(file = paste("ImprintSeq_Targets_Info", "txt", sep = "."), na.strings = c("", " ", "NA"), header = T, sep = "\t", quote = "", stringsAsFactors = F, check.names = T)

SubPheno = read.table(file = paste("SubPheno", "txt", sep = "."), na.strings = c("", " ", "NA"), header = T, sep = "\t", quote = "", stringsAsFactors = F, check.names = T)

#####################################################################

#####################################################################


###########
# Loading #
###########

dir = "path to your output folder"
dir = paste(paste(dir, paste(SD, "SD", sep = ""), sep = ""), "\\", sep = "")
setwd(dir)

load(file = paste(paste("MethylCal_data_corr", paste(Calibration_list, collapse = "_"), sep = "_"), "RData", sep = "."))
load(file = paste(paste("MethylCal_HyperHypomethylated", paste(Calibration_list, collapse = "_"), sep = "_"), "RData", sep = "."))
load(file = paste(paste("MethylCal_HyperHypomethylated_sign", paste(Calibration_list, collapse = "_"), sep = "_"), "RData", sep = "."))
load(file = paste(paste("MethylCal_HyperHypomethylated_value", paste(Calibration_list, collapse = "_"), sep = "_"), "RData", sep = "."))

#####################################################################

#####################################################################


####################
# Removing Targets #
####################

idx_targets = match(ImprintSeq_Targets_Info$Target, ImprintSeq_Targets$Target)
ImprintSeq_Targets = ImprintSeq_Targets[idx_targets, ]

#####################################################################

#####################################################################

utargets = unique(ImprintSeq_Targets$Target)
n_utargets = length(utargets)

utargets_OLD = utargets
n_utargets_OLD = n_utargets

ocal = c("CS0_L1",  "CS50_L1",  "CS100_L1", 
         "CS0_L2",  "CS50_L2",  "CS100_L2", 
         "CS0_L3",  "CS50_L3",  "CS100_L3", 
         "CS0o_L4", "CS50o_L4", "CS100o_L4", 
         "CS0n_L4", "CS50n_L4", "CS100n_L4", 
         "CS0_L5",  "CS50_L5",  "CS100_L5")

DNAset = c("Old", "old", "Old", 
           "New", "New", "New", 
           "New", "",    "New", 
           "Old", "Old", "Old", 
           "New", "New", "New", 
           "Old", "",    "")

ImprintSeq_Cal_methCpG_OLD = ImprintSeq_Cal_methCpG
ImprintSeq_Cal_methCpG = NULL

for (Calibration.ID in Calibration_list)
{
  
  if (Calibration.ID == "L1")
  {
    idx = which(ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[1] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[2] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[3])
  } else if (Calibration.ID == "L2") {
    idx = which(ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[4] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[5] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[6])
  } else if (Calibration.ID == "L3") {
    idx = which(ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[7] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[8] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[9])
  } else if (Calibration.ID == "L4o") {
    idx = which(ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[10] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[11] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[12])
  } else if (Calibration.ID == "L4n") {
    idx = which(ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[13] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[14] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[15])
  } else if (Calibration.ID == "L5") {
    idx = which(ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[16] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[17] | ImprintSeq_Cal_methCpG_OLD$Calibration_ID == ocal[18])
  }
  
  if (is.null(ImprintSeq_Cal_methCpG))
  {
    ImprintSeq_Cal_methCpG = ImprintSeq_Cal_methCpG_OLD[idx, -c(2, 3)]   # ***
  } else {
    ImprintSeq_Cal_methCpG = cbind(ImprintSeq_Cal_methCpG, ImprintSeq_Cal_methCpG_OLD[idx, -c(1 : 5)])   # ***
  }

}

n_ImprintSeq_Cal_methCpG = nrow(ImprintSeq_Cal_methCpG)
p_ImprintSeq_Cal_methCpG = ncol(ImprintSeq_Cal_methCpG)
n_ImprintSeq_methCpG = nrow(ImprintSeq_methCpG)
p_ImprintSeq_methCpG = ncol(ImprintSeq_methCpG)

idx = grep("C02", colnames(ImprintSeq_methCpG))
c_idx = setdiff((1 : p_ImprintSeq_methCpG), idx)   # ***
n_Control = length(idx)
n_Case = length(c_idx)
print(n_Control)
print(n_Case)
c_idx = c_idx[match(SubPheno[, 1], colnames(ImprintSeq_methCpG[, c_idx]))]
n_Case = length(c_idx)
print(n_Case)
ImprintSeq_methCpG = ImprintSeq_methCpG[, c(idx, c_idx)]

ImprintSeq_data = cbind(ImprintSeq_Cal_methCpG, ImprintSeq_methCpG)

#####################################################################

#####################################################################

col_names = colnames(ImprintSeq_methCpG)
Table_MethylCal_HyperHypomethylated_sign = matrix(NA, (n_Control + n_Case), n_utargets)

for (t in 1 : n_utargets)
{
  t_tmp = idx_targets[t]
  HyperHypomethylated = unlist(MethylCal_HyperHypomethylated[[t_tmp]])
  idx = match(HyperHypomethylated, col_names)
  MethylCal_Hyper_idx = unlist(MethylCal_HyperHypomethylated_sign[[t_tmp]]) == "+"
  MethylCal_Hypo_idx = unlist(MethylCal_HyperHypomethylated_sign[[t_tmp]]) == "-"
  Table_MethylCal_HyperHypomethylated_sign[idx[MethylCal_Hyper_idx], t] = 1
  Table_MethylCal_HyperHypomethylated_sign[idx[MethylCal_Hypo_idx], t] = -1
}

utargets = ImprintSeq_Targets_Info$Target.official.names
rownames(Table_MethylCal_HyperHypomethylated_sign) = col_names
colnames(Table_MethylCal_HyperHypomethylated_sign) = utargets

Table_MethylCal_HyperHypomethylated_value = matrix(NA, (n_Control + n_Case), n_utargets)

for (t in 1 : n_utargets)
{
  t_tmp = idx_targets[t]
  HyperHypomethylated = unlist(MethylCal_HyperHypomethylated[[t_tmp]])
  idx = match(HyperHypomethylated, col_names)
  idx_idx = !is.na(idx)
  Table_MethylCal_HyperHypomethylated_value[idx[idx_idx], t] = unlist(MethylCal_HyperHypomethylated_value[[t_tmp]])[idx_idx]
}

utargets = ImprintSeq_Targets_Info$Target.official.names
rownames(Table_MethylCal_HyperHypomethylated_value) = col_names
colnames(Table_MethylCal_HyperHypomethylated_value) = utargets

#####################################################################

#####################################################################

idx_C = grep("C02", colnames(ImprintSeq_methCpG))
FP = sum(rowSums(!is.na(Table_MethylCal_HyperHypomethylated_sign))[idx_C])
FPR = FP / (length(idx_C) * length(utargets))

TP = sum(rowSums(!is.na(Table_MethylCal_HyperHypomethylated_sign)))
FDR = FP / TP

print(round(c(FPR, (1 - pnorm(SD)) * 2), 4))
print(round(c(FPR, (1 - pnorm(SD) + (1 - pnorm(SD)) /2)), 4))
print(round(c(FP, TP, FDR), 4))

#####################################################################

#####################################################################

idx_Control = intersect(grep("C02", colnames(data_corr)), grep("corrected", colnames(data_corr)))
idx_Case = setdiff(grep("corrected", colnames(data_corr)), idx_Control)
Control_corrected = matrix(0, length(idx_Control), n_utargets)
Control_corrected_average = matrix(0, 1, n_utargets)
Case_corrected = matrix(0, length(idx_Case), n_utargets)

for (t in 1 : n_utargets)
{
  idx_targets = which(data_corr$Target == utargets_OLD[t])
  Control_corrected[, t] = apply(data_corr[idx_targets, idx_Control], 2, median, na.rm = T)
  Control_corrected_average[t] = mean(apply(data_corr[idx_targets, idx_Control], 2, median, na.rm = T))
  Case_corrected[, t] = apply(data_corr[idx_targets, idx_Case], 2, median, na.rm = T)
}

Control_corrected_average = Control_corrected_average + eps
Case_corrected = Case_corrected + eps
Table_MethylCal_HyperHypomethylated_value = Table_MethylCal_HyperHypomethylated_value + eps

Table_MethylCal_HyperHypomethylated_DAML = Table_MethylCal_HyperHypomethylated_value

for (t in 1 : n_utargets)
{
   Table_MethylCal_HyperHypomethylated_DAML[, t] = round(abs(Table_MethylCal_HyperHypomethylated_value[, t] - Control_corrected[t]), 3)
}

rownames(Control_corrected) = col_names[c(1 : length(idx_Control))]
colnames(Control_corrected) = utargets
colnames(Control_corrected_average) = utargets
rownames(Case_corrected) = col_names[-c(1 : length(idx_Control))]
colnames(Case_corrected) = utargets

#####################################################################

#####################################################################


##########
# Saving #
##########

dir = "path to your output folder"
dir = paste(paste(dir, paste(SD, "SD", sep = ""), sep = ""), "\\", sep = "")
setwd(dir)

write.table(Control_corrected, file = paste(paste("Table_MethylCal_Controls", paste(SD, "SD", sep = ""), sep = "_"), "xls", sep = "."), quote = F, row.names = T, col.names = T, sep ="\t")
write.table(Control_corrected_average, file = paste(paste("Table_MethylCal_Controls_Average", paste(SD, "SD", sep = ""), sep = "_"), "xls", sep = "."), quote = F, row.names = F, col.names = T, sep ="\t")
write.table(Case_corrected, file = paste(paste("Table_MethylCal_Cases", paste(SD, "SD", sep = ""), sep = "_"), "xls", sep = "."), quote = F, row.names = T, col.names = T, sep ="\t")
write.table(Table_MethylCal_HyperHypomethylated_DAML, file = paste(paste("Table_MethylCal_HyperHypomethylated_DAML", paste(SD, "SD", sep = ""), sep = "_"), "xls", sep = "."), quote = F, row.names = T, col.names = T, sep ="\t")

#####################################################################
