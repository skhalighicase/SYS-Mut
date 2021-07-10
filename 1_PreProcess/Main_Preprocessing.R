## Main Preprocessing 7/7/2021

graphics.off() 
rm(list=ls()) 

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TCGA2STAT")

## Set Path
CurrentDir=getwd()
setwd(CurrentDir)

## Library
library(readxl)
library(dplyr)
library(zoo)
library(stringi)
library(stringr)
library(RTCGAToolbox)

#--------------- Parameters
source("Init_Parameters.R")
Params<-Init_Parameters()
#--------------- 
CANCER_TYPE<<-Params$CANCER

##--------------- Read RNASeq, Methyl, CNA
source("Fun_Read_RNASeq.R")
Fun_Read_RNASeq(Params$RNAFileName_Pan)

## Read Methyl Featrues Per Cancer
source("Fun_Read_Methylation.R")
Fun_Read_Methylation(Params$MethFileName_Pan)

## Read CNV data of Pan_cancer
source("Fun_Read_CNV.R")
Fun_Read_CNV(Params$CNVFileName_Pan)

## Read Influence Network (New Network integrated from Super Pathway)
source("Fun_InfluenceNet.R")
Fun_InfluenceNet()

## Remove All NaN with out applying Z-Score and Calculate the Mean/Var/Max for the Expression values
source("Fun_RemoveallNAN_RNA_WO_Zscore.R")
Fun_RemoveallNAN_RNA_WO_Zscore(Params$Mutation_Pan)

## Remove All NaN
source("Fun_RemoveallNAN.R")
Fun_RemoveallNAN(Params$Mutation_Pan)

## Find the CNV and methylation for the Genes that have expression
source("CommonGenes_CNV_Methyl.R")
CommonGenes_CNV_Methyl()

## Find Common samplesre
source("CommonSamples.R")
CommonSamples(Params$Mutation_Pan)

## Apply Transformation to the CNA and Methyl Data
source("TrasnformCNA_Meth.R")
TrasnformCNA_Meth()

## Combine all data
source("Combine_AllData_WO_Zscor.R")
Combine_AllData_WO_Zscor()
