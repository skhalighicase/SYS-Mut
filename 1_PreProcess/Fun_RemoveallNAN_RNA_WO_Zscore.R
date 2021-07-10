# Zscore<-function(myRow){
#   Variance<- var(myRow,na.rm=TRUE)
#   return(Variance)
# }

#########
## This function doesnot apply Z-score. It first remove NA values and after that calculate Var mean and Max for the expression data
Fun_RemoveallNAN_RNA_WO_Zscore<-function(Mutation_Pan){
  Epsilon=1e-06
  load("./Data/ReadRNASeq.RData")
  load(Mutation_Pan)
  
  # # ## Normal Samples   
  # if(dim(GenRNASeqMat0)[2]){
  #   GenRNASeqMat0=GenRNASeqMat0[rowSums(is.na(GenRNASeqMat0)) != ncol(GenRNASeqMat0), ] ## Remove all NA
  #   Rem_NanIndex0=which(is.na(GenRNASeqMat0),arr.ind=TRUE)
  #   GenRNASeqMat0[Rem_NanIndex0]<-0
  #   Min_GenRNA0=min(GenRNASeqMat0)
  #   GenRNASeqMat0=GenRNASeqMat0-Min_GenRNA0
  #   GenRNASeqMat0=log2(GenRNASeqMat0+Epsilon)
  #   Var_Genes0=apply(GenRNASeqMat0,1,var)
  #   Mean_Genes0=apply(GenRNASeqMat0,1,mean)
  #   Max_Genes0=apply(GenRNASeqMat0,1,max)
  #   GenRNASeqMat0[which(is.nan(GenRNASeqMat0))]=0  ## should replace with the log of zero+epsilon
  # }
  ###### Samples and Patients of MAF file
  Mutation_Samples=data.frame(FullName_Sample_MAF=unique(Mutt_Current_Tissue$Tumor_Sample_Barcode))
  Mutation_Samples$MAF_Sample_Name=substr(Mutation_Samples$FullName_Sample_MAF,1,16)
  Mutation_Samples$MAF_Patient_Name=substr(Mutation_Samples$FullName_Sample_MAF,1,12)
  Mutation_Samples = Mutation_Samples[!duplicated(Mutation_Samples$MAF_Sample_Name),]
  
  #######Samples and Patients from RNA-Seq File
  RNA_Samples=data.frame(Fullname_RNA_Samples=colnames(GenRNASeqMat1))
  RNA_Samples$RNA_Patient_Name=substr(RNA_Samples$Fullname_RNA_Samples,1,12)
  RNA_Samples$RNA_Sample_Name=substr(RNA_Samples$Fullname_RNA_Samples,1,16)
  RNA_Samples = RNA_Samples[!duplicated(RNA_Samples$RNA_Sample_Name),]

  Overlapped_Samples=merge(RNA_Samples,Mutation_Samples,by.x="RNA_Patient_Name",by.y ="MAF_Patient_Name",all.x = T,all.y = T)
  Overlapped_Samples=Overlapped_Samples[complete.cases(Overlapped_Samples), ]
  Overlapped_Samples=Overlapped_Samples[!duplicated(Overlapped_Samples$RNA_Patient_Name),]
  
  GenRNASeqMat1=GenRNASeqMat1[,which(colnames(GenRNASeqMat1)%in%Overlapped_Samples$Fullname_RNA_Samples)]  

  ###### Temp Tumor samples
  if(ncol(GenRNASeqMat1)){
    GenRNASeqMat1=GenRNASeqMat1[rowSums(is.na(GenRNASeqMat1)) != ncol(GenRNASeqMat1), ]
    Rem_NanIndex1=which(is.na(GenRNASeqMat1),arr.ind=TRUE)
    GenRNASeqMat1[Rem_NanIndex1]<-0
    Min_GenRNA1=min(GenRNASeqMat1)
    GenRNASeqMat1=GenRNASeqMat1-Min_GenRNA1
    GenRNASeqMat1=log2(GenRNASeqMat1+Epsilon)
    Var_Genes=apply(GenRNASeqMat1,1,var)
    Mean_Genes=apply(GenRNASeqMat1,1,mean)
    Max_Genes=apply(GenRNASeqMat1,1,max)
    # GenRNASeqMat1=t(apply(GenRNASeqMat1, 1,Zscore))
    GenR=as.matrix(GenRNASeqMat1)
    GenR[which(is.nan(GenR))]=0 
    GenRNASeqMat1=as.data.frame(GenR)  ## should replace with the log of zero+epsilon
  }
  AllGene_Inf_CancerSamp=cbind(Var_Genes,Mean_Genes,Max_Genes)
  FileName=paste0("./Data/FinalRNASeq_WO_Zscor_",CANCER_TYPE,".RData")
  save(AllGene_Inf_CancerSamp,file=FileName)
}