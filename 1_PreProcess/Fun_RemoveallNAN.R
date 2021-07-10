# This function remove all nan genes.
Zscore<-function(myRow){
  mu <- mean(myRow,na.rm = TRUE) 
  sigma   <- sqrt (var(myRow,na.rm=TRUE))
  myRow <- (myRow - mu )/ sigma
}

Fun_RemoveallNAN<-function(Mutation_Pan){
  Epsilon=1e-04
  load("./Data/ReadRNASeq.RData")
  load(Mutation_Pan)

  Normal_Samples=colnames(GenRNASeqMat0)
  # ## Normal Samples
  if(ncol(GenRNASeqMat0)){
    GenRNASeqMat0=GenRNASeqMat0[rowSums(is.na(GenRNASeqMat0)) != ncol(GenRNASeqMat0), ] ## Remove all NA
    Rem_NanIndex0=which(is.na(GenRNASeqMat0),arr.ind=TRUE)
    GenRNASeqMat0[Rem_NanIndex0]<-0
    Min_GenRNA0=min(GenRNASeqMat0)
    GenRNASeqMat0=GenRNASeqMat0-Min_GenRNA0
    GenRNASeqMat0=log2(GenRNASeqMat0+Epsilon)
    if(!is.null(dim(GenRNASeqMat0)[2])){
      GenRNASeqMat0=t(apply(GenRNASeqMat0,1,Zscore))
    }
    GenRNASeqMat0[which(is.nan(GenRNASeqMat0))]=0  ## should replace with the log of zero+epsilon
  }
  GenRNASeqMat0=as.data.frame(GenRNASeqMat0)
  colnames(GenRNASeqMat0)=Normal_Samples
  
  ############################ Tumor Samples
  
  ###### Samples and Patients of MAF file
  Mutation_Samples=data.frame(FullName_Sample_MAF=unique(Mutt_Current_Tissue$Tumor_Sample_Barcode))
  Mutation_Samples$MAF_Sample_Name=substr(Mutation_Samples$FullName_Sample_MAF,1,16)
  Mutation_Samples$MAF_Patient_Name=substr(Mutation_Samples$FullName_Sample_MAF,1,12)
  Mutation_Samples = Mutation_Samples[!duplicated(Mutation_Samples$MAF_Sample_Name),]
  
  #######Samples and Patients from RNA-Seq File
  RNA_Samples=data.frame(Fullname_RNA_Samples=colnames(GenRNASeqMat1))
  RNA_Samples$RNA_Patient_Name=substr(RNA_Samples$Fullname_RNA_Samples,1,12)
  RNA_Samples$RNA_Sample_Name=substr(RNA_Samples$Fullname_RNA_Samples,1,16)
  RNA_Samples=RNA_Samples[!duplicated(RNA_Samples$RNA_Sample_Name),]
  
  Overlapped_Samples=merge(RNA_Samples,Mutation_Samples,by.x="RNA_Patient_Name",by.y ="MAF_Patient_Name",all.x = T,all.y = T)
  Overlapped_Samples=Overlapped_Samples[complete.cases(Overlapped_Samples), ]
  Overlapped_Samples=Overlapped_Samples[!duplicated(Overlapped_Samples$RNA_Patient_Name),]
  
  GenRNASeqMat1=GenRNASeqMat1[,which(colnames(GenRNASeqMat1)%in%Overlapped_Samples$Fullname_RNA_Samples)]  
  
  #################################
  ## Temp Tumor samples
  if(ncol(GenRNASeqMat1)){
    GenRNASeqMat1=GenRNASeqMat1[rowSums(is.na(GenRNASeqMat1)) != ncol(GenRNASeqMat1), ]
    Rem_NanIndex1=which(is.na(GenRNASeqMat1),arr.ind=TRUE)
    GenRNASeqMat1[Rem_NanIndex1]<-0
    Min_GenRNA1=min(GenRNASeqMat1)
    GenRNASeqMat1=GenRNASeqMat1-Min_GenRNA1
    GenRNASeqMat1=log2(GenRNASeqMat1+Epsilon)
    GenRNASeqMat1=t(apply(GenRNASeqMat1, 1,Zscore))
    GenRNASeqMat1[which(is.nan(GenRNASeqMat1))]=0  ## should replace with the log of zero+epsilon
  }
  
  save(GenRNASeqMat0,GenRNASeqMat1,file="./Data/FinalRNASeq.RData")
  
  ##
  load("./Data/ReadMethyl_MeanMaxMin.RData") ## Here i am calculating the NAN remove just for mean
  
  RemoveallNAN_Methylation<-function(GenMethylMat0,GenMethylMat1){
    GenMethylMat0=GenMethylMat0[(rowSums(is.na(GenMethylMat0)) != ncol(GenMethylMat0)), ]
    GenMethylMat1=GenMethylMat1[(rowSums(is.na(GenMethylMat1)) != ncol(GenMethylMat1)), ]
    
    Rem_NanIndex0=which(is.na(GenMethylMat0),arr.ind=TRUE)
    Rem_NanIndex1=which(is.na(GenMethylMat1),arr.ind=TRUE)
    
    GenMethylMat0[Rem_NanIndex0]<-0
    GenMethylMat1[Rem_NanIndex1]<-0
    
    NewGenMethylMat0=GenMethylMat0
    NewGenMethylMat1=GenMethylMat1
    
    Metheyl<-list(NewGenMethylMat0,GenMethylMat1)
    return(Metheyl)
  }
  
  ## Remove nan from Methylation data
  Methyl_MeanMat0=Meth_Mean[[1]]
  Methyl_MeanMat1=Meth_Mean[[2]]
  Metheyl_Mean=RemoveallNAN_Methylation(Methyl_MeanMat0,Methyl_MeanMat1)
  
  Methyl_MaxMat0=Meth_Max[[1]]
  Methyl_MaxMat1=Meth_Max[[2]]
  Metheyl_Max=RemoveallNAN_Methylation(Methyl_MaxMat0,Methyl_MaxMat1)
  
  Methyl_MinMat0=Meth_Min[[1]]
  Methyl_MinMat1=Meth_Min[[2]]
  Metheyl_Min=RemoveallNAN_Methylation(Methyl_MinMat0,Methyl_MinMat1)
  
  save(Metheyl_Mean,Metheyl_Max,Metheyl_Min,file="./Data/FinalMethyl_MeanMaxMin.RData")
  
  ####
  load("./Data/ReadCNV.RData")
  GenCNVMat0=GenCNVMat0[rowSums(is.na(GenCNVMat0)) != ncol(GenCNVMat0), ]
  GenCNVMat1=GenCNVMat1[rowSums(is.na(GenCNVMat1)) != ncol(GenCNVMat1), ]
  
  ## Estimating the value for the elements that we don't have information based on the other near genes.
  
  GenCNVMat0=na.locf(GenCNVMat0)
  GenCNVMat1=na.locf(GenCNVMat1)
  
  ## there is another function  $$na.approx(m)$$ which estimate a linear interpolation.
  
  save(GenCNVMat0,GenCNVMat1,file="./Data/FinalCNV.RData")
}


