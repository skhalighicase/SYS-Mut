## This function finds the cnv and methylation data of the genes that has expression
CommonGenes_CNV_Methyl<-function(){
  
  load("./Data/FinalCNV.RData")
  load("./Data/FinalMethyl_MeanMaxMin.RData")
  load("./Data/FinalRNASeq.RData")

  ###### 
  if(ncol(GenRNASeqMat0)){
    GeneNames_RNASeq0 <-rownames(GenRNASeqMat0)
    GeneNames_CNV0 <-rownames(GenCNVMat0)
    GeneNames_Methyl0 <-rownames(Metheyl_Mean[[1]])
    
    ## Common Genes of CNV and Methyl
    Com_CNV_RNA0<-which(GeneNames_CNV0%in%GeneNames_RNASeq0)
    Com_methyl_RNA0<-which(GeneNames_Methyl0%in%GeneNames_RNASeq0)
    
    ## Common CNV0
    Comm_CNV0=GenCNVMat0[Com_CNV_RNA0,]
    
    ## RNASeq
    Comm_RNAS0=GenRNASeqMat0
    
    ## Methy Mean Max Min 0
    Comm_Meth_Mean0=data.frame(Metheyl_Mean[[1]][Com_methyl_RNA0,])
    colnames(Comm_Meth_Mean0)=gsub("\\.","-",colnames(Comm_Meth_Mean0))
    Comm_Meth_Min0=data.frame(Metheyl_Min[[1]][Com_methyl_RNA0,])
    colnames(Comm_Meth_Min0)=gsub("\\.","-",colnames(Comm_Meth_Min0))
    Comm_Meth_Max0=data.frame(Metheyl_Max[[1]][Com_methyl_RNA0,])
    colnames(Comm_Meth_Max0)=gsub("\\.","-",colnames(Comm_Meth_Max0))
  }else{
    Com_CNV_RNA0=matrix(NA,0,0)
    Com_methyl_RNA0=matrix(NA,0,0)
    Comm_CNV0=matrix(NA,0,0)
    Comm_RNAS0=matrix(NA,0,0)
    Comm_Meth_Mean0=matrix(NA,0,0)
    Comm_Meth_Min0=matrix(NA,0,0)
    Comm_Meth_Max0=matrix(NA,0,0)
  }
  
  
  ######### For Cancer Samples ################
  GeneNames_RNASeq1 <-rownames(GenRNASeqMat1)
  GeneNames_CNV1 <-rownames(GenCNVMat1)
  GeneNames_Methyl1 <-rownames(Metheyl_Mean[[2]])
  
  ######### Find Common Genes between Mut and RNASeq
  
  ## Common Genes of CNV and Methyl
  Com_CNV_RNA1<-which(GeneNames_CNV1%in%GeneNames_RNASeq1)
  Com_methyl_RNA1<-which(GeneNames_Methyl1%in%GeneNames_RNASeq1)
  
  ## Common CNV1
  Comm_CNV1=GenCNVMat1[Com_CNV_RNA1,]
  
  ## Methy Mean Max Min 1
  Comm_Meth_Mean1=data.frame(Metheyl_Mean[[2]][Com_methyl_RNA1,])
  colnames(Comm_Meth_Mean1)=gsub("\\.","-",colnames(Comm_Meth_Mean1))
  Comm_Meth_Min1=data.frame(Metheyl_Min[[2]][Com_methyl_RNA1,])
  colnames(Comm_Meth_Min1)=gsub("\\.","-",colnames(Comm_Meth_Min1))
  Comm_Meth_Max1=data.frame(Metheyl_Max[[2]][Com_methyl_RNA1,])
  colnames(Comm_Meth_Max1)=gsub("\\.","-",colnames(Comm_Meth_Max1))
  
  ## RNASeq
  Comm_RNAS1=GenRNASeqMat1
  FileName=paste0("./Data/CommonGeneData_",CANCER_TYPE,".RData")
  save(Com_CNV_RNA0,Com_methyl_RNA0,Com_CNV_RNA1,Com_methyl_RNA1,Comm_CNV0,Comm_CNV1,Comm_Meth_Mean0,Comm_Meth_Mean1,Comm_Meth_Min0,Comm_Meth_Min1,Comm_Meth_Max0,Comm_Meth_Max1,Comm_RNAS0,Comm_RNAS1,file=FileName)
}