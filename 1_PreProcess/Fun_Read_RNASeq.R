## Read the RNAseq files  
# Last Update 9/14/2020
Fun_Read_RNASeq<-function (RNAFileName){
  
  ## Read RNASeq data 
  load(RNAFileName)
  
  Solid_Tumor_Sample_All=RNASeq2[,which(substr(colnames(RNASeq2),14,15)=="01")]
  Match_Norma_Sample_All=RNASeq2[,which(substr(colnames(RNASeq2),14,15)%in%c("11","10"))]

  GenRNASeqMat0=Match_Norma_Sample_All
  GenRNASeqMat1=Solid_Tumor_Sample_All
  save(GenRNASeqMat0,GenRNASeqMat1,file="./Data/ReadRNASeq.RData")
}

