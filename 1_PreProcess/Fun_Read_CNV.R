## Read the sCNA files 
## Last Update 9/14/2020

Fun_Read_CNV<-function (CNVFileName){
  
  load(CNVFileName)
  
  Solid_Tumor_Sample_All=CNV_RES[,which(substr(colnames(CNV_RES),14,15)=="01")]
  Match_Norma_Sample_All=CNV_RES[,which(substr(colnames(CNV_RES),14,15)%in%c("11","10"))]
  
  GenCNVMat0=Match_Norma_Sample_All
  GenCNVMat1=Solid_Tumor_Sample_All
  
  save(GenCNVMat0,GenCNVMat1,file="./Data/ReadCNV.RData")
}