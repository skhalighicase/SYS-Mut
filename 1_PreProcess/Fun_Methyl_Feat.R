Fun_Methyl_Feat<-function(Methyl){
  ##  Read Methyl data 
  Solid_Tumor_Sample_All=Methyl[,which(substr(colnames(Methyl),14,15)=="01")]
  Match_Norma_Sample_All=Methyl[,which(substr(colnames(Methyl),14,15)%in%c("11","10"))]
  
  GenMethylMat0=Match_Norma_Sample_All
  GenMethylMat1=Solid_Tumor_Sample_All
  
  Gen_Meth<-list(GenMethylMat0,GenMethylMat1)
  return(Gen_Meth)
}

