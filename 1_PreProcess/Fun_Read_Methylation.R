## Read the Methyl Features Min Max Mean 
# Last Update 9/14/2020

Fun_Read_Methylation<-function (MethFileName){
  
  Methyl_Feat_names<-load(MethFileName)
  
  source("Fun_Methyl_Feat.R")
  for(FeatIndx in 1:length(Methyl_Feat_names)){
    if (Methyl_Feat_names[FeatIndx]=="Min_Methy_F"){
      Meth_Min<-Fun_Methyl_Feat(Min_Methy_F)
    }else if(Methyl_Feat_names[FeatIndx]=="Max_Methy_F"){
      Meth_Max<-Fun_Methyl_Feat(Max_Methy_F)
    }else{
      Meth_Mean<-Fun_Methyl_Feat(Mean_Methy_F)
    }
  }
  save(Meth_Min,Meth_Max,Meth_Mean,file="./Data/ReadMethyl_MeanMaxMin.Rdata")
}