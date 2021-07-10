## Lat Update 6/9/2020
## Analysis of the result 

All_GI_Imp_ECDF_RelTG_TrueVal<-function(){
  
  #############
  filename_Rel=paste0("../4_SYSModel_Rel_TG/Analysis/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt")
  Res_RelTG=read.delim(filename_Rel,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Mutated_GI <-as.character(Res_RelTG$GI_ALL[Res_RelTG$Perce_Impact_Final>0])
  
  Res_RelTG_For_ECDS=Res_RelTG$Perce_Impact_Final[Res_RelTG$Perce_Impact_Final>0]
  #########
  par(mar=c(1,1,1,1))
  # create a vector of quantiles
  quants <-seq(0,1,length=81)[2:80]
  
  
  ecdf_Estimated = ecdf(Res_RelTG_For_ECDS)  
  Res_RelTG$Probability=ecdf_Estimated(Res_RelTG$Perce_Impact_Final)
  #Res_GI_WI_Prob=subset(Res_RelTG,select = c(GI_ALL,Perce_Impact_Final,Probability))
  
  write.csv(Res_RelTG,file="./Results/Res_GIImp_WI_Prob_ECDF_RelTG_RelVal.csv",row.names = F)
}