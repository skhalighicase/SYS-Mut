Average_RnG_WI_NewFilter_logRank<-function(){
  ############ Iterations of random replacement of the network
  Max_Iteration=10
  
  #############
  Res_RelTG=read.delim(Params$filename_Rel,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Mutated_GI <-as.character(Res_RelTG$GI_ALL[Res_RelTG$Perce_Impact_Final>0])
  
  Res_RelALLTG=read.delim(Params$filename_Rel_Can,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  ############ Random_TG
  Res_RndTG=read.delim(Params$filename_Rand,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  for(Fil_Idx in 2:Max_Iteration){
    filename=paste0("../5_SYSModel_Rnd_TG/RUNS/RUN_",Fil_Idx,"/Analysis/Results_NonSyn_Wild_Impacts_RUN_",Fil_Idx,".txt")
    Res_Temp_RndTG=read.delim(filename,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
    Curr_Res=dplyr::select(Res_Temp_RndTG,Hyper_GI,GI,Old_TG,TG,LogLiklhd)
    names(Curr_Res)[c(4,5)] <- c(paste0("TG_",Fil_Idx),paste0("LogLiklhd_",Fil_Idx))
    Res_RndTG=merge(Res_RndTG,Curr_Res,by = c("Hyper_GI","GI","Old_TG"),)
  }
  Res_RndTG$Avg_Log <- rowMeans(subset(Res_RndTG, select = c(LogLiklhd,LogLiklhd_2,LogLiklhd_3,LogLiklhd_4,LogLiklhd_5,LogLiklhd_6,LogLiklhd_7,LogLiklhd_8,LogLiklhd_9,LogLiklhd_10)), na.rm = TRUE)
  Res_RelALLTG$combine = paste(Res_RelALLTG$Hyper_GI, Res_RelALLTG$GI, Res_RelALLTG$TG,sep = "_")
  Res_RndTG$combine = paste(Res_RndTG$Hyper_GI, Res_RndTG$GI, Res_RndTG$Old_TG,sep = "_")
  Res_RelALLTG=Res_RelALLTG[which(Res_RelALLTG$combine%in%intersect(Res_RelALLTG$combine,Res_RndTG$combine)),]
  Res_RelALLTG$LogLiklhd=-log10(Res_RelALLTG$WPval)
  
  Res_RndTG=Res_RndTG[which(Res_RndTG$combine%in%intersect(Res_RelALLTG$combine,Res_RndTG$combine)),]
  setnames(Res_RndTG, old = c("TG","LogLiklhd"), new = c("TG_1","LogLiklhd_1"))

  Res_RndTG_Selected=dplyr::select(Res_RndTG,c("Hyper_GI","GI","Old_TG","TG_1","LogLiklhd_1","TG_2","LogLiklhd_2","TG_3","LogLiklhd_3","TG_4","LogLiklhd_4","TG_5","LogLiklhd_5","TG_6","LogLiklhd_6"
                                        ,"TG_7","LogLiklhd_7","TG_8","LogLiklhd_8","TG_9","LogLiklhd_9","TG_10","LogLiklhd_10"))
  
  ################
  All_Log_Likelihoods=merge(Res_RelALLTG,Res_RndTG_Selected,by.x = c("Hyper_GI","GI","TG"), by.y = c("Hyper_GI","GI","Old_TG"))
  
  Comparison_LogLike=dplyr::select(All_Log_Likelihoods,c("Hyper_GI","GI","TG","Parent_Genes","Conn_Level","NonSyn_Label","Wild_Label","All_label","Mut_Rate_GI_Universe","Mut_Rate_GI_Run","LogLiklhd","LogLiklhd_1","LogLiklhd_2","LogLiklhd_3","LogLiklhd_4","LogLiklhd_5","LogLiklhd_6","LogLiklhd_7","LogLiklhd_8","LogLiklhd_9","LogLiklhd_10"))
  Comp=subset(Comparison_LogLike,select = c("LogLiklhd","LogLiklhd_1","LogLiklhd_2","LogLiklhd_3","LogLiklhd_4","LogLiklhd_5","LogLiklhd_6","LogLiklhd_7","LogLiklhd_8","LogLiklhd_9","LogLiklhd_10"))
  
  ##############
  Comparison_LogLike$P_values_Comp_Rel_Rand_GI_TG=apply(Comp,1, function(A) {
    wilcox.test(as.numeric(A[2:11]),NULL,alternative = "less",mu = as.numeric(A[1]),conf.level = 0.95)$p.value
  })
  Comparison_LogLike$Impact_Scor_GI_TG= -log10(Comparison_LogLike$P_values_Comp_Rel_Rand_GI_TG)
  
  Output_filename=paste0("./Results/GI_To_TG_ImpScore_Pval_",CANCER_TYPE,".csv")
  write.csv(Comparison_LogLike, Output_filename,row.names = FALSE)
}