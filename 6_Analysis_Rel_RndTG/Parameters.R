Parameters<-function(){
    Par<-list(
    CANCER=CANCER_TYPE,  ### Name of the study
    filename_Rel="../4_SYSModel_Rel_TG/Analysis/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt",
    filename_Rel_Can=paste0("../4_SYSModel_Rel_TG/Analysis/Res_NonSyn_Wild_Impacts_Filtered_RelTG_NewFilter_",CANCER_TYPE,".txt"),
    filename_Rand="../5_SYSModel_Rnd_TG/RUNS/RUN_1/Analysis/Results_NonSyn_Wild_Impacts_RUN_1.txt",
    filename_Rand_Summ="../5_SYSModel_Rnd_TG/RUNS/RUN_1/Analysis/Results_NonSyn_Wild_Summry_GI_RUN_1.txt"
  )
  return(Par)
}

