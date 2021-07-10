## Last Modification 7/6/2021
## Initialize the Parameters

Init_Parameters<-function(){
  CANCER_TYPE="ACC"
  ALL_GIs=c("TP53","KRAS","NRAS")  ## read from a text file
  
  Par<-list(
    CANCER=CANCER_TYPE,  ### Name of the study
    RNAFileName_Pan=paste0("./InputData/Pan_RNASEQV2_HGNC_",CANCER_TYPE,".RData"),
    MethFileName_Pan=paste0("./InputData/Pan_Methyl_",CANCER_TYPE,".RData"),
    CNVFileName_Pan=paste0("./InputData/Pan_CNV_",CANCER_TYPE,".RData"),
    Mutation_Pan=paste0("./InputData/Mut_WO_Hyper_HGNC_",CANCER_TYPE,".RData"),
    ALL_GIs=ALL_GIs
  )
  return(Par)
}

