## Last Modification 6/6/2021
## Initialize the Parameters

Parameters<-function(){
  cores=detectCores()
  CANCER_TYPE="ACC"
  ALL_GIs=c("TP53","KRAS","NRAS")
  if(.Platform$OS.type == "unix") { ##Unix Machine
    Path=getwd()
    Splited_path=strsplit(Path,"/")
    Subgroup=tail(unlist(Splited_path), n=1)
    Subgroup_Leng=1
    } else { ## Windows Machine
    Subgroup=1
    Subgroup_Leng=1
  }
  Par<-list(
    cores=cores, 
    ALL_GIs=ALL_GIs,
    Subgroup=as.numeric(Subgroup),
    Subgroup_Leng=Subgroup_Leng,
    CANCER=CANCER_TYPE,  ### Name of the study
    MutFile_Name=paste0("../Data/Updated_Mut_",CANCER_TYPE,".RData"),  ## This file is created based on list of the RNAsaq samples.
    RNA_File_Name=paste0("../Data/Combin_Gen_Samp_RNA_",CANCER_TYPE,".RData"),
    Trans_File_Name=paste0("../Data/Combin_Gen_Samp_Trans_",CANCER_TYPE,".RData"),
    Sig_File_Name=paste0("../Data/Combin_Gen_Samp_Sig_",CANCER_TYPE,".RData"),
    SYSMut_InfNet="../Data/SYS_Net.RData",
    GGI_Net_Info="../Data/GGI_Net_Info.RData",
    Best_Parents_File=paste0("../Data/Results_TG_Parents_",CANCER_TYPE,".RData")
  )
  return(Par)
}

