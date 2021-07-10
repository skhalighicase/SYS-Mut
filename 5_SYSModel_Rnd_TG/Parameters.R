## Last Modification 6/6/2021
## Initialize the Parameters

Parameters<-function(){
  cores=detectCores()
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
    Iteration=1  ###### You need to change this value for each iteration
    cores=cores, 
    Subgroup=as.numeric(Subgroup),
    Subgroup_Leng=Subgroup_Leng,
    CANCER=CANCER_TYPE,  ### Name of the study
    MutFile_Name=paste0("../1_PreProcess/Data/Updated_Mut_",CANCER_TYPE,".RData"),  ## This file is created based on list of the RNAsaq samples.
    RNA_File_Name=paste0("../1_PreProcess/Data/Combin_Gen_Samp_RNA_",CANCER_TYPE,".RData"),
    Trans_File_Name=paste0("../1_PreProcess/Data/Combin_Gen_Samp_Trans_",CANCER_TYPE,".RData"),
    Sig_File_Name=paste0("../1_PreProcess/Data/Combin_Gen_Samp_Sig_",CANCER_TYPE,".RData"),
    SYSMut_InfNet="../1_PreProcess/Data/SYS_Net.RData",
    GGI_Net_Info="../2_Findconnection_Randomly/Results/GGI_Net_Info.RData",
    TG_PG_Connections ="../2_Findconnection_Randomly/Results/TG_PGs_Connec.RData",
    Best_Parents_File=paste0("../3_FindBestParents/Results/Results_TG_Parents_",CANCER_TYPE,".RData"),
    Rel_GI_Summry="../4_SYSModel_Rel_TG/Analysis/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt"
  )
  return(Par)
}

