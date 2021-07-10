## Last Modification 2/5/2020
## Initialize the Parameteres

Parameters<-function(){
  cores=detectCores()

  if(.Platform$OS.type == "unix") { ##Unix Machine
    Path=getwd()
    Splited_path=strsplit(Path,"/")
    Subgroup=tail(unlist(Splited_path), n=1)
    Subgroup_Leng=1
    } else {  ## Windows Machine
    Subgroup=1
    Subgroup_Leng=1
  }
  Par<-list(
    cores=cores, 
    Subgroup=as.numeric(Subgroup),
    Subgroup_Leng=Subgroup_Leng,
    MutFile_Name=paste0("../1_PreProcess/Data/Updated_Mut_",CANCER_TYPE,".RData"),  ## This file is created based on list of the RNAsaq samples.
    RNA_File_Name=paste0("../1_PreProcess/Data/Combin_Gen_Samp_RNA_",CANCER_TYPE,".RData"),
    Trans_File_Name=paste0("../1_PreProcess/Data/Combin_Gen_Samp_Trans_",CANCER_TYPE,".RData"),
    Sig_File_Name=paste0("../1_PreProcess/Data/Combin_Gen_Samp_Sig_",CANCER_TYPE,".RData"),
    SYS_Net="../1_PreProcess/Data/SYS_Net.RData",
    GI_TG_PG="../2_Findconnection_Randomly/Results/Final_GI_TG_PG.RData",
    GGI_Net_Info="../2_Findconnection_Randomly/Results/GGI_Net_Info.RData",
    WO_Zscore_File_Name=paste0("../1_PreProcess/Data/FinalRNASeq_WO_Zscor_",CANCER_TYPE,".RData")
  )
  return(Par)
}

