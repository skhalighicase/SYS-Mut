## Last Modification 2/5/2020

Normal_Ele<-function(TGs,RNA_GI,GI,Hyper_GI){
  TGs=intersect(TGs,Selec_GI_TGs)
  Selec_GI_TGs<<-Selec_GI_TGs[!Selec_GI_TGs %in% TGs]
  ## Initialization
  Seen_TGs=which(!is.na(match(TGs,Flags_Seen)))
  if(length(Seen_TGs)){
    TGs=TGs[-Seen_TGs]
    if (length(TGs)){
      Flags_Seen<<-c(Flags_Seen,TGs)
    }
  }else{
    Flags_Seen<<-c(Flags_Seen,TGs)
  }
  #####################
  ## Expression of Parents and CNV Methylation of Target Genes 
  if (length(TGs)){
    Result_TG=array(list(),length(TGs))
    for(TGIdx in 1:length(TGs)){
      ########### 
      TG<-TGs[TGIdx] ## RNA of the Current Target Gene 
      source("CVN_Meth_Parent_Calcu.R")
      Final_Data<-CVN_Meth_Parent_Calcu(TG,RNA_GI,GI,Hyper_GI)
      ####################################
      source("SYS_Model_oneCluster.R")
      PostMutation2<-SYS_Model_oneCluster(Final_Data$CNV_Meth_Par_U,Final_Data$RNA_GI_U,Final_Data$RNA_TG_U,Final_Data$PG_Weights,Final_Data$Miss_Com_Indx_U,Final_Data$weightsGI)
      ###################################
      Analysis_TG<-list(
        Hyper_GI=Hyper_GI,
        GI=GI,
        TG=TG,
        PG=as.character(Final_Data$PG),
        Res_NM_M=PostMutation2,
        NonSyn_Label=Final_Data$Non_Syn_length,
        Wild_Label=Final_Data$Wild_type_Length
      )
      Result_TG[[TGIdx]]=Analysis_TG
    }
    return(Result_TG)
  }else{
    return(NULL)
  }
}

