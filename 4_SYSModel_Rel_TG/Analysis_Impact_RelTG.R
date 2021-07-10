## Last Modification 6/6/2020

Analysis_Impact_RelTG<-function(){
  ######### Thresholds 
  Number_Sigma=1
  Curr_BhtThresh=0.0001
  Curr_BhtRatio=3
  Pval_Threshold=0.001  ### threshold for the result of wilcoxen test p-value of NM-M and control
  Mutated_Threshod=2    ## less than 2 sample should  be filtered
  TG_Numb_Threshold=8   ## minimum number of TGs for analysis 
  TG_Numb_Threshold_NewFilter=5   ## minimum number of TGs for analysis 
  Percent_Thresh_NewFilter=1
  
  ######### Read the Raw-Mutation Rate
  Raw_Mut_File <-paste0("../1_PreProcess/InputData/MutRt_Raw_",CANCER_TYPE,".txt")
  Raw_MutRate=read.delim(Raw_Mut_File,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  ########
  filename=paste0("./Results/Res_NonSyn_Wild_",CANCER_TYPE,".txt")
  Result=read.delim(filename,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Result=Result[!duplicated(Result),]
  Result$Percentage_Run=Result$Mut_Rate_GI_Run
  colnames(Result)[colnames(Result)=="Conn_level"]="Conn_Level"
  
  ######### Ratio of Bhattacharyya
  Result$Mean_Sigma_Cont=NA
  Result$Mean_Sigma_Cont=Result$Bht_C_Mean+(Number_Sigma*sqrt(Result$Bht_C_Var))  ### 
  
  Result$Ratio_Bha=NA
  Result$Ratio_Bha=(Result$Bht_Nm_M/Result$Mean_Sigma_Cont) ## 
  
  ################# Rule1 (If WPval is greater than a threshold then remove the connection)
  Result$Rule1a_absNsigma=NA
  WI_Impact_Idx=which(Result$WPval<Pval_Threshold)
  Result$Rule1a_absNsigma[WI_Impact_Idx]=1
  NA_IDX=which(is.na(Result$Rule1a_absNsigma))
  Result$Rule1a_absNsigma[NA_IDX]=0
  
  ################# Rule 2 (Filter out Noise)
  Result$Rule1b_Bht_NMM=NA
  Result$Rule1b_Bht_NMM=ifelse(Result$Bht_Nm_M>Curr_BhtThresh,1,0)  
  
  ################ Rule 3, Rule to detect high Vs.low Impact-only when rule to detect impact is true-Bht_Ratio>3 ifabs(Bhtratio>3)1 else 0
  Result$Rule2_Imp_Detection=NA
  Result$Rule2_Imp_Detection=(Result$Rule1a_absNsigma*Result$Rule1b_Bht_NMM)
  Impacted_Connections=which(Result$Rule2_Imp_Detection==1)
  
  Result$Rule2_HighLow_Imp=NA 
  Result$Rule2_HighLow_Imp[Impacted_Connections]=ifelse(Result$Ratio_Bha[Impacted_Connections]>Curr_BhtRatio,2,1)
  Result$Rule2_HighLow_Imp[which(is.na(Result$Rule2_HighLow_Imp))]=0
  
  ######### Add Raw-MutRates to the Result matrix
  Result_GI=Result[(which(is.na(Result$Hyper_GI))),]
  Result_HyperGI=Result[(which(!is.na(Result$Hyper_GI))),]
  
  Results_GI_WI_MutRat=merge(Result_GI,Raw_MutRate,by.x="GI",by.y = "Gene_Name",all.x = TRUE)
  Result_HyperGI_WI_MutRat=merge(Result_HyperGI,Raw_MutRate,by.x="Hyper_GI",by.y = "Gene_Name",all.x = TRUE)
  Result=rbind(Result_HyperGI_WI_MutRat,Results_GI_WI_MutRat)
  
  filename_Res=paste0("./Results/Res_NonSyn_Wild_Impacts_",CANCER_TYPE,".txt")
  write.table(Result,file=filename_Res,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ############### Find all GIs level one and level Two Which is euqal to All_GIs that we have.
  HyperGI=unique(Result$GI[which(is.na(Result$Hyper_GI))])
  GI_wt_Hyp_Idx=which(!is.na(Result$Hyper_GI))
  GI_wt_Hyper=unique(Result$Hyper_GI[GI_wt_Hyp_Idx])
  ALL_GI=unique(c(GI_wt_Hyper,HyperGI))
  
  ############### 
  K=0
  RES=array(list())
  
  for(GIIdx in 1:length(ALL_GI)){
    ############# TG Direct
    Current_GI=ALL_GI[GIIdx] 
    GI_Level1=which((Result$GI==Current_GI)&(is.na(Result$Hyper_GI)))
    GI_Level1_Complx=which((Result$Hyper_GI==Current_GI)&(Result$Conn_Level==1))
    Level1_IDX=c(GI_Level1,GI_Level1_Complx)
    if(length(Level1_IDX)){
      TG_Direct_Info=Result[Level1_IDX,]
      WO_Impct_B=which(TG_Direct_Info$Rule2_Imp_Detection==0)
      WI_Impct_B=which(TG_Direct_Info$Rule2_Imp_Detection==1)
      ## Percentage of Impacted
      TG_Dir_WO_B=length(WO_Impct_B)
      TG_Dir_WI_B=length(WI_Impct_B)
      Perce_Impacted_Dir=(TG_Dir_WI_B/(TG_Dir_WI_B+TG_Dir_WO_B))*100
      Perce_UnImpacted_Dir=(TG_Dir_WO_B/(TG_Dir_WI_B+TG_Dir_WO_B))*100
      ## High and low Impacts of level one
      WO_Impct_B1=which(TG_Direct_Info$Rule2_HighLow_Imp==0)
      WL_Impct_B1=which(TG_Direct_Info$Rule2_HighLow_Imp==1)
      WH_Impct_B1=which(TG_Direct_Info$Rule2_HighLow_Imp==2)
      ## Percentage of High impact
      TG_Dir_WO_B1=length(WO_Impct_B1)
      TG_Dir_WL_B1=length(WL_Impct_B1)
      TG_Dir_WH_B1=length(WH_Impct_B1)
      ##Percentage of Dir High impact
      Perce_Impacted_Dir_High1=(TG_Dir_WH_B1/(TG_Dir_WH_B1+TG_Dir_WO_B1+TG_Dir_WL_B1))*100
      ##Percentage of Dir Low impact
      Perce_Impacted_Dir_Low1=(TG_Dir_WL_B1/(TG_Dir_WH_B1+TG_Dir_WO_B1+TG_Dir_WL_B1))*100
      ##Percentage of NonSyn_Wild
      Percentage_NonSy_Wild=TG_Direct_Info$Percentage_Run[1]
      ##Mutrateof Universe
      MutRate_Univers=TG_Direct_Info$Mut_Rate_GI_Universe[1] 
      ##Raw_Mut_Rate
      Raw_MutRate=TG_Direct_Info$MutRat_Raw[1]
      ##No. Mutated Samples
      No_Mut_Samples=TG_Direct_Info$NonSyn_Label[1]
      ## No. Mutated Universe
      No_Mut_Samp_Univers=TG_Direct_Info$Mutated_Sampl_Unniverse[1]
      ## All Samples that run
      All_Samples=TG_Direct_Info$All_label[1]
      ## All Samples Universe
      All_Samples_Univers=TG_Direct_Info$All_label_Universe[1]
    }else{
      TG_Dir_WO_B=0
      TG_Dir_WI_B=0
      Perce_Impacted_Dir=0
      Perce_UnImpacted_Dir=0
      TG_Dir_WO_B1=0
      TG_Dir_WH_B1=0
      TG_Dir_WL_B1=0 
      Perce_Impacted_Dir_High1=0
      Perce_Impacted_Dir_Low1=0
      Percentage_NonSy_Wild=0
      MutRate_Univers=0
      Raw_MutRate=0
      No_Mut_Samples=0
      No_Mut_Samp_Univers=0
      All_Samples=0
      All_Samples_Univers=0
    }
    ############# TG Indirect
    Level2_IDX=which((Result$Hyper_GI==Current_GI)&(Result$Conn_Level==2))
    if(length(Level2_IDX)){
      
      TG_InDirect_Info=Result[Level2_IDX,]
      ###
      WO_Impct_B2=which(TG_InDirect_Info$Rule2_Imp_Detection==0)
      WI_Impct_B2=which(TG_InDirect_Info$Rule2_Imp_Detection==1)
      
      ## Percentage of Impacted
      TG_InDir_WO_B=length(WO_Impct_B2)
      TG_InDir_WI_B=length(WI_Impct_B2)
      Perce_Impacted_InDir=(TG_InDir_WI_B/(TG_InDir_WI_B+TG_InDir_WO_B))*100 
      Perce_UnImpacted_InDir=(TG_InDir_WO_B/(TG_InDir_WI_B+TG_InDir_WO_B))*100 
      
      ## High and low Impacts of level one
      WO_Impct_B21=which(TG_InDirect_Info$Rule2_HighLow_Imp==0)
      WL_Impct_B2=which(TG_InDirect_Info$Rule2_HighLow_Imp==1)
      WH_Impct_B2=which(TG_InDirect_Info$Rule2_HighLow_Imp==2)
      ##
      ## Percentage of High impact
      TG_InDir_WO_B2=length(WO_Impct_B21)
      TG_InDir_WL_B2=length(WL_Impct_B2)
      TG_InDir_WH_B2=length(WH_Impct_B2)
      ##Percentage of High impact
      Perce_Impacted_InDir_High2=(TG_InDir_WH_B2/(TG_InDir_WH_B2+TG_InDir_WL_B2+TG_InDir_WO_B2))*100
      ##Percentage of Low impact
      Perce_Impacted_InDir_Low2=(TG_InDir_WL_B2/(TG_InDir_WH_B2+TG_InDir_WL_B2+TG_InDir_WO_B2))*100
      Percentage_NonSy_Wild=TG_InDirect_Info$Percentage_Run[1]
      
      MutRate_Univers=TG_InDirect_Info$Mut_Rate_GI_Universe[1] 
      
      Raw_MutRate=TG_InDirect_Info$MutRat_Raw[1]
      
      No_Mut_Samp_Univers=TG_InDirect_Info$Mutated_Sampl_Unniverse[1]
      
      No_Mut_Samples=TG_InDirect_Info$NonSyn_Label[1]
      ## All Samples that run
      All_Samples=TG_InDirect_Info$All_label[1]
      
      All_Samples_Univers=TG_InDirect_Info$All_label_Universe[1]
      
    }else{
      TG_InDir_WO_B=0
      TG_InDir_WI_B=0
      Perce_Impacted_InDir=0
      Perce_UnImpacted_InDir=0
      TG_InDir_WO_B2=0
      TG_InDir_WL_B2=0
      TG_InDir_WH_B2=0
      Perce_Impacted_InDir_High2=0
      Perce_Impacted_InDir_Low2=0
    }
    
    ####################################
    ### TG_All(Direct and Indirect)_Bha  
    TG_All_WO_B=TG_Dir_WO_B+TG_InDir_WO_B
    TG_All_WI_B=TG_Dir_WI_B+TG_InDir_WI_B
    Percent_Imp_Final_B=(TG_All_WI_B/(TG_All_WI_B+TG_All_WO_B))*100
    Percent_UnImp_Final_B=(TG_All_WO_B/(TG_All_WI_B+TG_All_WO_B))*100
    
    ####################################
    TG_All_WO_B12=TG_Dir_WO_B1+TG_InDir_WO_B2
    TG_All_WH_B12=TG_Dir_WH_B1+TG_InDir_WH_B2
    TG_All_WL_B12=TG_Dir_WL_B1+TG_InDir_WL_B2
    Percent_HighImp_Final_B=(TG_All_WH_B12/(TG_All_WH_B12+TG_All_WL_B12+TG_All_WO_B12))*100
    Percent_LowImp_Final_B=(TG_All_WL_B12/(TG_All_WH_B12+TG_All_WL_B12+TG_All_WO_B12))*100
    
    GI_TG_total=TG_All_WO_B+TG_All_WI_B
    K=K+1
    
    ANAlysis_Inf0=list(
      GI_ALL=Current_GI,
      Percentage_NonSy_Wild=Percentage_NonSy_Wild,
      MutRate_Univers=MutRate_Univers,
      Raw_MutRate=Raw_MutRate,
      No_Mut_Samples=No_Mut_Samples,
      No_Mut_Samp_Univers=No_Mut_Samp_Univers,
      All_Samples=All_Samples,
      All_Samples_Univers=All_Samples_Univers,
      Perce_Impact_Dir=Perce_Impacted_Dir,
      Perce_Impact_InDir=Perce_Impacted_InDir,
      Perce_Impact_Final=Percent_Imp_Final_B,
      TG_Dir_WI=TG_Dir_WI_B,
      TG_Dir_WO=TG_Dir_WO_B,
      TG_InDir_WI=TG_InDir_WI_B,
      TG_InDir_WO=TG_InDir_WO_B,
      TG_All_WI=TG_All_WI_B,
      TG_All_WO=TG_All_WO_B,
      ##########
      Perce_Dir_High=Perce_Impacted_Dir_High1,
      Perce_Dir_Low=Perce_Impacted_Dir_Low1,
      
      Perce_InDir_High=Perce_Impacted_InDir_High2,
      Perce_InDir_Low=Perce_Impacted_InDir_Low2,
      
      Percent_All_High=Percent_HighImp_Final_B,
      Percent_All_Low=Percent_LowImp_Final_B,
      
      TG_Dir_WH=TG_Dir_WH_B1,
      TG_Dir_WL=TG_Dir_WL_B1,
      TG_InDir_WH=TG_InDir_WH_B2,
      TG_InDir_WL=TG_InDir_WL_B2,
      ##
      TG_All_WH=TG_All_WH_B12,
      TG_All_WL=TG_All_WL_B12,
      TG_All2_WO=TG_All_WO_B12,
      
      GI_TG_total=GI_TG_total
      ##
    )
    RES[[K]]=ANAlysis_Inf0
  }
  
  ###############
  datalist = list()
  M=0
  for(i in 1:length(RES)){
    GI_temp=RES[[i]]
    Data_Analysis=do.call(cbind.data.frame, GI_temp)
    M=M+1
    datalist[[M]] <- Data_Analysis
  }
  Big_data = do.call(rbind, datalist)
  
  Summary_Name=paste0("./Analysis/Res_NonSyn_Wild_Summry_GI_",CANCER_TYPE,".txt")
  write.table(Big_data,file=Summary_Name,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ################ Filterout based Bhat
  Filtered_BigData_B=Big_data
  
  Filtered_BigData_B=Filtered_BigData_B[which(Filtered_BigData_B$No_Mut_Samples>=Mutated_Threshod),]
  Filtered_BigData_B=Filtered_BigData_B[which((Filtered_BigData_B$TG_All_WI+Filtered_BigData_B$TG_All_WO)>=TG_Numb_Threshold),]  ## minimum number of TGs 
  Filt_Sum_Name=paste0("./Analysis/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG.txt")
  write.table(Filtered_BigData_B,file=Filt_Sum_Name,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ###########################
  Result=Result[which((Result$Hyper_GI%in%Filtered_BigData_B$GI_ALL)|(is.na(Result$Hyper_GI))&(Result$GI%in%Filtered_BigData_B$GI_ALL)),]
  filename_Res=paste0("./Analysis/Res_NonSyn_Wild_Impacts_Filtered_RelTG_",CANCER_TYPE,".txt")
  write.table(Result,file=filename_Res,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ###########################
  Impcts=list(
    WO_Imp=length(which(Result$Rule2_HighLow_Imp==0)),
    WH_Imp=length(which(Result$Rule2_HighLow_Imp==2)),
    WL_Imp=length(which(Result$Rule2_HighLow_Imp==1)))
  
  filename_Imp=paste0("./Analysis/Res_NonSyn_Wild_Summary_Impacts_",CANCER_TYPE,".txt")
  write.table(Impcts,file=filename_Imp,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ##########################################################################
  ######################### Filterout Based on new Rules
  
  ################ Filterout based Bhat
  Filtered_BigData_c=Big_data
  
  Filtered_BigData_c=Filtered_BigData_c[which(Filtered_BigData_c$No_Mut_Samples>=Mutated_Threshod),]  ####Filter1=
  Filtered_BigData_c=Filtered_BigData_c[which((Filtered_BigData_c$GI_TG_total)>=TG_Numb_Threshold_NewFilter),]  ## minimum number of TGs 
  
  Index_of_Filterd_Elements=which((Filtered_BigData_c$No_Mut_Samples==Mutated_Threshod)&(Filtered_BigData_c$Percentage_NonSy_Wild<Percent_Thresh_NewFilter))
  if(length(Index_of_Filterd_Elements)){
    Filtered_BigData_c=Filtered_BigData_c[-Index_of_Filterd_Elements,]  
  }
  
  Filt_Sum_Name=paste0("./Analysis/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt")
  write.table(Filtered_BigData_c,file=Filt_Sum_Name,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ###########################
  Result_New=Result[which((Result$Hyper_GI%in%Filtered_BigData_c$GI_ALL)|(is.na(Result$Hyper_GI))&(Result$GI%in%Filtered_BigData_c$GI_ALL)),]
  filename_Res=paste0("./Analysis/Res_NonSyn_Wild_Impacts_Filtered_RelTG_NewFilter_",CANCER_TYPE,".txt")
  write.table(Result_New,file=filename_Res,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ###########################
  Impcts=list(
    WO_Imp=length(which(Result_New$Rule2_HighLow_Imp==0)),
    WH_Imp=length(which(Result_New$Rule2_HighLow_Imp==2)),
    WL_Imp=length(which(Result_New$Rule2_HighLow_Imp==1)))
  
  filename_Imp=paste0("./Analysis/Res_NonSyn_Wild_Summary_Impacts_NewFilter_",CANCER_TYPE,".txt")
  write.table(Impcts,file=filename_Imp,append=FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
}