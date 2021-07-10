## ZScore

############
Combine_AllData_WO_Zscor<-function(){
  Tran_CNVFileName=paste0("./Data/Transformed_CNV_Meth_",CANCER_TYPE,".RData")
  Comm_File_Name=paste0("./Data/CommonGenesamples_",CANCER_TYPE,".RData")
  load (Tran_CNVFileName)
  load (Comm_File_Name)
  
  ## Load Influence Network (New Influence Network)
  load("./Data/SYS_Net.RData")
  ##  Common Sample Gene Data 0, Replace NaN with zero.
  if(!is.null(dim(Comm_Sam_RNAS0)[2])){
    if(dim(Comm_Sam_RNAS0)[2]){
      Nor_Comm_CNV0=Comm_Sam_CNV0#t(apply(Comm_Sam_CNV0, 1,Zscore))
      Nor_Comm_CNV0[which(is.na(Nor_Comm_CNV0),arr.ind=TRUE)]<-0
      Nor_Comm_Meth_Max0=Comm_Sam_Meth_Max0#t(apply(Comm_Sam_Meth_Max0, 1,Zscore))
      Nor_Comm_Meth_Max0[which(is.na(Nor_Comm_Meth_Max0),arr.ind=TRUE)]<-0
      Nor_Comm_Meth_Min0=Comm_Sam_Meth_Min0#t(apply(Comm_Sam_Meth_Min0, 1,Zscore))
      Nor_Comm_Meth_Min0[which(is.na(Nor_Comm_Meth_Min0),arr.ind=TRUE)]<-0
      Nor_Comm_Meth_Mean0=Comm_Sam_Meth_Mean0#t(apply(Comm_Sam_Meth_Mean0, 1,Zscore))
      Nor_Comm_Meth_Mean0[which(is.na(Nor_Comm_Meth_Mean0),arr.ind=TRUE)]<-0
      Nor_Comm_RNA0=Comm_Sam_RNAS0#t(apply(Comm_Sam_RNAS0, 1,Zscore))
      Nor_Comm_RNA0[which(is.na(Nor_Comm_RNA0),arr.ind=TRUE)]<-0
      
      ###################################################################
      ## Remove Unnecessary Genes that we don't have in influence Network
      ######################
      Init_Genes_N=rownames(Nor_Comm_RNA0)
      Exist_GI_N=match(Init_Genes_N,Influ_net$GeneA)
      Nor_Comm_RNA0_GI=Nor_Comm_RNA0[-which(is.na(Exist_GI_N)),]
      
      Nor_Comm_Data0<-list(
        Nor_Comm_CNV0=Nor_Comm_CNV0,
        Nor_Comm_Meth_Max0=Nor_Comm_Meth_Max0,
        Nor_Comm_Meth_Min0=Nor_Comm_Meth_Min0,
        Nor_Comm_Meth_Mean0=Nor_Comm_Meth_Mean0,
        Nor_Comm_RNA0=Nor_Comm_RNA0,
        Nor_Comm_RNA0_GI=Nor_Comm_RNA0_GI
      )
      ######################
      ##  Transform data sigmoid
      Nor_Sig_CNV0=SigmoidCNV$SigmoidCNV0#t(apply(SigmoidCNV$SigmoidCNV0, 1,Zscore))
      Nor_Sig_CNV0[which(is.na(Nor_Sig_CNV0),arr.ind=TRUE)]<-0
      Nor_Sig_Meth_Max0=SigmoidMeth$SigmoidMeth_Max0#t(apply(SigmoidMeth$SigmoidMeth_Max0, 1,Zscore))
      Nor_Sig_Meth_Max0[which(is.na(Nor_Sig_Meth_Max0),arr.ind=TRUE)]<-0
      Nor_Sig_Meth_Min0=SigmoidMeth$SigmoidMeth_Min0#t(apply(SigmoidMeth$SigmoidMeth_Min0, 1,Zscore))
      Nor_Sig_Meth_Min0[which(is.na(Nor_Sig_Meth_Min0),arr.ind=TRUE)]<-0
      Nor_Sig_Meth_Mean0=SigmoidMeth$SigmoidMeth_Mean0#t(apply(SigmoidMeth$SigmoidMeth_Mean0, 1,Zscore))
      Nor_Sig_Meth_Mean0[which(is.na(Nor_Sig_Meth_Mean0),arr.ind=TRUE)]<-0
      Nor_ThreshCNV0=ThreshCNV$ThreshCNV0#t(apply(ThreshCNV$ThreshCNV0, 1,Zscore))
      Nor_ThreshCNV0[which(is.na(Nor_ThreshCNV0),arr.ind=TRUE)]<-0
      Nor_ThreshMeth_Max0=ThreshMeth$ThreshMeth_Max0#t(apply(ThreshMeth$ThreshMeth_Max0, 1,Zscore))
      Nor_ThreshMeth_Max0[which(is.na(Nor_ThreshMeth_Max0),arr.ind=TRUE)]<-0
      Nor_ThreshMeth_Min0=ThreshMeth$ThreshMeth_Min0#t(apply(ThreshMeth$ThreshMeth_Min0, 1,Zscore))
      Nor_ThreshMeth_Min0[which(is.na(Nor_ThreshMeth_Min0),arr.ind=TRUE)]<-0
      Nor_ThreshMeth_Mean0=ThreshMeth$ThreshMeth_Mean0#t(apply(ThreshMeth$ThreshMeth_Mean0, 1,Zscore))
      Nor_ThreshMeth_Mean0[which(is.na(Nor_ThreshMeth_Mean0),arr.ind=TRUE)]<-0
      
      Nor_Comm_Trans_Data0<-list(
        Nor_Sig_CNV0=Nor_Sig_CNV0,
        Nor_Sig_Meth_Max0=Nor_Sig_Meth_Max0,
        Nor_Sig_Meth_Min0=Nor_Sig_Meth_Min0,
        Nor_Sig_Meth_Mean0=Nor_Sig_Meth_Mean0,
        Nor_ThreshCNV0=Nor_ThreshCNV0,
        Nor_ThreshMeth_Max0=Nor_ThreshMeth_Max0,
        Nor_ThreshMeth_Min0=Nor_ThreshMeth_Min0,
        Nor_ThreshMeth_Mean0=Nor_ThreshMeth_Mean0
      )
    }else{
      Nor_Comm_Data0=list(NA)
      Nor_Comm_Trans_Data0=list(NA)
    }
  }else{
    Nor_Comm_Data0=list(NA)
    Nor_Comm_Trans_Data0=list(NA)
  }
  
  Combine_File_Name=paste0("./Data/Combin_Gen_Samp_Normal_",CANCER_TYPE,".RData")
  save(Nor_Comm_Data0,Nor_Comm_Trans_Data0,file=Combine_File_Name)
  rm(Nor_Comm_Data0,Nor_Comm_Trans_Data0)
  
  #############################################################
  ## Common Sample Gene Data 1, Replace NaN with zero. 
  
  Nor_Comm_CNV1=Comm_Sam_CNV1#t(apply(Comm_Sam_CNV1, 1,Zscore))
  Nor_Comm_CNV1[which(is.na(Nor_Comm_CNV1),arr.ind=TRUE)]<-0
  Nor_Comm_Meth_Max1=Comm_Sam_Meth_Max1#t(apply(Comm_Sam_Meth_Max1,1,Zscore))
  Nor_Comm_Meth_Max1[which(is.na(Nor_Comm_Meth_Max1),arr.ind=TRUE)]<-0
  Nor_Comm_Meth_Min1=Comm_Sam_Meth_Min1#t(apply(Comm_Sam_Meth_Min1,1,Zscore))
  Nor_Comm_Meth_Min1[which(is.na(Nor_Comm_Meth_Min1),arr.ind=TRUE)]<-0
  Nor_Comm_Meth_Mean1=Comm_Sam_Meth_Mean1#t(apply(Comm_Sam_Meth_Mean1,1,Zscore))
  Nor_Comm_Meth_Mean1[which(is.na(Nor_Comm_Meth_Mean1),arr.ind=TRUE)]<-0
  Nor_Comm_RNA1=Comm_Sam_RNAS1#t(apply(Comm_Sam_RNAS1,1,Zscore))
  Nor_Comm_RNA1[which(is.na(Nor_Comm_RNA1),arr.ind=TRUE)]<-0
  
  ###################################################################
  ## Remove Unnecessary Genes that we don't have in influence Network
  ######################
  Init_Genes=rownames(Nor_Comm_RNA1)
  Exist_GI=match(Init_Genes,Influ_net$GeneA)
  Nor_Comm_RNA1_GI=Nor_Comm_RNA1[-which(is.na(Exist_GI)),]
  
  #####################
  Nor_Comm_Data1<-list(
    Nor_Comm_RNA1=Nor_Comm_RNA1,
    Nor_Comm_RNA1_GI=Nor_Comm_RNA1_GI
  )
  Combine_File_Name2=paste0("./Data/Combin_Gen_Samp_RNA_",CANCER_TYPE,".RData")
  save(Nor_Comm_Data1,file=Combine_File_Name2)
  
  Nor_Comm_Meth1_CNV1<-list(
    Nor_Comm_CNV1=Nor_Comm_CNV1,
    Nor_Comm_Meth_Max1=Nor_Comm_Meth_Max1,
    Nor_Comm_Meth_Min1=Nor_Comm_Meth_Min1,
    Nor_Comm_Meth_Mean1=Nor_Comm_Meth_Mean1
    
  )
  Combine_File_Name3=paste0("./Data/Combin_Gen_Samp_CNV_Methyl_",CANCER_TYPE,".RData")
  save(Nor_Comm_Meth1_CNV1,file=Combine_File_Name3)
  
  #####################
  Nor_Sig_CNV1=SigmoidCNV$SigmoidCNV1#t(apply(SigmoidCNV$SigmoidCNV1, 1,Zscore))
  Nor_Sig_CNV1[which(is.na(Nor_Sig_CNV1),arr.ind=TRUE)]<-0
  Nor_Sig_Meth_Max1=SigmoidMeth$SigmoidMeth_Max1#t(apply(SigmoidMeth$SigmoidMeth_Max1, 1,Zscore))
  Nor_Sig_Meth_Max1[which(is.na(Nor_Sig_Meth_Max1),arr.ind=TRUE)]<-0
  Nor_Sig_Meth_Min1=SigmoidMeth$SigmoidMeth_Min1#t(apply(SigmoidMeth$SigmoidMeth_Min1, 1,Zscore))
  Nor_Sig_Meth_Min1[which(is.na(Nor_Sig_Meth_Min1),arr.ind=TRUE)]<-0
  Nor_Sig_Meth_Mean1=SigmoidMeth$SigmoidMeth_Mean1#t(apply(SigmoidMeth$SigmoidMeth_Mean1, 1,Zscore))
  Nor_Sig_Meth_Mean1[which(is.na(Nor_Sig_Meth_Mean1),arr.ind=TRUE)]<-0
  
  Nor_Comm_Sig_Data1<-list(
    Nor_Sig_CNV1=Nor_Sig_CNV1,
    Nor_Sig_Meth_Max1=Nor_Sig_Meth_Max1,
    Nor_Sig_Meth_Min1=Nor_Sig_Meth_Min1,
    Nor_Sig_Meth_Mean1=Nor_Sig_Meth_Mean1
  )
  Comb_Samp_Name4=paste0("./Data/Combin_Gen_Samp_Sig_",CANCER_TYPE,".RData")
  save(Nor_Comm_Sig_Data1,file=Comb_Samp_Name4)
  rm(Nor_Comm_Sig_Data1,Nor_Sig_Meth_Mean1,Nor_Sig_Meth_Max1,Nor_Sig_CNV1,Nor_Sig_Meth_Min1)
  
  Nor_ThreshCNV1=ThreshCNV$ThreshCNV1#t(apply(ThreshCNV$ThreshCNV1,1,Zscore))
  Nor_ThreshCNV1[which(is.na(Nor_ThreshCNV1),arr.ind=TRUE)]<-0
  Nor_ThreshMeth_Max1=ThreshMeth$ThreshMeth_Max1#t(apply(ThreshMeth$ThreshMeth_Max1, 1,Zscore))
  Nor_ThreshMeth_Max1[which(is.na(Nor_ThreshMeth_Max1),arr.ind=TRUE)]<-0
  Nor_ThreshMeth_Min1=ThreshMeth$ThreshMeth_Min1#t(apply(ThreshMeth$ThreshMeth_Min1, 1,Zscore))
  Nor_ThreshMeth_Min1[which(is.na(Nor_ThreshMeth_Min1),arr.ind=TRUE)]<-0
  Nor_ThreshMeth_Mean1=ThreshMeth$ThreshMeth_Mean1#t(apply(ThreshMeth$ThreshMeth_Mean1, 1,Zscore))
  Nor_ThreshMeth_Mean1[which(is.na(Nor_ThreshMeth_Mean1),arr.ind=TRUE)]<-0
  
  Nor_Comm_Trans_Data1<-list(
    Nor_ThreshCNV1=Nor_ThreshCNV1,
    Nor_ThreshMeth_Max1=Nor_ThreshMeth_Max1,
    Nor_ThreshMeth_Min1=Nor_ThreshMeth_Min1,
    Nor_ThreshMeth_Mean1=Nor_ThreshMeth_Mean1
  )
  Comb_Samp_Name5=paste0("./Data/Combin_Gen_Samp_Trans_",CANCER_TYPE,".RData")
  save(Nor_Comm_Trans_Data1,file=Comb_Samp_Name5)
}