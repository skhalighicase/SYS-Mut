## Last Modification 5/31/2021
#### This function find all the required data for the Model.

Normal_Ele_Rnd<-function(OldTGs,RNA_GI,GI,Hyper_GI){
  library(fpc)
  Seen_OldTGs=which(OldTGs%in%Flags_Seen_OldTG)
  if(length(Seen_OldTGs)){
    OldTGs=OldTGs[-Seen_OldTGs]
    if (length(OldTGs)){
      Flags_Seen_OldTG<<-c(Flags_Seen_OldTG,OldTGs)
    }
  }else{
    Flags_Seen_OldTG<<-c(Flags_Seen_OldTG,OldTGs)
  }
  ###### Find Random Targets ######
  if(length(OldTGs)){
    source("Gene_Replace.R")
    #TGs<-Gene_Replace(OldTGs,GI,Hyper_GI)
    TGsAll<-Gene_Replace(OldTGs,GI,Hyper_GI)
    TGs<-TGsAll$TGs
  }else{
    TGs<-OldTGs
    TGsAll=data.frame(TGs=OldTGs,Old_TGs=OldTGs)
  }
  ##################### 
  ## Expression of Parents and CNV Methylation of Target Genes 
  if (length(TGs)){
    Result_TG=array(list(),length(TGs))
    for(TGIdx in 1:length(TGs)){ 
      ########### 
      TG<-TGs[TGIdx]  ## RNA of the Current Target Gene 
      TG_Par_Info=All_Info_TG_Par[which(All_Info_TG_Par$TG%in%TG),]
      ##########################
      ## We will update the list of the samples for control. We eliminate the samples that have mutation on all related genes. parents and TG and GI. 
      SamplesAll<-data.frame("Samples_Name"=rownames(as.data.frame(RNA_GI)))
      ## Mutation Info ## A function to read and Lable the Missense Mutations
      source("Mutation_Labels_Clean_Wildtype.R")
      Mutation_Info<-Mutation_Labels_Clean_Wildtype(GI,Hyper_GI,TG_Par_Info,TG,SamplesAll) 
      Mutation_Rate_RUN=(length(which(Mutation_Info$Mut_Info==2))/nrow(Mutation_Info))*100
      if(is.na(Hyper_GI)){
        Mutated_Sampl_Unniverse =length(unique(Current_Mut_Data$Tumor_Sample_Barcode[which(Current_Mut_Data$Hugo_Symbol%in%GI)]))
        Mutation_Rate=(Mutated_Sampl_Unniverse/length(unique(Current_Mut_Data$Tumor_Sample_Barcode)))*100
      }else{
        Mutated_Sampl_Unniverse =length(unique(Current_Mut_Data$Tumor_Sample_Barcode[which(Current_Mut_Data$Hugo_Symbol%in%Hyper_GI)]))
        Mutation_Rate=(Mutated_Sampl_Unniverse/length(unique(Current_Mut_Data$Tumor_Sample_Barcode)))*100
      }
      ## RNA of the Current Target Gene 
      source("CVN_Meth_Par_Cal_Rnd.R")
      Final_Data<-CVN_Meth_Par_Cal_Rnd(TG,RNA_GI,Mutation_Info,GI,Hyper_GI,TG_Par_Info)
      MutIndex=which(Final_Data$Miss_Com_Indx_U==2)
      
      if (length(MutIndex)>1){
        ############################################
        ### Model NonSyn Vs.Wild-Type
        source("SYS_Model.R")
        PostMutation2<-SYS_Model(Final_Data$CNV_Meth_Par_U,Final_Data$RNA_GI_U,Final_Data$RNA_TG_U,Final_Data$PG_Weights,Final_Data$Miss_Com_Indx_U,Final_Data$weightsGI)
        
        ################################# Control
        ## Calculating ratio of mutated to Non-Mutated samples
        Maxfold=2##20
        Result1=array(list(),dim=Maxfold)
        # for(Foldidx in 1:Maxfold){
        Result1 <-foreach(Foldidx=1:Maxfold)%dopar%{
          Wild_Com_Indx_U=array(NA,dim=Final_Data$Wild_type_Length)
          Group2_IDX=sample(Final_Data$Wild_type_Length, Final_Data$Numb_Grp_2, replace = FALSE)
          Wild_Com_Indx_U[Group2_IDX]=2 
          Wild_Com_Indx_U[is.na(Wild_Com_Indx_U)]=1
          source("SYS_Model.R")
          PostMutation_Wild<-SYS_Model(Final_Data$CNV_Meth_Par_Wild,Final_Data$RNA_GI_Wild,Final_Data$RNA_TG_Wild,Final_Data$PG_Weights,Wild_Com_Indx_U,Final_Data$weightsGI)
          Result1[[Foldidx]]<-PostMutation_Wild
        }
      }else{
        PostMutation2=NA
        Result1=NA
      }
      ###################################
      if(length(which(!is.na(match(Connection_Flag,Activa_Type))))){
        Conn_level=2
      }else{
        Conn_level=1
      }
      Old_TGG=TGsAll$Old_TGs[which(TGsAll$TGs%in%TG)]
      ###################################
      Analysis_TG<-list(
        Hyper_GI=Hyper_GI,
        GI=GI,
        Old_TG=Old_TGG,
        TG=TG,
        PG=as.character(Final_Data$PG),
        Res_NM_M=PostMutation2,
        Res_Cont=Result1,
        NonSyn_Label=Final_Data$Non_Syn_length,
        Wild_Label=Final_Data$Wild_type_Length,
        All_label=Final_Data$All_Mut_Indx,
        Conn_level=Conn_level,
        Mut_Rate_GI_Universe=Mutation_Rate,
        Mut_Rate_GI_Run=Mutation_Rate_RUN,
        Mutated_Sampl_Unniverse=Mutated_Sampl_Unniverse
      )
      Result_TG[[TGIdx]]=Analysis_TG
    }
    return(Result_TG)
  }else{
    return(NULL)
  }
}