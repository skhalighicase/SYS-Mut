## Last Modification 9/15/2020
#### This function find all the required data for the Model.

Normal_Ele_Rel<-function(TGs,RNA_GI,GI,Hyper_GI){
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
      TG<-TGs[TGIdx]  ## RNA of the Current Target Gene 
      TG_Par_Info=All_Info_TG_Par[which(!is.na(match(All_Info_TG_Par$TG,TG))),]
      ##########################
      ## We will update the list of the samples for control. We eliminate the samples that have mutation on all related genes. parents and TG and GI. 
      SamplesAll<-data.frame("Samples_Name"=rownames(as.data.frame(RNA_GI)))
      ## Mutation Info ## A function to read and Lable the Missense Mutations
      source("Mutation_Labels_Clean_Wildtype.R")
      Mutation_Info<-Mutation_Labels_Clean_Wildtype(GI,Hyper_GI,TG,SamplesAll,TG_Par_Info) 
      
      ############# Mutation Rate
      # Mutation_Rate_RUN=(length(which(Mutation_Info$Mut_Info==2))/nrow(Mutation_Info))*100
      # 
      # if(is.na(Hyper_GI)){
      #   Mutated_Sampl_Unniverse =length(unique(Current_Mut_Data$Tumor_Sample_Barcode[which(Current_Mut_Data$Hugo_Symbol%in%GI)]))
      #   Mutation_Rate=(Mutated_Sampl_Unniverse/length(unique(Current_Mut_Data$Tumor_Sample_Barcode)))*100
      # }else{
      #   Mutated_Sampl_Unniverse =length(unique(Current_Mut_Data$Tumor_Sample_Barcode[which(Current_Mut_Data$Hugo_Symbol%in%Hyper_GI)]))
      #   Mutation_Rate=(Mutated_Sampl_Unniverse/length(unique(Current_Mut_Data$Tumor_Sample_Barcode)))*100
      # }
      
      ############# Mutation Rate 
      All_Sampl_Universe=nrow(Mutation_Info)
      Mutated_Sampl_Unniverse=length(which(Mutation_Info$Mut_Info==2))
      Mutation_Rate=(Mutated_Sampl_Unniverse/All_Sampl_Universe)*100
      
      ##########################
      source("CVN_Meth_Par_Cal_Rel.R")
      Final_Data<-CVN_Meth_Par_Cal_Rel(TG,RNA_GI,Mutation_Info,GI,Hyper_GI,TG_Par_Info)
      MutIndex=which(Final_Data$Miss_Com_Indx_U==2)  ## Nonsyn_samples
      Mutation_Rate_RUN=(Final_Data$Non_Syn_length/Final_Data$All_Mut_Indx)*100
      if (length(MutIndex)>1){
        ############################################
        ### Model NonSyn Vs.Wild-Type
        source("SYS_Model.R")
        PostMutation2<-SYS_Model(Final_Data$CNV_Meth_Par_U,Final_Data$RNA_GI_U,Final_Data$RNA_TG_U,Final_Data$PG_Weights,Final_Data$Miss_Com_Indx_U,Final_Data$weightsGI)
        ################################# Control
        ## Calculating ratio of mutated to Non-Mutated samples
        Maxfold=2 ## Maximum Fold
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
      ###################################
      Analysis_TG<-list(
        Hyper_GI=Hyper_GI,
        GI=GI,
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
        Mutated_Sampl_Unniverse=Mutated_Sampl_Unniverse,
        All_label_Universe=All_Sampl_Universe
      )
      Result_TG[[TGIdx]]=Analysis_TG
    }
    return(Result_TG)
  }else{
    return(NULL)
  }
}