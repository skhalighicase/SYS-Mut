## Last Modification 2/5/2020
CVN_Meth_Parent_Calcu<-function(TG,RNA_GI,GI,Hyper_GI){
  RNA_TG<-RNA[TG,]
  ## CNV data 
  CNV_Gen_List=rownames(CNV_Methyl$CNV)
  if(TG%in%CNV_Gen_List){
    CNV_TG<-CNV_Methyl$CNV[TG,] ## Trasform 
    Common_CNV_RNA=intersect(names(CNV_TG),names(RNA_TG))
  }else{
    Common_CNV_RNA=names(RNA_TG)
    CNV_TG<-array(NA,length(Common_CNV_RNA))
    names(CNV_TG)=Common_CNV_RNA
  }
  ## Methyl Data
  Methyl_Gen_List=rownames(CNV_Methyl$Meth_Mean)
  if (TG%in%Methyl_Gen_List){
    Meth_TG <-CNV_Methyl$Meth_Mean[TG,]
    Common_Meth_RNA=intersect(names(Meth_TG),names(RNA_TG))
  }else{
    Common_Meth_RNA=names(RNA_TG)
    Meth_TG<-array(NA,length(Common_Meth_RNA))
    names(Meth_TG)=Common_Meth_RNA
  }
  Universe_Samples=intersect(Common_CNV_RNA,Common_Meth_RNA)
  
  ########
  CNV_TG=CNV_TG[which(!is.na(match(names(CNV_TG),Universe_Samples)))]
  CNV_TG=CNV_TG[Universe_Samples]
  Meth_TG=Meth_TG[which(!is.na(match(names(Meth_TG),Universe_Samples)))]
  Meth_TG=Meth_TG[Universe_Samples]
  
  RNA_TG=RNA_TG[which(!is.na(match(names(RNA_TG),Universe_Samples)))]
  RNA_TG=RNA_TG[Universe_Samples]
  ## Update RNA
  if (length(Universe_Samples)!=dim(RNA)[2]){
    Updated_RNA_Sampl=RNA[,which(!is.na(match(colnames(RNA),Universe_Samples)))]
  }else{
    Updated_RNA_Sampl=RNA
  }
  Updated_RNA_Sampl=Updated_RNA_Sampl[,Universe_Samples]
  
  ###################################
  ## Parent Calculation for Real Target Genes
  source("Parents_Calc_Old.R")
  PG<-Parents_Calc_Old(TG,GI,Hyper_GI)
  ## RNA of the parent Genes
  if(length(PG)){
    RNA_PG<-Updated_RNA_Sampl[PG,]
  }else{
    RNA_PG<-array(NA,dim(Updated_RNA_Sampl)[2])
    names(RNA_PG)=colnames(Updated_RNA_Sampl)
  }
  #############################
  ## Appending the CNV and Parents and Methyl data
  CNV_Meth_Par<-rbind(CNV_TG,Meth_TG,RNA_PG)   ## adding Cis regulators and parents
  
  #######################
  if (is.matrix(CNV_Meth_Par)){
    CNV_Meth_Par=CNV_Meth_Par[rowSums(is.na(CNV_Meth_Par)) != ncol(CNV_Meth_Par),]
    if(is.matrix(CNV_Meth_Par)){
      if (!(dim(CNV_Meth_Par)[1])){
        CNV_Meth_Par=array(0,length(RNA_TG))
      }
    }
  }
  
  if (is.matrix(CNV_Meth_Par)){
    PG_Weights<-array(0,dim(CNV_Meth_Par)[1])
  }else{
    PG_Weights<-array(0,1)
    CNV_Meth_Par=array(0,length(RNA_TG))
  }
  weightsGI=0
  Mutation_Info=Mutation_Info[which(Mutation_Info$Samples_Name%in%Universe_Samples),]
  
  ## Find the Noise and Synon genes
  Syn_Noise_Indx=which((Mutation_Info$Mut_Info==3)|(Mutation_Info$Mut_Info==4))
  
  ## Updating data based on the type of mutation
  if(length(Syn_Noise_Indx)){
    if(is.matrix(CNV_Meth_Par)){
      CNV_Meth_Par_U=CNV_Meth_Par[,-Syn_Noise_Indx]
    }else{
      CNV_Meth_Par_U=CNV_Meth_Par[-Syn_Noise_Indx]
    }
    RNA_GI=RNA_GI[which(!is.na(match(names(RNA_GI),Universe_Samples)))]
    RNA_GI=RNA_GI[Universe_Samples]
    
    RNA_GI_U=RNA_GI[-Syn_Noise_Indx]
    RNA_TG_U=RNA_TG[-Syn_Noise_Indx]
    Miss_Com_Indx_U=Mutation_Info$Mut_Info[-Syn_Noise_Indx]
  }else{
    CNV_Meth_Par_U=CNV_Meth_Par
    RNA_GI_U=RNA_GI
    RNA_TG_U=RNA_TG
    Miss_Com_Indx_U=Mutation_Info$Mut_Info
  }
  
  Non_Syn_length=length(which(Miss_Com_Indx_U==2))
  WildType_Idx=which(Miss_Com_Indx_U==1)
  Wild_type_Length=length(WildType_Idx)
  TRN_All_Mu_Wild=Wild_type_Length+Non_Syn_length
  
  ############################# 
  NonSyn_Wild_Ratio=Non_Syn_length/TRN_All_Mu_Wild
  Numb_Grp_2=floor(Wild_type_Length*NonSyn_Wild_Ratio)
  
  ## Updating data based on the type of mutation for control group
  if (length(WildType_Idx)){
    if(is.matrix(CNV_Meth_Par_U)){
      CNV_Meth_Par_Wild=CNV_Meth_Par_U[,WildType_Idx]
    }else{
      CNV_Meth_Par_Wild=CNV_Meth_Par_U[WildType_Idx]
    }
    RNA_GI_Wild=RNA_GI_U[WildType_Idx]
    RNA_TG_Wild=RNA_TG_U[WildType_Idx]
  }else{
    CNV_Meth_Par_Wild=NA
    RNA_GI_Wild=NA
    RNA_TG_Wild=NA
  }
  ###############################
  Miss_Com_Indx_U[which(Miss_Com_Indx_U==2)]=1
  
  Final_Data<-list(
    RNA_GI_U=RNA_GI_U,
    RNA_TG_U=RNA_TG_U,
    Miss_Com_Indx_U=Miss_Com_Indx_U,
    CNV_Meth_Par_U=CNV_Meth_Par_U,
    weightsGI=weightsGI,
    PG_Weights=PG_Weights,
    PG=PG,
    RNA_GI_Wild=RNA_GI_Wild,
    RNA_TG_Wild=RNA_TG_Wild,
    CNV_Meth_Par_Wild=CNV_Meth_Par_Wild,
    Numb_Grp_2=Numb_Grp_2,
    Wild_type_Length=Wild_type_Length,
    TRN_All_Mu_Wild=TRN_All_Mu_Wild,
    Non_Syn_length=Non_Syn_length
  )
  return(Final_Data)
}