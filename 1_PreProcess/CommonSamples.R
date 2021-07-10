## This function finds the common samples.
CommonSamples<-function(Mutation_Pan){
  File_Name=paste0("./Data/CommonGeneData_",CANCER_TYPE,".RData")
  load(File_Name)
  load(Mutation_Pan)

  ## Find common Normal samples 
  ### Find duplicates of Methylation data
  if(length(Comm_RNAS0)){
    SampNames_Methyl0 <-colnames(Comm_Meth_Max0)
    SampNames_Methyl0_Cod=substr(SampNames_Methyl0,1,16)
    Dup_Metyl0=which(duplicated(SampNames_Methyl0_Cod))
    if(length(Dup_Metyl0)!=0){
      Comm_Meth_Max0 =Comm_Meth_Max0[,-Dup_Metyl0]
      Comm_Meth_Mean0=Comm_Meth_Mean0[,-Dup_Metyl0]
      Comm_Meth_Min0=Comm_Meth_Min0[,-Dup_Metyl0]
      SampNames_Methyl0 <-colnames(Comm_Meth_Max0)
      SampNames_Methyl0_Cod=substr(SampNames_Methyl0,1,16)
    }
    ### Find duplicates of CNV data
    SampNames_CNV0 <-colnames(Comm_CNV0)
    SampNames_CNV0_Cod=substr(SampNames_CNV0,1,16)
    Dup_CNV0=which(duplicated(SampNames_CNV0_Cod))
    if(length(Dup_CNV0)!=0){
      Comm_CNV0=Comm_CNV0[,-Dup_CNV0]
      SampNames_CNV0 <-colnames(Comm_CNV0)
      SampNames_CNV0_Cod=substr(SampNames_CNV0,1,16)
    }
    ### Find duplicates of RNA data
    SampNames_RNASeq0<-colnames(Comm_RNAS0)
    SampNames_RNASeq0_Cod=substr(SampNames_RNASeq0,1,16)
    Dup_RNA0=which(duplicated(SampNames_RNASeq0_Cod))
    if(length(Dup_RNA0!=0)){
      Comm_RNAS0=Comm_RNAS0[,-Dup_RNA0]
      SampNames_RNASeq0<-colnames(Comm_RNAS0)
      SampNames_RNASeq0_Cod=substr(SampNames_RNASeq0,1,16)
    }
    
    Com_Samp0<-intersect(intersect(SampNames_CNV0_Cod,SampNames_RNASeq0_Cod),SampNames_Methyl0_Cod)
    
    ## CNV0
    Com_CNV_RNA0_Idx<-which(SampNames_CNV0_Cod%in%Com_Samp0)
    Comm_Sam_CNV0=Comm_CNV0[,Com_CNV_RNA0_Idx]
    
        ## Methy0
    Com_Methyl_RNA0_Idx<-which(SampNames_Methyl0_Cod%in%Com_Samp0)
    
    Comm_Sam_Meth_Max0=Comm_Meth_Max0[,Com_Methyl_RNA0_Idx]
    Comm_Sam_Meth_Mean0=Comm_Meth_Mean0[,Com_Methyl_RNA0_Idx]
    Comm_Sam_Meth_Min0=Comm_Meth_Min0[,Com_Methyl_RNA0_Idx]
    
    ## RNASeq0
    Com_RNA0_Idx<-which(SampNames_RNASeq0_Cod%in%Com_Samp0)
    Comm_Sam_RNAS0=Comm_RNAS0[,Com_RNA0_Idx]
  }else{
    Comm_Sam_RNAS0=matrix(NA,0,0)
    Comm_Sam_Meth_Min0=matrix(NA,0,0)
    Comm_Sam_Meth_Mean0=matrix(NA,0,0)
    Comm_Sam_Meth_Max0=matrix(NA,0,0)
    Comm_Sam_CNV0=matrix(NA,0,0)
  }
  
  ########################
  ## Find common Cancer samples
  Methyl_Samples=data.frame(FullName_Methyl_Samples=colnames(Comm_Meth_Max1))
  Methyl_Samples$Methyl_Sample_Name=substr(Methyl_Samples$FullName_Methyl_Samples,1,16)
  Methyl_Samples$Methyl_Patient_Name=substr(Methyl_Samples$FullName_Methyl_Samples,1,12)
  Methyl_Samples = Methyl_Samples[!duplicated(Methyl_Samples$Methyl_Sample_Name),]
  
  
  CNV_Samples=data.frame(FullName_CNV_Samples=colnames(Comm_CNV1))
  CNV_Samples$CNV_Sample_Name=substr(CNV_Samples$FullName_CNV_Samples,1,16)
  CNV_Samples$CNV_Patient_Name=substr(CNV_Samples$FullName_CNV_Samples,1,12)
  CNV_Samples = CNV_Samples[!duplicated(CNV_Samples$CNV_Sample_Name),]
  
  
  RNA_Samples=data.frame(FullName_RNA_Samples=colnames(Comm_RNAS1))
  RNA_Samples$RNA_Sample_Name=substr(RNA_Samples$FullName_RNA_Samples,1,16)
  RNA_Samples$RNA_Patient_Name=substr(RNA_Samples$FullName_RNA_Samples,1,12)
  RNA_Samples = RNA_Samples[!duplicated(RNA_Samples$RNA_Sample_Name),]
  
  ###### Samples and Patients of MAF file
  Mutation_Samples=data.frame(FullName_Sample_MAF=unique(Mutt_Current_Tissue$Tumor_Sample_Barcode))
  Mutation_Samples$MAF_Sample_Name=unique(substr(Mutation_Samples$FullName_Sample_MAF,1,16))
  Mutation_Samples$MAF_Patient_Name=unique(substr(Mutation_Samples$FullName_Sample_MAF,1,12))
  Mutation_Samples = Mutation_Samples[!duplicated(Mutation_Samples$MAF_Sample_Name),]

  Overlapped_Samples_RNA=merge(RNA_Samples,Mutation_Samples,by.x="RNA_Patient_Name",by.y ="MAF_Patient_Name",all.x = T,all.y = T)
  Com_Samp_RNA_Mut_1=Overlapped_Samples_RNA[complete.cases(Overlapped_Samples_RNA), ]
  Com_Samp_RNA_Mut_1 = Com_Samp_RNA_Mut_1[!duplicated(Com_Samp_RNA_Mut_1$RNA_Patient_Name),]
  
  Overlapped_Samples_CNV=merge(CNV_Samples,Mutation_Samples,by.x="CNV_Patient_Name",by.y ="MAF_Patient_Name",all.x = T,all.y = T)
  Com_Samp_CNV_Unierse_1=Overlapped_Samples_CNV[complete.cases(Overlapped_Samples_CNV), ]
  Com_Samp_CNV_Unierse_1 = Com_Samp_CNV_Unierse_1[!duplicated(Com_Samp_CNV_Unierse_1$CNV_Patient_Name),]
  
  Overlapped_Samples_Methy=merge(Methyl_Samples,Mutation_Samples,by.x="Methyl_Patient_Name",by.y ="MAF_Patient_Name",all.x = T,all.y = T)
  Com_Samp_Methyl_Unierse_1=Overlapped_Samples_Methy[complete.cases(Overlapped_Samples_Methy), ]
  Com_Samp_Methyl_Unierse_1 = Com_Samp_Methyl_Unierse_1[!duplicated(Com_Samp_Methyl_Unierse_1$Methyl_Patient_Name),]
  
  #######
  Comm_Sam_CNV1=Comm_CNV1[,which(colnames(Comm_CNV1)%in%as.character(Com_Samp_CNV_Unierse_1$FullName_CNV_Samples))] 
  Comm_Sam_RNAS1=Comm_RNAS1[,which(colnames(Comm_RNAS1)%in%as.character(Com_Samp_RNA_Mut_1$FullName_RNA_Samples))] 
  Comm_Sam_Meth_Max1=Comm_Meth_Max1[,which(colnames(Comm_Meth_Max1)%in%as.character(Com_Samp_Methyl_Unierse_1$FullName_Methyl_Samples))] 
  Comm_Sam_Meth_Mean1=Comm_Meth_Mean1[,which(colnames(Comm_Meth_Mean1)%in%as.character(Com_Samp_Methyl_Unierse_1$FullName_Methyl_Samples))] 
  Comm_Sam_Meth_Min1=Comm_Meth_Min1[,which(colnames(Comm_Meth_Min1)%in%as.character(Com_Samp_Methyl_Unierse_1$FullName_Methyl_Samples))] 
  
  ## Mutation
  Current_Mut_Data=Mutt_Current_Tissue[which(Mutt_Current_Tissue$Tumor_Sample_Barcode%in%as.character(Com_Samp_RNA_Mut_1$FullName_Sample_MAF)),] 
  
  File_Name_Mut=paste0("./Data/Updated_Mut_",CANCER_TYPE,".RData")
  save(Current_Mut_Data,file=File_Name_Mut)
  
  File_Name_Save=paste0("./Data/CommonGenesamples_",CANCER_TYPE,".RData")
  save(Comm_Sam_RNAS0,Comm_Sam_Meth_Min0,Comm_Sam_Meth_Mean0,Comm_Sam_Meth_Max0,Comm_Sam_CNV0,Comm_Sam_RNAS1,Comm_Sam_Meth_Min1,Comm_Sam_Meth_Mean1,Comm_Sam_Meth_Max1,Comm_Sam_CNV1,file=File_Name_Save)
}