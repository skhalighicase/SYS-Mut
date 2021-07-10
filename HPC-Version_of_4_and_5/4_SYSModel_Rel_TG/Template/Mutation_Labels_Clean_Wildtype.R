## Last Modification 6/30/2020
## Find missence mutations

Mutation_Labels_Clean_Wildtype<-function (GI,Hyper_GI,TG,SamplesAll,TG_Par_Info){
  if(is.na(Hyper_GI)){
    MutData_Temp=Current_Mut_Data[which(Current_Mut_Data$Hugo_Symbol%in%GI),]
  }else{
    MutData_Temp=Current_Mut_Data[which(Current_Mut_Data$Hugo_Symbol%in%Hyper_GI),]
  }
  
  NonSynonumous_Labels=c("Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","Splice_Site","Frame_Shift_Ins","Translation_Start_Site","In_Frame_Ins","Nonstop_Mutation")
  Synonumous_Labes=c("Silent","3'UTR","RNA","Intron","3'Flank","5'UTR","5'Flank")
  
  if(nrow(MutData_Temp)){
    Mut_Sample_Info=data.frame("Samples_Name"=unique(MutData_Temp$Tumor_Sample_Barcode),"Mut_Info"=NA)
    
    NonSyn_all=MutData_Temp[which(MutData_Temp$Variant_Classification%in%NonSynonumous_Labels),]
    NonSyn_Temp=NonSyn_all[which((NonSyn_all$t_depth>=30)&(NonSyn_all$n_depth>=30)&(NonSyn_all$t_freq>=0.1)&(NonSyn_all$n_freq<=0.05)),]
    Mut_Sample_Info$Mut_Info[which(Mut_Sample_Info$Samples_Name%in%unique(NonSyn_Temp$Tumor_Sample_Barcode))]=2   ###Non-Syn_Labels
    
    ### Synanoumus samples
    Syn_Samples=unique(MutData_Temp$Tumor_Sample_Barcode[which(MutData_Temp$Variant_Classification%in%Synonumous_Labes)])
    Mut_Sample_Info$Mut_Info[which(Mut_Sample_Info$Samples_Name%in%Syn_Samples)]=3  ##Syn_Labels
    Mut_Sample_Info$Mut_Info[is.na(Mut_Sample_Info$Mut_Info)]=4   ### Non_Syn_Noise
    Mut_Sample_Info$Samples_Name=substr(Mut_Sample_Info$Samples_Name,1,16)
    Mutation_Info<-merge(SamplesAll,Mut_Sample_Info,all.x=TRUE)
    Mutation_Info$Mut_Info[is.na(Mutation_Info$Mut_Info)]<-NA
  }else{
    SamplesAll$Mut_Info=NA
    Mutation_Info<-SamplesAll
    Mutation_Info$Mut_Info[is.na(Mutation_Info$Mut_Info)]<-NA
  }
  #### If you want to apply the clean Wild-type idea use this peace of code if not comment it
  PG=TG_Par_Info$Parent_Genes
  Connected_Genes=c(PG$P_Name,TG)
  MutData_Connec_Gen=Current_Mut_Data[which(Current_Mut_Data$Hugo_Symbol%in%Connected_Genes),]
  NonSyn_Connect=MutData_Connec_Gen[which(MutData_Connec_Gen$Variant_Classification%in%NonSynonumous_Labels),]
  NonSyn_Samples_Conn=unique(substr(NonSyn_Connect$Tumor_Sample_Barcode,1,16))
  Wild_type_GI=Mutation_Info[is.na(Mutation_Info$Mut_Info),]
  NonSyn_Connect=Wild_type_GI$Samples_Name[which(!is.na(match(Wild_type_GI$Samples_Name,NonSyn_Samples_Conn)))]
  Mutation_Info$Mut_Info[which(Mutation_Info$Samples_Name%in%NonSyn_Connect)]=3
  Mutation_Info$Mut_Info[is.na(Mutation_Info$Mut_Info)]=1
  Wild_ty_Numb=length(which(Mutation_Info$Mut_Info%in%1))
  
  if(Wild_ty_Numb==0){
    Sampl_Not_Clean=which(Mutation_Info$Samples_Name%in%NonSyn_Connect)
    Candid_replac_Samples=sample(Sampl_Not_Clean,3,replace=FALSE)  ## When we don't have enough sample for wild-type we randomly select samples that with the mutated parents.
    Mutation_Info$Mut_Info[Candid_replac_Samples]=1
  }
  return(Mutation_Info)
}