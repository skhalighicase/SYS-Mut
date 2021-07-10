## Last Modification 2/5/2020
## Find missence mutations

Fun_Mutation_Labels<-function (GI){
  
  MutData_Temp=Current_Mut_Data[which(Current_Mut_Data$Hugo_Symbol%in%GI),]
  NonSynonumous_Labels=c("Missense_Mutation","Frame_Shift_Del","Nonsense_Mutation","In_Frame_Del","In_Frame_Ins","Splice_Site","Frame_Shift_Ins","Translation_Start_Site","Nonstop_Mutation")
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
    Mutation_Info$Mut_Info[is.na(Mutation_Info$Mut_Info)]<-1
  }else{
    SamplesAll$Mut_Info=NA
    Mutation_Info<-SamplesAll
    Mutation_Info$Mut_Info[is.na(Mutation_Info$Mut_Info)]<-1
  }
  return(Mutation_Info)
}