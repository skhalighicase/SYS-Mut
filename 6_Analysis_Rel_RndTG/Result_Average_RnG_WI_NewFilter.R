
Result_Average_RnG_WI_NewFilter<-function(){
  ############ Iterations of random replacement of the network
  Max_Iteration=10
  #############
  Res_RelTG=read.delim(Params$filename_Rel,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Mutated_GI <-as.character(Res_RelTG$GI_ALL[Res_RelTG$Perce_Impact_Final>0])
  
  ############
  Res_RndTG=read.delim(Params$filename_Rand_Summ,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Res_temp=Res_RndTG[which(Res_RndTG$GI_ALL%in%Mutated_GI),]
  
    for(Fil_Idx in 2:Max_Iteration){
    filename=paste0("../5_SYSModel_Rnd_TG/RUNS/RUN_",Fil_Idx,"/Analysis/Results_NonSyn_Wild_Summry_GI_RUN_",Fil_Idx,".txt")
    Res_Temp_RndTG=read.delim(filename,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
    Result=Res_Temp_RndTG[which(Res_Temp_RndTG$GI_ALL%in%Res_temp$GI_ALL),]
    Res_temp$Perce_Impact_Final=Res_temp$Perce_Impact_Final+Result$Perce_Impact_Final
    Res_temp$Perce_Impact_Dir=Res_temp$Perce_Impact_Dir+Result$Perce_Impact_Dir
    Res_temp$Perce_Impact_InDir=Res_temp$Perce_Impact_InDir+Result$Perce_Impact_InDir
    Res_temp$Perce_Dir_High=Res_temp$Perce_Dir_High+Result$Perce_Dir_High
    Res_temp$Perce_Dir_Low=Res_temp$Perce_Dir_Low+Result$Perce_Dir_Low
    Res_temp$Perce_InDir_High=Res_temp$Perce_InDir_High+Result$Perce_InDir_High
    Res_temp$Perce_InDir_Low=Res_temp$Perce_InDir_Low+Result$Perce_InDir_Low
    Res_temp$Percent_All_High=Res_temp$Percent_All_High+Result$Percent_All_High
    Res_temp$Percent_All_Low=Res_temp$Percent_All_Low+Result$Percent_All_Low
    Res_temp$TG_All2_WO=Res_temp$TG_All2_WO+Result$TG_All2_WO
    Res_temp$TG_Dir_WI=Res_temp$TG_Dir_WI+Result$TG_Dir_WI
    Res_temp$TG_Dir_WO=Res_temp$TG_Dir_WO+Result$TG_Dir_WO
    Res_temp$TG_Dir_WH=Res_temp$TG_Dir_WH+Result$TG_Dir_WH
    Res_temp$TG_Dir_WL=Res_temp$TG_Dir_WL+Result$TG_Dir_WL
    Res_temp$TG_InDir_WI=Res_temp$TG_InDir_WI+Result$TG_InDir_WI
    Res_temp$TG_InDir_WO=Res_temp$TG_InDir_WO+Result$TG_InDir_WO
    Res_temp$TG_InDir_WH=Res_temp$TG_InDir_WH+Result$TG_InDir_WH
    Res_temp$TG_InDir_WL=Res_temp$TG_InDir_WL+Result$TG_InDir_WL
    Res_temp$TG_All_WI=Res_temp$TG_All_WI+Result$TG_All_WI
    Res_temp$TG_All_WO=Res_temp$TG_All_WO+Result$TG_All_WO
    Res_temp$TG_All_WH=Res_temp$TG_All_WH+Result$TG_All_WH
    Res_temp$TG_All_WL=Res_temp$TG_All_WL+Result$TG_All_WL
    Res_temp$GI_TG_total=Res_temp$GI_TG_total+Result$GI_TG_total
  }
  
  #### Average Calculation
  Res_temp$Perce_Impact_Final=Res_temp$Perce_Impact_Final/Max_Iteration
  Res_temp$Perce_Impact_Dir=Res_temp$Perce_Impact_Dir/Max_Iteration
  Res_temp$Perce_Impact_InDir=  Res_temp$Perce_Impact_InDir/Max_Iteration
  Res_temp$Perce_Dir_High=Res_temp$Perce_Dir_High/Max_Iteration
  Res_temp$Perce_Dir_Low=Res_temp$Perce_Dir_Low/Max_Iteration
  Res_temp$Perce_InDir_High=Res_temp$Perce_InDir_High/Max_Iteration
  Res_temp$Perce_InDir_Low=Res_temp$Perce_InDir_Low/Max_Iteration
  Res_temp$Percent_All_High=Res_temp$Percent_All_High/Max_Iteration
  Res_temp$Percent_All_Low=Res_temp$Percent_All_Low/Max_Iteration
  Res_temp$TG_All2_WO=Res_temp$TG_All2_WO/Max_Iteration
  Res_temp$TG_Dir_WI=Res_temp$TG_Dir_WI/Max_Iteration
  Res_temp$TG_Dir_WO=Res_temp$TG_Dir_WO/Max_Iteration
  Res_temp$TG_Dir_WH=Res_temp$TG_Dir_WH/Max_Iteration
  Res_temp$TG_Dir_WL=Res_temp$TG_Dir_WL/Max_Iteration
  Res_temp$TG_InDir_WI=Res_temp$TG_InDir_WI/Max_Iteration
  Res_temp$TG_InDir_WO=Res_temp$TG_InDir_WO/Max_Iteration
  Res_temp$TG_InDir_WH=Res_temp$TG_InDir_WH/Max_Iteration
  Res_temp$TG_InDir_WL=Res_temp$TG_InDir_WL/Max_Iteration
  Res_temp$TG_All_WI=Res_temp$TG_All_WI/Max_Iteration
  Res_temp$TG_All_WO=Res_temp$TG_All_WO/Max_Iteration
  Res_temp$TG_All_WH=Res_temp$TG_All_WH/Max_Iteration
  Res_temp$TG_All_WL=Res_temp$TG_All_WL/Max_Iteration
  Res_temp$GI_TG_total=Res_temp$GI_TG_total/Max_Iteration
  ##
  
  Output_filename=paste0("./Results/Avg_Perc_Impacted_RndTG_",CANCER_TYPE,".txt")
  write.table(Res_temp, Output_filename, sep="\t",row.names = FALSE)
}