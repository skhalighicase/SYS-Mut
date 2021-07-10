Calc_Results_Measure<-function(){
  setpath="../Results"
  setwd(setpath)
  k=0
  datalist = list()
  files <- list.files(pattern = "\\.RData$")
  
  for (Fil_Idx in 1:length(files)){
    filename=files[Fil_Idx]
    load(filename)
    
    for (ResIndex in 1:length(RES)){
      Current_RES=RES[[ResIndex]]
      if (length(Current_RES)){
        for(i in 1:length(Current_RES)){
          Curr_GI=Current_RES[[i]]
          if(length(Curr_GI)){
            for(J in 1:length(Curr_GI)){
              Curr_Conection=Curr_GI[[J]]
              Data_Analysis=do.call(cbind.data.frame, Curr_Conection[1:4])
              Parent_Genes=length(unlist(Curr_Conection[5]))
              Data_Analysis=cbind(Data_Analysis,as.data.frame(Parent_Genes))
              samples_Wild_NonSyn=do.call(cbind.data.frame, Curr_Conection[8:14])
              Data_Analysis=cbind(Data_Analysis,samples_Wild_NonSyn)
              
              if (is.na(Curr_Conection[[6]])){
                Measures=list()
                Measures$Bht_Nm_M=0
                Measures$WPval=0
                Measures$Bht_C_Mean=0
                Measures$Bht_C_Var=0
              }else{
                ### Measurments Calculation
                source("../Template/Bhatt_NonSyn_Wild.R")
                Measures<-Bhatt_NonSyn_Wild(Curr_Conection[[6]],Curr_Conection[[7]])
              }
              Measur_Analysis=do.call(cbind.data.frame, Measures)
              Data_Analysis=cbind(Data_Analysis,Measur_Analysis)
              ###
              k=k+1
              datalist[[k]] <- Data_Analysis
              Data_Analysis=NULL
              Measur_Analysis=NULL
            }
            Curr_GI=NULL
          }
        }
      }
    }
  }
  
  Results_Impacts = do.call(rbind, datalist)
  write.table(Results_Impacts, "Results_NonSyn_Wild.txt", sep="\t",row.names = FALSE)
  write.table(Results_Impacts, "../Analysis/Results_NonSyn_Wild.txt", sep="\t",row.names = FALSE)  ### Write a copy of this file in Analysis folder
}