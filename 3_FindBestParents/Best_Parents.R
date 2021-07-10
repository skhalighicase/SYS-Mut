## Last Modification 2/5/2020
##  Merge the Best Parents 
Best_Parents<-function(){
  setpath="./Results"
  setwd(setpath)
  Max_Leng_Parent = 100  ### we have considered Maximum 100 confounder in the analysis 
  k = 0
  TG_Par_Info = list()
  TG_Par_List = list()
  files <- list.files(pattern = "\\.RData$")
  for (Fil_Idx in 1:length(files)) {
    filename = files[Fil_Idx]
    load(filename)
    for (ResIndex in 1:length(RES)) {
      Current_RES = RES[[ResIndex]]
      if (length(Current_RES)) {
        for (i in 1:length(Current_RES)) {
          Curr_GI = Current_RES[[i]]
          if (length(Curr_GI)) {
            for (J in 1:length(Curr_GI)) {
              Curr_Conection = Curr_GI[[J]]
              Parent_leng = length(unlist(Curr_Conection[4]))
              if (Parent_leng) {
                Parent_Names = Curr_Conection$PG
                Mu_Beta_Par = Curr_Conection$Res_NM_M$BetaInfo$Mu_Beta
                Diff_Length = length(Mu_Beta_Par) - length(Parent_Names)
                if (Diff_Length) {
                  Mu_Beta_Par = Mu_Beta_Par[-(1:Diff_Length)]
                }
                if (is.na(Curr_Conection$Hyper_GI)) {
                  Parent_Names = c(Curr_Conection$GI,Parent_Names)
                }else{
                  Parent_Names = c(Curr_Conection$Hyper_GI,Parent_Names)
                }
                Parent_leng = Parent_leng + 1
                Mu_Beta_Par = c(Curr_Conection$Res_NM_M$MutInfo$Mu_Phi2,Mu_Beta_Par)
                Parents = data.frame(Parent_Names, Mu_Beta_Par)
                if ((length(Parents$Parent_Names)) > Max_Leng_Parent) {
                  Parents <- Parents[order(-abs(Parents$Mu_Beta_Par)), ]
                  Parents = Parents[c(1:Max_Leng_Parent), ]
                  Remove_Idx = which(abs(Parents$Mu_Beta_Par) < 0.01)
                  if (length(Remove_Idx)) {
                    Parents = Parents[-Remove_Idx, ]
                  }
                }
              }
              else {
                if (is.na(Curr_Conection$Hyper_GI)) {
                  Parent_Names = Curr_Conection$GI
                }
                else {
                  Parent_Names = Curr_Conection$Hyper_GI
                }
                Parent_leng = Parent_leng + 1
                Parents = data.frame(Parent_Names)
              }
              TG_Par_Info = list(TG = as.character(Curr_Conection[3]), 
                                 Parent_leng = Parent_leng, P_Name = as.character(Parents$Parent_Names))
              k = k + 1
              TG_Par_List[[k]] <- TG_Par_Info
              Parents = NULL
              Parent_Names = NULL
              Mu_Beta_Par = NULL
              TG_Par_Info = NULL
            }
            Curr_GI = NULL
          }
        }
      }
    }
  }
  
  datalist = list()
  M = 0
  for (Fil_Idx in 1:length(TG_Par_List)) {
    TG_Info = TG_Par_List[[Fil_Idx]]
    Data_Analysis1 = do.call(cbind.data.frame, TG_Info[1:2])
    Data_Analysis1$Parent_Genes = TG_Info[3]
    M = M + 1
    datalist[[M]] <- Data_Analysis1
    Data_Analysis1 = NULL
  }
  All_Info_TG_Par = do.call(rbind, datalist)
  All_Info_TG_Par <- All_Info_TG_Par[order(All_Info_TG_Par$Parent_leng),]
  
  #### Assign the real parent names of each TG
  
  for (I in 1:length(All_Info_TG_Par$TG)){
    Curr_Gen=as.character(All_Info_TG_Par$TG[I])
    All_Info_TG_Par$Parent_leng[I]=Final_GI_TG_PG$PG_Num[which(!is.na(match(Final_GI_TG_PG$TG,Curr_Gen)))]
  }
  File_Name=paste0("Results_TG_Parents_",CANCER_TYPE,".RData")
  save(All_Info_TG_Par, file = File_Name)
}