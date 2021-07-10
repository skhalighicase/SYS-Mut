## Last Modification 2/5/2020
Find_AllConnect_Gen<-function(Final_GGI_Net,Final_AllGI,Final_AllTG){
  ################## Load All possible connections between ALL GIs and TGs
  
  TG_Parents=array(list(),dim=length(Final_AllTG))
  
  for(TG_Idx in 1:length(Final_AllTG)){
    TG=Final_AllTG[TG_Idx]
    PGs_All=Final_GGI_Net$GI_All[which(Final_GGI_Net$TGs_All==TG)]
    PGs=unique(PGs_All)
    PG_Num=length(PGs)
    Analysis_TG<-list(
      TG=TG,
      PG_Num=PG_Num,
      PG=PGs
    )
    TG_Parents[[TG_Idx]]=Analysis_TG
    Analysis_TG=NULL
    PGs=NULL
    PGs_All=NULL
  }
  
  ######################
  datalist = list()
  k=0
  for (Fil_Idx in 1:length(TG_Parents)){
    TG_Info=TG_Parents[[Fil_Idx]]
    Data_Analysis=do.call(cbind.data.frame, TG_Info[1:2])
    Data_Analysis$Parent_Genes=TG_Info[3]
    k=k+1
    datalist[[k]] <- Data_Analysis
    Data_Analysis=NULL
  }
  Result = do.call(rbind,datalist)
  TG_Parents_Connections <- Result[order(Result$PG_Num),] 
  
  save(TG_Parents_Connections,file="./Results/TG_PGs_Connec.RData")
  
  ##################### Find the Gene Connections 2/5/2020 (Find_ALL_GI_TG_PG)
  
  for (Indx in 1:length(TG_Parents_Connections$TG)) {
    TG_Parents_Connections$GI[Indx]=TG_Parents_Connections$Parent_Genes[Indx]$PG[1]
  }
  
  Final_GI_TG_PG=TG_Parents_Connections[,c(4,1,2,3)]
  Saved_File_name=paste0("./Results/Final_GI_TG_PG.RData") ## this file save all TGs of each GI.
  save(Final_GI_TG_PG,file=Saved_File_name)
}

