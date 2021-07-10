######## Last Modification 5/31/2021 ------------###  
######## Finding the best genes for random replacement 

Gene_Replace<-function(OldTGs,GI,Hyper_GI){
  if(is.na(Hyper_GI)){
    All_TGs_GI=Final_GGI_Net$TGs_All[which(Final_GGI_Net$GI_All%in%GI)]
  }else{
    All_TGs_GI=Final_GGI_Net$TGs_All[which(Final_GGI_Net$GI_All%in%Hyper_GI)] ## find all real TGs of GI
  }
  ############# Take all the information of the curent TG
  TG_PG_Old=TG_Parents_Connections[which(TG_Parents_Connections$TG%in%OldTGs),]
  Numb_TGs=length(OldTGs)
  ####------------------ Updated list of TGs
  TG_Par_Conn_Updated=TG_Parents_Connections[-which(TG_Parents_Connections$TG%in%All_TGs_GI),]  ## Find the parents of real TGs
  TGs=array(NA,Numb_TGs)
  TGs_OLD_Replace <- data.frame(TGs = character(),          # Specify empty vectors in data.frame
                                Old_TGs = character(),
                                stringsAsFactors = FALSE)
  RdTGIdx=1
  while(RdTGIdx <= Numb_TGs){
    TG_PG_Numb<-TG_PG_Old$PG_Num[RdTGIdx]
    #Parents_Curr_OldTG<-as.character(unlist(TG_PG_Old$Parent_Genes[RdTGIdx])) ###?? its better to select all parents of the current TGs not just current Gene
    Parents_Curr_OldTG<-unique(as.character(unlist(TG_PG_Old$Parent_Genes)))
    Indx_Parent_Old=which(TG_Par_Conn_Updated$TG%in%Parents_Curr_OldTG)
    if(length(Indx_Parent_Old)){
      TG_Par_Conn_Updated_Curr=TG_Par_Conn_Updated[-Indx_Parent_Old,]  ## Find the parents of real TGs
    }else{
      TG_Par_Conn_Updated_Curr=TG_Par_Conn_Updated  ## Find the parents of real TGs
    }
    ######
    Potential_replace_Genes=Closest(TG_Par_Conn_Updated_Curr$PG_Num,TG_PG_Numb)
    if(length(Potential_replace_Genes)>=20){
      Poten_TGsIdx=which(TG_Par_Conn_Updated_Curr$PG_Num%in%Potential_replace_Genes)
    }else if(length(Potential_replace_Genes)){
      Poten_TGsIdx=which(TG_Par_Conn_Updated_Curr$PG_Num%in%Potential_replace_Genes)
      Max_IDX=length(TG_Par_Conn_Updated_Curr$TG)
      ####---------------------------
      while(length(Poten_TGsIdx)<20){
        NewIDX=head(Poten_TGsIdx,n=1)-1
        if(length(NewIDX)&(NewIDX>0)){
          Poten_TGsIdx=c(NewIDX,Poten_TGsIdx)
        }else{
          NewIDX=tail(Poten_TGsIdx,n=1)+1
          if((NewIDX< Max_IDX)){
            Poten_TGsIdx=c(Poten_TGsIdx,NewIDX)
          }
        }
        NewIDX=tail(Poten_TGsIdx,n=1)+1
        if((NewIDX< Max_IDX)){
          Poten_TGsIdx=c(Poten_TGsIdx,NewIDX)
        }else{
          NewIDX=head(Poten_TGsIdx, n=1)-1
          if(length(NewIDX)&(NewIDX>0)){
            Poten_TGsIdx=c(NewIDX,Poten_TGsIdx)
          }
        }
      }
    }
    Candid_replac_Gen=as.character(sample(TG_Par_Conn_Updated_Curr$TG[Poten_TGsIdx],1,replace=FALSE))
    
    ####-------------- If we have seen this candidate gene do not use it, else add it to the list
    Seen_TGs=which(Candid_replac_Gen%in%Flags_Seen)
    if(!length(Seen_TGs)){
      Flags_Seen<<-c(Flags_Seen,Candid_replac_Gen)
      TGs[RdTGIdx]<-Candid_replac_Gen
      TGs_OLD_Replace[RdTGIdx, ] <- c(as.character(Candid_replac_Gen),as.character(TG_PG_Old$TG[RdTGIdx]))
      RdTGIdx=RdTGIdx+1
    }
    TG_Par_Conn_Updated=TG_Par_Conn_Updated[-which(TG_Par_Conn_Updated$TG==Candid_replac_Gen),]
  }
  return(TGs_OLD_Replace)
}