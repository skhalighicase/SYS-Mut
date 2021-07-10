## Last Modification 6/6/2020
## Find Parents of a real Target Gene

Parents_Calc<-function(GI,Hyper_GI,TG_Par_Info){
  ## finding the parents
  PGAll<-as.character(unlist(TG_Par_Info$Parent_Genes))
 
  if(is.na(Hyper_GI)){
    GIindx=which(!is.na(match(PGAll,GI))) 
    Temp_GI=GI
  }else{
    GIindx=which(!is.na(match(PGAll,Hyper_GI)))
    Temp_GI=Hyper_GI
  }
  
  if(length(GIindx)){
    PG<-PGAll[-which(!is.na(match(PGAll,Temp_GI)))] ## Check if we have this Gene in the network
  }else{
    PG<-PGAll
  }
  return(PG)
}
