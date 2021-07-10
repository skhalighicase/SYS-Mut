## Last Modification 2/5/2020
## Find Parents of a real Target Gene

Parents_Calc_Old<-function(TG,GI,Hyper_GI){
  ## Finding the parents
  IndxTGs<-which(!is.na(match(Final_GI_TG_PG$TG,TG)))## Check if we have this Gene in the network
  
  PGAll=as.character(unlist(Final_GI_TG_PG$Parent_Genes[IndxTGs]))
  
  if(is.na(Hyper_GI)){
    GIindx=which(!is.na(match(PGAll,GI))) 
    Temp_GI=GI
  }else{
    GIindx=which(!is.na(match(PGAll,Hyper_GI)))
    Temp_GI=Hyper_GI
  }
  if(length(GIindx)){
    PG_Curr<-PGAll[-which(!is.na(match(PGAll,Temp_GI)))]## Check if we have this Gene in the network
  }else{
    PG_Curr<-PGAll
  }
  
  if(length(PG_Curr)>100){  ## Maximum length of the parents that we can consider in the regression.
    PG<-intersect(rownames(AllGene_Info),PG_Curr) ## we have filtered some genes that they have small expression(Filtered applied in the main script) so here maybe we eliminate some of the PG_Currents
  }else{
    PG<-PG_Curr
  }
  return(PG)
}
