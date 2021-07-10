## Last Modification 9/12/2019

Normal_Ele<-function(TGs,GI,Hyper_GI){
  ## Initialization
  library(fpc)
  Seen_TGs=which(!is.na(match(TGs,Flags_Seen_GI)))
  if(length(Seen_TGs)){
    TGs=TGs[-Seen_TGs]
    if (length(TGs)){ 
      Flags_Seen_GI<<-c(Flags_Seen_GI,TGs)
    }
  }else{
    Flags_Seen_GI<<-c(Flags_Seen_GI,TGs)
  }
}

