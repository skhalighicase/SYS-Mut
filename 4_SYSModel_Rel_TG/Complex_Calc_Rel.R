## Last Modification 6/6/2020
## Find the connections type Complex

Complex_Calc_Rel<-function(Complex,RNA_GI,Hyper_GI){
  Complex_Eles=as.vector(Complex$Complex_Elements)
  for(Comp_ElIdx in 1:length(Complex_Eles)){ #l
    GI3<-Complex_Eles[Comp_ElIdx]
    Connection_Flag<<-c(Connection_Flag,as.character(Complex$Con_Type[Comp_ElIdx]))
    source("TGs_Find.R")
    TGs_GI3<-TGs_Find(GI3)
    if (length(TGs_GI3)){
      for (TG_Ele_Indx3 in 1:length(TGs_GI3)){## It may have two groups of target Genes, -t and others 1for -t 2 for others
        ###############################################
        switch(TG_Ele_Indx3, 
               "1"={
                 if (!is.null(TGs_GI3[[TG_Ele_Indx3]])){
                   TGs3=as.vector(TGs_GI3[[1]]$Gene_Elements)
                   source("Normal_Ele_Rel.R")
                   Result[[Res_IDX]]<<-Normal_Ele_Rel(TGs3,RNA_GI,GI3,Hyper_GI)
                   Res_IDX<<-Res_IDX+1
                 }
               },
               "2"={ ## Running the model in all Complex cases
                 if (!is.null(TGs_GI3[[TG_Ele_Indx3]])){
                   Complex_Eles3=as.vector(TGs_GI3[[TG_Ele_Indx3]]$Complex_Elements)
                   Seen_Complex=which(!is.na(match(Complex_Eles3,Flag_Complex_Seen)))
                   if(length(Seen_Complex)){
                     Complex_Eles3=Complex_Eles3[-Seen_Complex]
                     if (length(Complex_Eles3)){
                       Flag_Complex_Seen<<-c(Flag_Complex_Seen,Complex_Eles3)
                       source("Complex_Calc_Rel.R")
                       Complex_Calc_Rel(TGs_GI3[[TG_Ele_Indx3]],RNA_GI,Hyper_GI)
                     }
                   }else{
                     Flag_Complex_Seen<<-c(Flag_Complex_Seen,Complex_Eles3)
                     source("Complex_Calc_Rel.R")
                     Complex_Calc_Rel(TGs_GI3[[TG_Ele_Indx3]],RNA_GI,Hyper_GI)
                   }
                 }
               }
        )
      }
    }
    Connection_Flag<<-Connection_Flag[-length(Connection_Flag)]
  }
}