### Find the Gene Connections 2/5/2020
graphics.off() 
rm(list=ls())

library(stringi)
library(stringr)
library(lattice)
library(gridExtra)
library(gapminder)
library(grid)
library(fpc)
library(DescTools)

################## Update and Setpath
CurrentDir=getwd()
setwd(CurrentDir)

source("../1_PreProcess/Init_Parameters.R")
Params<-Init_Parameters()

#--------------- 
CANCER_TYPE<<-Params$CANCER

################## Load Combined Data
FileName_1=paste0("../1_PreProcess/Data/Combin_Gen_Samp_RNA_",CANCER_TYPE,".RData")
load(FileName_1)
RNA_Data_GI=Nor_Comm_Data1$Nor_Comm_RNA1_GI
RNA_Data_AllGenes=Nor_Comm_Data1$Nor_Comm_RNA1
rm(Nor_Comm_Data1)

################## Load Influence Network
load("../1_PreProcess/Data/SYS_Net.RData")
Influ_net<<-Influ_net

##################  Define our Universe for GI and TG
Exist_TG<<-rownames(RNA_Data_AllGenes)
All_GIs=rownames(RNA_Data_GI)

RES=array(list(),dim=length(All_GIs))  ## Save the results
################## Connection Types, For all Genes that have Expression data We will do the analysis :: All_GIs
Trans_Types<<-c("-t","-t|","-t>")
Activa_Type<<-c("-a|","-a>")

##################  For All GIs
for(GIIdx in 1:length(All_GIs)){ #
  GI<-sample(All_GIs,1)
  ############## Find Possible Target Genes of GI (including Complex, TG and family)
  Connection_Flag<<-NULL
  source("TGs_Find.R")
  TGs_GI<-TGs_Find(GI)
  
  ############## If GI does not have any target gene goes to else
  if(length(TGs_GI)){ 
    Flags_Seen_GI<<-GI
    for(TG_Ele_Indx in 1:length(TGs_GI)){ 
      ######################
      switch(TG_Ele_Indx,
             "1"={ ### Running the model for all Normal Genes
               if (!is.null(TGs_GI[[TG_Ele_Indx]])){
                 Hyper_GI<-NA
                 TGs=as.vector(TGs_GI[[TG_Ele_Indx]]$Gene_Elements)
                 ##### Expression GI or Hyper_GI
                 source("Normal_Ele.R")
                 Normal_Ele(TGs,GI,Hyper_GI)
               }
             },
             "2"={ ### Running the model in all Complex cases
               if (!is.null(TGs_GI[[TG_Ele_Indx]])){
                 Hyper_GI<-GI
                 Flag_Complex_Seen<<-as.vector(TGs_GI[[TG_Ele_Indx]]$Complex_Elements)
                 source("Complex_Caclutation.R")
                 Complex_Caclutation(TGs_GI[[TG_Ele_Indx]],Hyper_GI)
               }
             }
      )
    }
    GI_All=Flags_Seen_GI[1]
    TGs_All=Flags_Seen_GI[-1]
    if(!length(TGs_All)){
      TGs_All=NA
    }
    RES[[GIIdx]]=data.frame(GI_All,TGs_All,stringsAsFactors = FALSE)
    GI_All=NULL
    TGs_All=NULL
    Flag_Complex_Seen=NULL
    Flags_Seen_GI=NULL
  }
  All_GIs<-All_GIs[-which(All_GIs==GI)]
}
GGI_Net=do.call(rbind.data.frame,RES)

Final_GGI_Net=GGI_Net[-which(is.na(GGI_Net$TGs_All)),]
Final_AllGI=unique(Final_GGI_Net$GI_All)
Final_AllTG=unique(Final_GGI_Net$TGs_All)

save(Final_GGI_Net,Final_AllGI,Final_AllTG,file="./Results/GGI_Net_Info.RData") ####

source(Find_AllConnect_Gen.R)
Find_AllConnect_Gen(Final_GGI_Net,Final_AllGI,Final_AllTG)