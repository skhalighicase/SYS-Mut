## Last Modification 2/5/2020
## SYS_Mut Find Best Parents

graphics.off() 
rm(list=ls())
start_time <- Sys.time()

################## Setup parallel backend to use many processors
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
#################
library(stringi)
library(stringr)
library(RTCGAToolbox)
library(lattice)
library(gridExtra)
library(gapminder)
library(grid)
library(fpc)
library(DescTools)

################## Update and Setpath
CurrentDir=getwd()
setwd(CurrentDir)

###############
start_time <- Sys.time()

source("../1_PreProcess/Init_Parameters.R")
InitParams<-Init_Parameters()

#--------------- 
CANCER_TYPE<<-InitParams$CANCER

source("Parameters.R")
Params<-Parameters()

#--------------- Setup parallel back-end to use many processors
cl <- makeCluster(Params$cores) #not to overload your computer
registerDoParallel(cl)

################## Load files 
## Load mutation File for the current cancer type
load(Params$MutFile_Name)

################### Load Combined Data
load(Params$RNA_File_Name)
colnames(Nor_Comm_Data1$Nor_Comm_RNA1)=substr(colnames(Nor_Comm_Data1$Nor_Comm_RNA1),1,16)
RNA<<-Nor_Comm_Data1$Nor_Comm_RNA1
colnames(Nor_Comm_Data1$Nor_Comm_RNA1_GI)=substr(colnames(Nor_Comm_Data1$Nor_Comm_RNA1_GI),1,16)
RNA_Data_GI=Nor_Comm_Data1$Nor_Comm_RNA1_GI

################## Load CNV data ## Load Methylation data
load(Params$Trans_File_Name)
load(Params$Sig_File_Name)

CNV_Methyl=list()
CNV_Methyl$CNV=Nor_Comm_Trans_Data1$Nor_ThreshCNV1
colnames(Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1)=substr(colnames(Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1),1,16)
CNV_Methyl$Meth_Mean=Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1
CNV_Methyl<<-CNV_Methyl
rm(Nor_Comm_Trans_Data1,Nor_Comm_Data1,Nor_Comm_Sig_Data1)

################## Load Influence Network (New Influence Network)
load(Params$SYS_Net)
Influ_net<<-Influ_net

################# All TGs and their PGs and the Selected GIs connections
load(Params$GI_TG_PG)
Final_GI_TG_PG<<-Final_GI_TG_PG

## loading the final GI that we found in the connection network, Note this file returns all the connections between each two gene
load(Params$GGI_Net_Info) 

################# Filtering some parents Genes (Just as the parents of the target genes).
load(Params$WO_Zscore_File_Name) ############
AllGene_Info<-as.data.frame(AllGene_Inf_CancerSamp) ## If class(AllGene_Inf_CancerSamp) is not dataframe we should change it to
rm(AllGene_Inf_CancerSamp)
Pot_Idx=(which(!((AllGene_Info$Mean_Genes<0.5)&(AllGene_Info$Var_Genes<1)))) #Parents with mean expression value<0.5 and with Var<1 should remove from the list.
AllGene_Info<<-AllGene_Info[Pot_Idx,]

################# Selected GIs and all TGs
All_GIs<-as.character(unique(Final_GI_TG_PG$GI))
Current_Mut_Data<<-Current_Mut_Data
Exist_TG<<-as.character(unique(Final_GI_TG_PG$TG))

################
SamplesAll<<-data.frame("Samples_Name"=colnames(RNA_Data_GI))

################# Group gene analysis
MaxleGenInnterest<-length(All_GIs)

################# Connection Types
## For all Genes that have Expression data We will do the analysis :: All_GIs
Trans_Types<<-c("-t","-t|","-t>")
Activa_Type<<-c("-a|","-a>")

################ Runing the model for all the GIs that they have TG based on the influence netwrok

for(GIIdx in 1:length(All_GIs)){
  RES=array(list(),dim=1)
  Res_IDX<<-1
  Result<<-array(list())
  GI<-All_GIs[GIIdx]
  Flags_Seen<<-GI
  Connection_Flag<<-NULL

  ############## Selected TGs
  Selec_GI_TGs<<-as.character(Final_GI_TG_PG$TG[which(!is.na(match(Final_GI_TG_PG$GI,GI)))])

  ############## Mutation Info ## A function to read and Lable the NonSyn-Mutations
  source("Fun_Mutation_Labels.R")
  Mutation_Info<<-Fun_Mutation_Labels(GI)  ##
  ############## Find Target Genes of GI
  source("TGs_Find.R")
  TGs_GI<-TGs_Find(GI)

  ############## If GI does not have any target gene goes to else
  if(length(TGs_GI)){
    for(TG_Ele_Indx in 1:length(TGs_GI)){
      ######################
      switch(TG_Ele_Indx,
             "1"={  ## Running the model for all Normal Genes
               if (!is.null(TGs_GI[[TG_Ele_Indx]])){
                 Hyper_GI<-NA
                 TGs=as.vector(TGs_GI[[TG_Ele_Indx]]$Gene_Elements)
                 ## Expression GI or Hyper_GI
                 RNA_GI<-RNA_Data_GI[GI,]
                 source("Normal_Ele.R")
                 Result[[Res_IDX]]<-Normal_Ele(TGs,RNA_GI,GI,Hyper_GI)
                 Res_IDX<-Res_IDX+1
               }
             },
             "2"={  ## Running the model in all Complex cases
               if (!is.null(TGs_GI[[TG_Ele_Indx]])){
                 Hyper_GI<-GI
                 RNA_GI<-RNA_Data_GI[Hyper_GI,]
                 Flag_Complex_Seen<<-as.vector(TGs_GI[[TG_Ele_Indx]]$Complex_Elements)
                 source("Complex_Caclutation.R")
                 Complex_Caclutation(TGs_GI[[TG_Ele_Indx]],RNA_GI,Hyper_GI)
               }
             }
      )
    }
  }
  RES[[GIIdx]]=Result
  Result=NULL
  Flag_Complex_Seen=NULL
  Selec_GI_TGs=NULL
  Saved_File_name=paste0("./Results/Res_BestPar_",CANCER_TYPE,"_",GIIdx,".RData")
  save(RES,file=Saved_File_name)
}

source("Best_Parents.R")
Best_Parents()

stopCluster(cl)
closeAllConnections()
end_time <- Sys.time()
print(end_time - start_time)