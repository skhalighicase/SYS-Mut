## Last Modification 6/6/2021
## SYS-Mut Real-TG Running the model

graphics.off() 
rm(list=ls())
start_time <- Sys.time()

#--------------- Setup parallel back end to use many processors
library(iterators)
library(parallel)
library(foreach)
library(doParallel)

#----------------
library(stringi)
library(stringr)
library(RTCGAToolbox)
library(lattice)
library(gridExtra)
library(gapminder)
library(grid)
library(fpc)
library(DescTools)
library(data.table)
#---------------- Set-path
CurrentDir=getwd()
setwd(CurrentDir)

#--------------- Parameters
source("Parameters.R")
Params<-Parameters()

#--------------- 
CANCER_TYPE<<-Params$CANCER

#--------------- Setup parallel back-end to use many processors
cl <- makeCluster(Params$cores) #not to overload your computer
registerDoParallel(cl)

################## Load files ###
## Load mutation File for the current cancer type
load (Params$MutFile_Name)
Current_Mut_Data<<-Current_Mut_Data

################### Load Combined Data
load(Params$RNA_File_Name)
colnames(Nor_Comm_Data1$Nor_Comm_RNA1)=substr(colnames(Nor_Comm_Data1$Nor_Comm_RNA1),1,16)
RNA<<-Nor_Comm_Data1$Nor_Comm_RNA1
colnames(Nor_Comm_Data1$Nor_Comm_RNA1_GI)=substr(colnames(Nor_Comm_Data1$Nor_Comm_RNA1_GI),1,16)
RNA_Data_GI=Nor_Comm_Data1$Nor_Comm_RNA1_GI

################### load CNV data ## Load Methylation data
load(Params$Trans_File_Name)
load(Params$Sig_File_Name)

CNV_Methyl=list()
CNV_Methyl$CNV=Nor_Comm_Trans_Data1$Nor_ThreshCNV1
colnames(Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1)=substr(colnames(Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1),1,16)
CNV_Methyl$Meth_Mean=Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1
CNV_Methyl<<-CNV_Methyl
rm(Nor_Comm_Trans_Data1,Nor_Comm_Data1,Nor_Comm_Sig_Data1)

################# Load Influence Network (New Influence Network)
load(Params$SYSMut_InfNet)
Influ_net<<-Influ_net

############## All Connections between each TG and It's Parent Genes
load(Params$TG_PG_Connections)  ### All Parents without Filter TG_Parent_Connection
TG_Parents_Connections<<-TG_Parents_Connections

############## Loading the final GI that we found in the connection network
load(Params$GGI_Net_Info) ## All GI,TG and All Connections #Final_GG_Net
Final_GGI_Net<<-Final_GGI_Net
Exist_TG<<-Final_AllTG
rm(Final_AllGI,Final_AllTG)

############## Load Parent Information (Best Parents)
load(Params$Best_Parents_File)
All_Info_TG_Par<<-All_Info_TG_Par #(Best Parent)

############## All GIs (universe) for the Random TGs. We got this list of GIs from the Run of Filtered real TG
Mutated_Genes=read.table(Params$Rel_GI_Summry, header = TRUE, sep = "\t")

All_GIs <-as.character(Mutated_Genes$GI_ALL[Mutated_Genes$Perce_Impact_Final>0])

MaxleGenInnterest<-length(All_GIs)

##----------------- Connection Types
## For all Genes that have Expression data We will do the analysis :: All_GIs
Trans_Types<<-c("-t","-t|","-t>")
Activa_Type<<-c("-a|","-a>")
if(length(All_GIs)>0){
  #############################
  for(GIIdx in 1:length(All_GIs)){ ####
    RES=array(list(),dim=1)
    Res_IDX<<-1
    Result<<-array(list())
    GI<-All_GIs[GIIdx]
    Flags_Seen<<-GI
    Flags_Seen_OldTG<<-GI
    Connection_Flag<<-NULL
    
    ############## Find Target Genes of GI
    source("TGs_Find.R")
    TGs_GI<-TGs_Find(GI)
    
    ### If GI does not have any target gene goes to else
    if(length(TGs_GI)){
      for(TG_Ele_Indx in 1:length(TGs_GI)){
        ######################
        switch(TG_Ele_Indx,
               "1"={  ## Running the model for all Normal Genes
                 if (!is.null(TGs_GI[[TG_Ele_Indx]])){
                   Hyper_GI<-NA
                   TGs=as.character(TGs_GI[[TG_Ele_Indx]]$Gene_Elements)
                   ## Expression GI or Hyper_GI
                   RNA_GI<-RNA_Data_GI[GI,]
                   source("Normal_Ele_Rnd.R")
                   Result[[Res_IDX]]<-Normal_Ele_Rnd(TGs,RNA_GI,GI,Hyper_GI)
                   Res_IDX<-Res_IDX+1
                 }
               },
               "2"={  ## Running the model in all Complex cases
                 if (!is.null(TGs_GI[[TG_Ele_Indx]])){
                   Hyper_GI<-GI
                   RNA_GI<-RNA_Data_GI[Hyper_GI,]
                   Flag_Complex_Seen<<-as.vector(TGs_GI[[TG_Ele_Indx]]$Complex_Elements)
                   source("Complex_Calc_Rnd.R")
                   Complex_Calc_Rnd(TGs_GI[[TG_Ele_Indx]],RNA_GI,Hyper_GI)
                 }
               }
        )
      }
    }
    RES[[GIIdx]]=Result
    Result=NULL
    Flag_Complex_Seen=NULL
    Flags_Seen=NULL
    Flags_Seen_OldTG=NULL
    Saved_File_name=paste0("../Results/Res_NonSyn_Wild_",GIIdx,"_",GI,".RData")
    save(RES,file=Saved_File_name)
  }
  
  #####################
  source("Calc_Results_Measure.R")
  Calc_Results_Measure()
  
  setpath="../Template/"
  setwd(setpath)
  
  source("Analysis_Impacts_New.R")
  Analysis_Impacts_New()
  
}else{
  print("None of GIs have Impact")
}

stopCluster(cl)
closeAllConnections()
end_time <- Sys.time()
print(end_time - start_time)