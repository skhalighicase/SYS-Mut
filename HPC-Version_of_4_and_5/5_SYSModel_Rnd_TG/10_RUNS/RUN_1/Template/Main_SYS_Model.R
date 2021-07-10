## Last Modification 5/31/2021
## SYS_Mut Random-TG
## to do: 
## to do: 
## to do: 
graphics.off() 
rm(list=ls())
start_time <- Sys.time()

################## Setup parallel backend to use many processors
library(iterators)
library(parallel)
library(foreach)
library(doParallel)
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

############### Cancer_Name From Folder
split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x))) #Function to return the Tissue name via the forder name
Folder_names=list.dirs(path = CurrentDir, full.names = TRUE, recursive = TRUE)
Splitted_Path=split_path(Folder_names[1])
Cancer_Name=Splitted_Path[5]
CANCER_TYPE=strsplit(Cancer_Name,"_")[[1]][2]

################## Parameters
source("Parameters.R")
Params<-Parameters()

################## Setup parallel backend to use many processors
cl <- makeCluster(Params$cores) #not to overload your computer
registerDoParallel(cl)

################## Load files ###
## Load mutation File for the current cancer type
MutFile_Name=paste0("../Data/Updated_Mut_",CANCER_TYPE,".RData")  ## This file is created based on list of the RNAsaq samples.
load (MutFile_Name)
Current_Mut_Data<<-Current_Mut_Data

################### Load Combined Data
RNA_File_Name=paste0("../Data/Combin_Gen_Samp_RNA_",CANCER_TYPE,".RData")
load(RNA_File_Name)
colnames(Nor_Comm_Data1$Nor_Comm_RNA1)=substr(colnames(Nor_Comm_Data1$Nor_Comm_RNA1),1,16)
RNA<<-Nor_Comm_Data1$Nor_Comm_RNA1
colnames(Nor_Comm_Data1$Nor_Comm_RNA1_GI)=substr(colnames(Nor_Comm_Data1$Nor_Comm_RNA1_GI),1,16)
RNA_Data_GI=Nor_Comm_Data1$Nor_Comm_RNA1_GI

################### load CNV data ## Load Methylation data
Trans_File_Name=paste0("../Data/Combin_Gen_Samp_Trans_",CANCER_TYPE,".RData")
load(Trans_File_Name)
Sig_File_Name=paste0("../Data/Combin_Gen_Samp_Sig_",CANCER_TYPE,".RData")
load(Sig_File_Name)

CNV_Methyl=list()
CNV_Methyl$CNV=Nor_Comm_Trans_Data1$Nor_ThreshCNV1
colnames(Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1)=substr(colnames(Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1),1,16)
CNV_Methyl$Meth_Mean=Nor_Comm_Sig_Data1$Nor_Sig_Meth_Mean1
CNV_Methyl<<-CNV_Methyl
rm(Nor_Comm_Trans_Data1,Nor_Comm_Data1,Nor_Comm_Sig_Data1)

################# Load Influence Network (New Influence Network)
load("../Data/SYS_Net.RData")
Influ_net<<-Influ_net

############## All Connections between each TG and It's Parent Genes 
load("../Data/TG_PGs_Connec.RData")  ### All Parents without Filter TG_Parent_Connection  
TG_Parents_Connections<<-TG_Parents_Connections

############## Load Parent Information (Best Parents)
Best_Parents_File_Name=paste0("../Data/Results_TG_Parents_",CANCER_TYPE,".RData")
load(Best_Parents_File_Name)
All_Info_TG_Par<<-All_Info_TG_Par # (Best Parent)

############## Loading the final GI that we found in the connection network
load("../Data/GGI_Net_Info.RData") ## All GI,TG and All Connections #Final_GG_Net
Final_GGI_Net<<-Final_GGI_Net
Exist_TG<<-Final_AllTG
rm(Final_AllGI,Final_AllTG)

############## All GIs (universe) for the Random TGs. We got this list of GIs from the Run of Filtered real TG
Gen_File_Name=paste0("../Data/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt")
Mutated_Genes=read.table(Gen_File_Name, header = TRUE, sep = "\t") 

All_GIs <-as.character(Mutated_Genes$GI_ALL[Mutated_Genes$Perce_Impact_Final>0])
###----------------- Group gene analysis (Devide into subgroups of GIs)
GenInterest=(as.numeric(Params$Subgroup)-1)*Params$Subgroup_Leng+1  
MaxleGenInnterest=as.numeric(Params$Subgroup)*Params$Subgroup_Leng

if(MaxleGenInnterest>length(All_GIs)){
  MaxleGenInnterest<-length(All_GIs)
}

###----------------- Connection Types
## For all Genes that have Expression data We will do the analysis :: All_GIs
Trans_Types<<-c("-t","-t|","-t>")
Activa_Type<<-c("-a|","-a>")

#############################
RES=array(list(),dim=MaxleGenInnterest-GenInterest+1)
for(GIIdx in GenInterest:MaxleGenInnterest){ #### length(All_GIs)  
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
}
Saved_File_name=paste0("../Results/Res_NonSyn_Wild_",GenInterest,"_",MaxleGenInnterest,".RData")
save(RES,file=Saved_File_name)
stopCluster(cl)
closeAllConnections()
end_time <- Sys.time()
print(end_time - start_time)
