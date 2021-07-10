########## Generate some folder with a copy of the Template folder
Gen_File_Name="./Data/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt"
Mutated_Genes=read.table(Gen_File_Name, header = TRUE, sep = "\t") 
All_Mutated_GIs <-length(Mutated_Genes$GI_ALL[Mutated_Genes$Perce_Impact_Final>0])
#Maxfolders=All_Mutated_GIs 
Maxfolders=ceiling(All_Mutated_GIs/10)


CurrentPath=getwd()
Folder_names=list.dirs(path = CurrentPath, full.names = TRUE, recursive = TRUE)
Path_Orig_Folder=paste0(Folder_names,"/RUNMain.slurm")

for(i in 2:(Maxfolders+1)){
    Command=paste("sbatch", Path_Orig_Folder[i])
    system(Command)
  }
