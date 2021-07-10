########## Generate some folder with a copy of the Template folder
load("./Data/GGI_Net_Info.RData")
Maxfolders=ceiling(length(Final_AllGI)/2)

CurrentPath=getwd()
Folder_names=list.dirs(path = CurrentPath, full.names = TRUE, recursive = TRUE)
Path_Orig_Folder=paste0(Folder_names,"/RUNMain.slurm")

for(i in 2:(Maxfolders+1)){
    Command=paste("sbatch", Path_Orig_Folder[i])
    system(Command)
  }
