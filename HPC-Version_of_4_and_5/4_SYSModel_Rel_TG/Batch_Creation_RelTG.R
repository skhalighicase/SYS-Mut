########## Generate some folder with a copy of the Template folder
Origindir <- c("/Template")

load("./Data/GGI_Net_Info.RData")
Maxfolders=ceiling(length(Final_AllGI)/2)

Currentpath=getwd()
Path_Orig_Folder=paste0(Currentpath,Origindir)
filestocopy <- list.files(path=Path_Orig_Folder,pattern = "\\.*$")

for(i in 1:Maxfolders){
  targetdir <- paste0(Currentpath,"/",i)
  if(dir.exists(paths = targetdir) == FALSE){
    dir.create(targetdir)
  }
  
  lapply(filestocopy, function(x) file.copy(paste (Path_Orig_Folder, x , sep = "/"),paste (targetdir,x, sep = "/"), recursive = FALSE,  copy.mode = TRUE))
  
  loadSlurm=paste0(targetdir,"/RUNMain.slurm")
  slurm_File_Read <- readLines(loadSlurm)
  slurm_File_Read[2]=paste0(slurm_File_Read[2],i)
  Splited_path=strsplit(loadSlurm,"/")
  Subgroup=tail(unlist(Splited_path), n=3)[-3]
  Cur_Folder=paste0(Currentpath,"/",i,"/Main_SYS_Model.R")
  Cur_Folder_CD=paste0(Currentpath,"/",i,"/")
  slurm_File_Read[18]=paste(slurm_File_Read[18],Cur_Folder_CD)
  slurm_File_Read[19]=paste(slurm_File_Read[19],Cur_Folder)
  writeLines(slurm_File_Read, loadSlurm)
}

#### For Best_Parent_Analysis
Res_Dir <- paste0(Currentpath,"/Results")

lapply("RUNMain.slurm" , function(x) file.copy(paste (Path_Orig_Folder, x , sep = "/"),paste (Res_Dir,x, sep = "/"), recursive = FALSE,  copy.mode = TRUE))

loadSlurm=paste0(Res_Dir,"/RUNMain.slurm") 
slurm_File_Read <- readLines(loadSlurm)
slurm_File_Read[2]=paste0(slurm_File_Read[2],"Analysis")
Splited_path=strsplit(loadSlurm,"/")
Subgroup=tail(unlist(Splited_path), n=3)[-3]
Cur_Folder=paste0(Currentpath,"/Results/Calc_Results_Measures.R")
Cur_Folder_CD=paste0(Currentpath,"/Results/")
slurm_File_Read[18]=paste(slurm_File_Read[18],Cur_Folder_CD)
slurm_File_Read[19]=paste(slurm_File_Read[19],Cur_Folder)
writeLines(slurm_File_Read, loadSlurm)