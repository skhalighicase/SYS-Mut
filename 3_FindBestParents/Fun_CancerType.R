## Compare and convert to the approved Gene

Fun_CancerType<-function(){
  split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x))) #Function to return the Tissue name via the forder name
  Folder_names=list.dirs(path = CurrentDir, full.names = TRUE, recursive = TRUE)
  Splitted_Path=split_path(Folder_names[1])
  Cancer_Name=Splitted_Path[3]
  CANCER_TYPE=strsplit(Cancer_Name,"_")[[1]][2]
  
  return(CANCER_TYPE)
}