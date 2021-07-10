Charts_Rnd_Real_WI_NewFilter<-function(){
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(data.table)
  ################## Load Mutation File
  Mutaiton_file=paste0("../1_PreProcess/InputData/Mut_WO_Hyper_HGNC_",CANCER_TYPE,".RData") 
  load(Mutaiton_file)
  Mutation_Rates=Mutt_Current_Tissue
  rm(Mutt_Current_Tissue)
  
  ################## Sort Filtered summary by Gene Names
  Summryfiles_Rl <-paste0("../4_SYSModel_Rel_TG/Analysis/Res_NonSyn_Wild_Summry_GI_Filtered_RelTG_NewFilter.txt")
  Filter_Summ_Rl=read.delim(Summryfiles_Rl,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Filter_Summ_Rl <- Filter_Summ_Rl[order(Filter_Summ_Rl$GI_ALL),] 
  Filtered_GenLst=Filter_Summ_Rl$GI_ALL
  write.table(Filter_Summ_Rl, "./Results/Sort_Filt_Summ_RL_NewFilter.txt", sep="\t",row.names = FALSE)
  
  #################  Sort Filtered summary by Gene Names
  Summryfiles_Rn <-paste0("./Results/Avg_Perc_Impacted_RndTG_",CANCER_TYPE,".txt")
  Filter_Summ_Rn=read.delim(Summryfiles_Rn,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Filter_Summ_Rn=Filter_Summ_Rn[which(!is.na(match(Filter_Summ_Rn$GI_ALL,Filtered_GenLst))),]

  Filter_Summ_Rn <- Filter_Summ_Rn[order(Filter_Summ_Rn$GI_ALL),]
  write.table(Filter_Summ_Rn, "./Results/Sort_Filt_Summ_RN.txt", sep="\t",row.names = FALSE)
  
  ########################################################################
  ################# RANDOM TARGET GENES ### Plot All Random_TG for each GI, 
  ## ALL Mutation Rates##
  Filter_Summ_All=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers,Perce_Impact_Final,GI_TG_total)
  
  Filter_Summ_All$TG_All_Categ=NA
  Filter_Summ_All$TG_All_Categ[which(Filter_Summ_All$GI_TG_total<10)]=10
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=10)&(Filter_Summ_All$GI_TG_total<25))]=25
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=25)&(Filter_Summ_All$GI_TG_total<50))]=50
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=50)&(Filter_Summ_All$GI_TG_total<100))]=100
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=100))]=200
  
  Filter_Summ_All$TG_All_Categ=as.factor(Filter_Summ_All$TG_All_Categ)
  
  nbaplot_all=ggplot(Filter_Summ_All, aes(x=MutRate_Univers, y=Perce_Impact_Final,label=GI_ALL))+geom_point(aes(size=TG_All_Categ,color=TG_All_Categ))+
    coord_cartesian(ylim = c(0, 100))+ xlab("Mutation Rate of All Mutated Genes") + ylab("Percentage of Impacted Random_TGs")+
    labs(size="Number of Downstream Target Genes",colour="Number of Downstream Target Genes")+
    theme(legend.position="top")
  
  file_Plot=paste0("./Results/Percent_Rnd_TG_All.png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot_all)
  
  ##################################################################################################################
  ############### Plot All Random_TG for each GI, Less than 5% 
  Filter_Summ_All=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers,Perce_Impact_Final,GI_TG_total)
  Filter_Summ_All=Filter_Summ_All[-which(Filter_Summ_All$MutRate_Univers>5),]
  # Filter_Summ_All=Filter_Summ_All[-which(Filter_Summ_All$Perce_Impact_Final<25),]
  
  Filter_Summ_All$TG_All_Categ=NA
  Filter_Summ_All$TG_All_Categ[which(Filter_Summ_All$GI_TG_total<10)]=10
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=10)&(Filter_Summ_All$GI_TG_total<25))]=25
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=25)&(Filter_Summ_All$GI_TG_total<50))]=50
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=50)&(Filter_Summ_All$GI_TG_total<100))]=100
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=100))]=200
  #### 
  Filter_Summ_All$TG_All_Categ=as.factor(Filter_Summ_All$TG_All_Categ)
  
  #### Change the point size
  nbaplot5=ggplot(Filter_Summ_All, aes(x=MutRate_Univers, y=Perce_Impact_Final,label=GI_ALL))+geom_point(aes(size=TG_All_Categ,color=TG_All_Categ))+coord_cartesian(ylim = c(0, 100))+ 
    xlab("Mutation Rate of All Mutated Genes, Less 5%") + ylab("Percentage of Impacted Random_TGs")
  
  file_Plot=paste0("./Results/Percent_Rnd_TG_LessThan_5",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot5)
  
  ##################################################################################################################
  #### Plot All Direct Random_TG of each GI,  
  Filter_Summ_Dir=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers,Perce_Impact_Dir,TG_Dir_WI,TG_Dir_WO)
  Filter_Summ_Dir$TG_Dir_All= Filter_Summ_Dir$TG_Dir_WI+Filter_Summ_Dir$TG_Dir_WO
  # Filter_Summ_Dir=Filter_Summ_Dir[-which(Filter_Summ_Dir$Perce_Impact_Dir<25),]
  
  Filter_Summ_Dir$TG_Dir_Categ=NA
  Filter_Summ_Dir$TG_Dir_Categ[which(Filter_Summ_Dir$TG_Dir_All<10)]=10
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=10)&(Filter_Summ_Dir$TG_Dir_All<25))]=25
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=25)&(Filter_Summ_Dir$TG_Dir_All<50))]=50
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=50)&(Filter_Summ_Dir$TG_Dir_All<100))]=100
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=100))]=200
  Filter_Summ_Dir$TG_Dir_Categ=as.factor(Filter_Summ_Dir$TG_Dir_Categ)
  
  nbaplot_Dir=ggplot(Filter_Summ_Dir, aes(x=MutRate_Univers, y=Perce_Impact_Dir,label=GI_ALL))+geom_point(aes(size=TG_Dir_Categ,color=TG_Dir_Categ))+coord_cartesian(ylim = c(0, 100))+
    xlab("Mutation Rate of All Mutated Genes") +  ylab("Percentage of Impacted Direct TG")
  
  file_Plot=paste0("./Results/Percent_Rnd_TG_All_Direct",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot_Dir)
  
  #########################################################################################
  #### Plot All Direct Random_TG of each GI less 5%,  
  
  Filter_Summ_Dir=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers,Perce_Impact_Dir,TG_Dir_WI,TG_Dir_WO)
  Filter_Summ_Dir$TG_Dir_All= Filter_Summ_Dir$TG_Dir_WI+Filter_Summ_Dir$TG_Dir_WO
  Filter_Summ_Dir=Filter_Summ_Dir[-which(Filter_Summ_Dir$MutRate_Univers>5),]
  # Filter_Summ_Dir=Filter_Summ_Dir[-which(Filter_Summ_Dir$Perce_Impact_Dir<25),]
  
  Filter_Summ_Dir$TG_Dir_Categ=NA
  Filter_Summ_Dir$TG_Dir_Categ[which(Filter_Summ_Dir$TG_Dir_All<10)]=10
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=10)&(Filter_Summ_Dir$TG_Dir_All<25))]=25
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=25)&(Filter_Summ_Dir$TG_Dir_All<50))]=50
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=50)&(Filter_Summ_Dir$TG_Dir_All<100))]=100
  Filter_Summ_Dir$TG_Dir_Categ[which((Filter_Summ_Dir$TG_Dir_All>=100))]=200
  Filter_Summ_Dir$TG_Dir_Categ=as.factor(Filter_Summ_Dir$TG_Dir_Categ)
  
  nbaplot_Dir5=ggplot(Filter_Summ_Dir, aes(x=MutRate_Univers, y=Perce_Impact_Dir,label=GI_ALL))+geom_point(aes(size=TG_Dir_Categ,color=TG_Dir_Categ))+coord_cartesian(ylim = c(0, 100))+
    xlab("Mutation Rate of All Mutated Genes, Less Than 5%") +  ylab("Percentage of Impacted Direct TG")
  
  ############################ 
  file_Plot=paste0("./Results/Percent_Rnd_TG_All_LessThann_5_Direct",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot_Dir5)
  
  ######################################################################################
  ## Plot All InDirect Random_TG of each GI,  
  
  Filter_Summ_InInDir=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers, Perce_Impact_InDir,TG_InDir_WI,TG_InDir_WO)
  Filter_Summ_InInDir$TG_InDir_All= Filter_Summ_InInDir$TG_InDir_WI+Filter_Summ_InInDir$TG_InDir_WO
  
  Filter_Summ_InInDir$TG_Dir_Categ=NA
  Filter_Summ_InInDir$TG_Dir_Categ[which(Filter_Summ_InInDir$TG_InDir_All<10)]=10
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=10)&(Filter_Summ_InInDir$TG_InDir_All<25))]=25
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=25)&(Filter_Summ_InInDir$TG_InDir_All<50))]=50
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=50)&(Filter_Summ_InInDir$TG_InDir_All<100))]=100
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=100))]=200
  Filter_Summ_InInDir$TG_Dir_Categ=as.factor(Filter_Summ_InInDir$TG_Dir_Categ)
  
  # Change the point size
  nbaplot_Indir=ggplot(Filter_Summ_InInDir, aes(x=MutRate_Univers, y=Perce_Impact_InDir,label=GI_ALL))+geom_point(aes(size=TG_Dir_Categ,color=TG_Dir_Categ))+coord_cartesian(ylim = c(0, 100))+
    xlab("Mutation Rate of All Mutated Genes") +  ylab("Percentage of Impacted InDirect Random TG, Level 2 ")
  
  file_Plot=paste0("./Results/Percent_Rnd_TG_All_Indirect",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot_Indir)
  
  #################################
  #### Plot InDirect Random_TG of each GI, Just rare mutation<5%
  
  Filter_Summ_InInDir=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers, Perce_Impact_InDir, TG_InDir_WI,TG_InDir_WO)
  Filter_Summ_InInDir$TG_InDir_All= Filter_Summ_InInDir$TG_InDir_WI+Filter_Summ_InInDir$TG_InDir_WO
  Filter_Summ_InInDir=Filter_Summ_InInDir[-which(Filter_Summ_InInDir$MutRate_Univers>5),]
  # Filter_Summ_InInDir=Filter_Summ_InInDir[-which(Filter_Summ_InInDir$Perce_Impact_InDir<25),]
  
  Filter_Summ_InInDir$TG_Dir_Categ=NA
  Filter_Summ_InInDir$TG_Dir_Categ[which(Filter_Summ_InInDir$TG_InDir_All<10)]=10
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=10)&(Filter_Summ_InInDir$TG_InDir_All<25))]=25
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=25)&(Filter_Summ_InInDir$TG_InDir_All<50))]=50
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=50)&(Filter_Summ_InInDir$TG_InDir_All<100))]=100
  Filter_Summ_InInDir$TG_Dir_Categ[which((Filter_Summ_InInDir$TG_InDir_All>=100))]=200
  Filter_Summ_InInDir$TG_Dir_Categ=as.factor(Filter_Summ_InInDir$TG_Dir_Categ)
  
  # Change the point size
  nbaplot_Indir5=ggplot(Filter_Summ_InInDir, aes(x=MutRate_Univers, y=Perce_Impact_InDir,label=GI_ALL))+geom_point(aes(size=TG_Dir_Categ,color=TG_Dir_Categ))+coord_cartesian(ylim = c(0, 100))+
    xlab("Mutation Rate of All Mutated Genes, Less 5%") +  ylab("Percentage of Impacted InDirect Random TG, Level 2 ")
  
  file_Plot=paste0("./Results/Percent_Rnd_TG_All_LessThann_5_Indirect",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot_Indir5)
  
  ########################################################
  ## Plot Highimpact All-TG
  
  Filter_Summ_All=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers, Percent_All_High, GI_TG_total)
  
  Filter_Summ_All$TG_All_Categ=NA
  Filter_Summ_All$TG_All_Categ[which(Filter_Summ_All$GI_TG_total<10)]=10
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=10)&(Filter_Summ_All$GI_TG_total<25))]=25
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=25)&(Filter_Summ_All$GI_TG_total<50))]=50
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=50)&(Filter_Summ_All$GI_TG_total<100))]=100
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=100))]=200
  
  Filter_Summ_All$TG_All_Categ=as.factor(Filter_Summ_All$TG_All_Categ)
  
  # Change the point size
  nbaplot3=ggplot(Filter_Summ_All, aes(x=MutRate_Univers, y=Percent_All_High,label=GI_ALL))+geom_point(aes(size=TG_All_Categ,color=TG_All_Categ))+coord_cartesian(ylim = c(0, 100))+
    xlab("Mutation Rate of All Mutated Genes") +  ylab("Percentage of High Impacte Direct and Indirect TG")
  
  file_Plot=paste0("./Results/Percent_HighImpact_Rnd_All",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot3)
  
  ##########################################################
  #### Plot Highimpact All-TG, Less 5%
  
  Filter_Summ_All=dplyr::select(Filter_Summ_Rn,GI_ALL,MutRate_Univers, Percent_All_High, GI_TG_total)
  Filter_Summ_All=Filter_Summ_All[-which(Filter_Summ_All$MutRate_Univers>5),]
  
  Filter_Summ_All$TG_All_Categ=NA
  Filter_Summ_All$TG_All_Categ[which(Filter_Summ_All$GI_TG_total<10)]=10
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=10)&(Filter_Summ_All$GI_TG_total<25))]=25
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=25)&(Filter_Summ_All$GI_TG_total<50))]=50
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=50)&(Filter_Summ_All$GI_TG_total<100))]=100
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=100))]=200
  
  Filter_Summ_All$TG_All_Categ=as.factor(Filter_Summ_All$TG_All_Categ)
  
  # Change the point size
  nbaplot3_5=ggplot(Filter_Summ_All, aes(x=MutRate_Univers, y=Percent_All_High,label=GI_ALL))+geom_point(aes(size=TG_All_Categ,color=TG_All_Categ))+coord_cartesian(ylim = c(0, 100))+
    xlab("Mutation Rate of All Mutated Genes, Less 5%") +  ylab("Percentage of High Impacte Direct and Indirect TG")
  
  file_Plot=paste0("./Results/HighImpact_Rnd_LessThan_5",CANCER_TYPE,".png")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 4000 / dpi
  height.calc <- 2000 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot3_5)
}