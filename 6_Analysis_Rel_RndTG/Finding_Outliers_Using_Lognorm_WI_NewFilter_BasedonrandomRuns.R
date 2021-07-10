## Lat Update 6/9/2020
Finding_Outliers_Using_Lognorm_WI_NewFilter_BasedonrandomRuns<-function(){
  browser()
  ################# Average Results of Random Run to compare and to estimate the threshold line
  Summryfiles_Rn <-paste0("./Results/Avg_Perc_Impacted_RndTG_",CANCER_TYPE,".txt")
  Res_Rnd_TG=read.delim(Summryfiles_Rn,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  ################## Find the Threshold based on 95% of Lognormal null distribution
  source("Thresh_Outlier.R")
  Fited_and_Data_quant<-Thresh_Outlier(Res_Rnd_TG)
  Threshold=Fited_and_Data_quant$data_quants[which(!is.na(match(rownames(Fited_and_Data_quant),"95%")))]
  ## Manual Threshold for Ourlier Selection
  
  # Threshold=16.5
  
  #########################################################
  ##--------- Plot all the datapoints of Real TGs--------##   
  Summryfiles <-"./Results/Sort_Filt_Summ_RL_NewFilter.txt"
  Filter_Summ_Rl=read.delim(Summryfiles,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  Filter_Summ_All=dplyr::select(Filter_Summ_Rl,GI_ALL,Raw_MutRate,Perce_Impact_Final,GI_TG_total)
  
  ##-----------------------------------------------------##
  ## Group the points based on the number of TGs (Make a Factor type)
  Filter_Summ_All$TG_All_Categ=NA
  Filter_Summ_All$TG_All_Categ[which(Filter_Summ_All$GI_TG_total<10)]=10
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=10)&(Filter_Summ_All$GI_TG_total<25))]=25
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=25)&(Filter_Summ_All$GI_TG_total<50))]=50
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=50)&(Filter_Summ_All$GI_TG_total<100))]=100
  Filter_Summ_All$TG_All_Categ[which((Filter_Summ_All$GI_TG_total>=100))]=200
  Filter_Summ_All$TG_All_Categ=as.factor(Filter_Summ_All$TG_All_Categ)
  
  Maximum=90
  ##--------------------- Plot --------------------------##
  library(scales)
  squish_trans <- function(from, to, factor) {
    
    trans <- function(x) {
      
      # get indices for the relevant regions
      isq <- x > from & x < to
      ito <- x >= to
      
      # apply transformation
      x[isq] <- from + (x[isq] - from)/factor
      x[ito] <- from + (to - from)/factor + (x[ito] - to)
      
      return(x)
    }
    
    inv <- function(x) {
      
      # get indices for the relevant regions
      isq <- x > from & x < from + (to - from)/factor
      ito <- x >= from + (to - from)/factor
      
      # apply transformation
      x[isq] <- from + (x[isq] - from) * factor
      x[ito] <- to + (x[ito] - (from + (to - from)/factor))
      
      return(x)
    }
    
    # return the transformation
    return(trans_new("squished", trans, inv))
  }
  
  #################
  nbaplot_All <- ggplot(data=Filter_Summ_All,aes(x=Raw_MutRate, y=Perce_Impact_Final,label="Name"))+
    geom_point(data=Filter_Summ_All,aes(size=TG_All_Categ,color=TG_All_Categ))+coord_cartesian(ylim = c(0, 100),xlim = c(0,100))+
    xlab("Mutation Rate of GI") + ylab("Percentage of Impacted TG, (Outlier Genes)")+
    geom_hline(aes(yintercept = Threshold),linetype="dotted", color = "blue", size=0.75)+
    annotate("text",Maximum-7,Threshold+2 , vjust = -1,label = "Outlier Threshold")+theme_bw()+theme(panel.grid = element_blank())+
    labs(size="Number of Downstream Target Genes",colour="Number of Downstream Target Genes")+
    theme(legend.position="top")
  
  nbaplot_All=nbaplot_All + scale_x_continuous(trans = squish_trans(5, 100,10),
                                               breaks = seq(0, 100, by =5))
  print(nbaplot_All)
  ##----------------- Save as a Figure-------------------##
  file_Plot=paste0("./Results/All_Genes_With_Thresholds_NewFilter_",CANCER_TYPE,"_RndTG.pdf")
  text.factor<-3
  dpi <- text.factor * 100
  width.calc <- 2000 / dpi
  height.calc <- 1500 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc,  units = 'in', plot = nbaplot_All)
  
  ####################################################
  ##--------------     Find and list the Outlier Genes
  Filter_Summ_Rl$Outlier=NA
  ## Select Threshold for outlier selection
  Curre_Thresh=Threshold
  for(I in 1:length(Filter_Summ_Rl$Raw_MutRate)){
    Filter_Summ_Rl$Outlier[I]=ifelse(Filter_Summ_Rl$Perce_Impact_Final[I]>=Curre_Thresh, 1,0)
  }
  Filter_Summ_Rl=Filter_Summ_Rl[order(-Filter_Summ_Rl$Perce_Impact_Final),]
  
  write.table(Filter_Summ_Rl, "./Results/Filter_Summ_Rl_NewFilter_RndTG.txt", sep="\t",row.names = FALSE)
  
  ############ save the Candidate Genes in Excel sheet
  Summryfiles_Cand_exl_all <-paste0("./Results/Filter_Summ_Rl_NewFilter_RndTG.xlsx")
  write_xlsx(Filter_Summ_Rl,Summryfiles_Cand_exl_all)
  
  ############# Save the Candidate Genes #####################
  Candidate_GenesInfo=Filter_Summ_Rl[which(Filter_Summ_Rl$Outlier==1),]
  Candidate_GenesInfo=Candidate_GenesInfo[order(-Candidate_GenesInfo$Perce_Impact_Final),]
  
  Summryfiles_Cand <-paste0("./Results/Candidate_GenesInfo_NewFilter_",CANCER_TYPE,"_RndTG.txt")
  write.table(Candidate_GenesInfo, Summryfiles_Cand, sep="\t",row.names = FALSE)
  
  ############ save the Candidate Genes in Excel sheet
  Summryfiles_Cand_exl <-paste0("./Results/Candidate_GenesInfo_NewFilter_",CANCER_TYPE,"_RndTG.xlsx")
  write_xlsx(Candidate_GenesInfo,Summryfiles_Cand_exl)
  
  ############ Read
  ## Sort Filtered summary by Gene Names
  Summryfiles_Cand <-paste0("./Results/Candidate_GenesInfo_NewFilter_",CANCER_TYPE,"_RndTG.txt")
  Candid_Gen_Summ=read.delim(Summryfiles_Cand,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  Filter_Summ_Candid=dplyr::select(Candid_Gen_Summ,GI_ALL,Raw_MutRate, Perce_Impact_Final, GI_TG_total)
  
  Filter_Summ_Candid$TG_All_Categ=NA
  Filter_Summ_Candid$TG_All_Categ[which(Filter_Summ_Candid$GI_TG_total<10)]=10
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=10)&(Filter_Summ_Candid$GI_TG_total<25))]=25
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=25)&(Filter_Summ_Candid$GI_TG_total<50))]=50
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=50)&(Filter_Summ_Candid$GI_TG_total<100))]=100
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=100))]=200
  
  Filter_Summ_Candid$TG_All_Categ=as.factor(Filter_Summ_Candid$TG_All_Categ)
  
  # Change the point size
  nbaplot3_detailed=ggplot(Filter_Summ_Candid, aes(x=Raw_MutRate, y=Perce_Impact_Final,label="GI_ALL"))+
    geom_point(data=Filter_Summ_Candid,aes(size=TG_All_Categ,color=TG_All_Categ))+coord_cartesian(ylim = c(0, 100),xlim = c(0,100))+
    xlab("Mutation Rate of GI") + ylab("Percentage of Impacted TG, (Outlier Genes)")+theme_bw()+theme(panel.grid = element_blank())+
    theme(legend.position="top")
  
  nbaplot3_detailed=nbaplot3_detailed + scale_x_continuous(trans = squish_trans(5, 100,10),
                                                           breaks = seq(0, 100, by =5))
  
  Labled_Plot_highimp_Det=nbaplot3_detailed + geom_label_repel(aes(label = GI_ALL),
                                                               box.padding   = 0.1,
                                                               point.padding = 0.1,
                                                               segment.color = 'grey50')+ 
    labs(size="Number of Downstream Target Genes",colour="Number of Downstream Target Genes")
  
  
  file_Plot=paste0("./Results/Labeled_Candidate_Genes_AllMutRate_NewFilter_",CANCER_TYPE,"_RndTG.pdf")
  text.factor<-2
  dpi <- text.factor * 100
  width.calc <- 2300 / dpi
  height.calc <- 1500 / dpi
  ggsave(filename = file_Plot, dpi = dpi, width = width.calc, height = height.calc , units = 'in', plot = Labled_Plot_highimp_Det)
  
  #######################################################################
  ######-----------Labeled Candidate genes, less than 5% --------------##
  ## Sort Filtered summary by Gene Names
  Summryfiles_Cand <-paste0("./Results/Candidate_GenesInfo_NewFilter_",CANCER_TYPE,"_RndTG.txt")
  Candid_Gen_Summ=read.delim(Summryfiles_Cand,header = TRUE, sep="\t",stringsAsFactors = FALSE,check.names = FALSE)
  
  Filter_Summ_Candid=dplyr::select(Candid_Gen_Summ,GI_ALL,Raw_MutRate, Perce_Impact_Final, GI_TG_total)
  Filter_Summ_Candid=Filter_Summ_Candid[-which(Filter_Summ_Candid$Raw_MutRate>5),]
  
  Filter_Summ_Candid$TG_All_Categ=NA
  Filter_Summ_Candid$TG_All_Categ[which(Filter_Summ_Candid$GI_TG_total<10)]=10
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=10)&(Filter_Summ_Candid$GI_TG_total<25))]=25
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=25)&(Filter_Summ_Candid$GI_TG_total<50))]=50
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=50)&(Filter_Summ_Candid$GI_TG_total<100))]=100
  Filter_Summ_Candid$TG_All_Categ[which((Filter_Summ_Candid$GI_TG_total>=100))]=200
  
  Filter_Summ_Candid$TG_All_Categ=as.factor(Filter_Summ_Candid$TG_All_Categ)
  
  # Change the point size
  nbaplot4_detailed=ggplot(Filter_Summ_Candid, aes(x=Raw_MutRate, y=Perce_Impact_Final,label="GI_ALL"))+
    geom_point(data=Filter_Summ_Candid,aes(size=TG_All_Categ,color=TG_All_Categ))+coord_cartesian(ylim = c(15, 100),xlim = c(0,5))+
    xlab("Mutation Rate of GI, Less than 5%") + ylab("Percentage of Impacted TG, (Outlier Genes)")+theme_bw()+theme(panel.grid = element_blank())+
    theme(legend.position="top")
  
  nbaplot4_detailed=nbaplot4_detailed + scale_y_continuous(trans = squish_trans(50, 100,10),
                                                           breaks = seq(15, 100, by =10))
  
  Labled_Plot_highimp_Det5=nbaplot4_detailed + geom_label_repel(aes(label = GI_ALL),
                                                                box.padding   = 0.1,
                                                                point.padding = 0.1,
                                                                segment.color = 'grey50')+ 
    labs(size="Number of Downstream Target Genes",colour="Number of Downstream Target Genes")
  
  file_Plot=paste0("./Results/Labeled_Candidate_Genes_5Percent_NewFilter_",CANCER_TYPE,"_RndTG.pdf")
  text.factor<-2
  dpi <- text.factor * 100
  width.calc <- 2300 / dpi
  height.calc <- 1500 / dpi
  ggsave(filename = file_Plot,  dpi = dpi, width = width.calc,  height = height.calc ,  units = 'in', plot = Labled_Plot_highimp_Det5)
}