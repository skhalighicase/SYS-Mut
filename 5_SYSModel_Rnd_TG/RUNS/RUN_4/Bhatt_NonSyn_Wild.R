## Last Modification 6/6/2020
#### Calculates bhattacharyya.dist between NonSyn_Wild distributions

Bhatt_NonSyn_Wild<-function(NonSyn_Wild,Control){
  Measur_All=list()
  Measur_All$Bht_Nm_M<-(as.numeric(bhattacharyya.dist(NonSyn_Wild$Mu_Phi2[1],NonSyn_Wild$Mu_Phi2[2],as.matrix((NonSyn_Wild$Var_Phi2[1])^2),as.matrix((NonSyn_Wild$Var_Phi2[2])^2))))

  ####
  ## Measure Control
  Maxfold=length(Control)
  Contrl_Tmp_Bh=array(list(),dim=Maxfold)
  Contrl_Tmp_Bh<-foreach(Foldidx=1:Maxfold)%dopar%{
  # for(Foldidx in 1:Maxfold){
    library(fpc)
    Curr_Control=Control[[Foldidx]]
    Contrl_Tmp_Bh[[Foldidx]]<-(bhattacharyya.dist(Curr_Control$Mu_Phi2[1],Curr_Control$Mu_Phi2[2],as.matrix((Curr_Control$Var_Phi2[1])^2),as.matrix((Curr_Control$Var_Phi2[2])^2)))
  }
  Contr_Bha_samples=as.vector(unlist(Contrl_Tmp_Bh))
  Res_Wtest=wilcox.test(Contr_Bha_samples,NULL,alternative = "less",mu = as.numeric(Measur_All$Bht_Nm_M),conf.level = 0.95) 
  Measur_All$WPval=Res_Wtest$p.value
  Measur_All$Bht_C_Mean=mean(unlist(Contrl_Tmp_Bh))
  Measur_All$Bht_C_Var= var(unlist(Contrl_Tmp_Bh)) 
  
  return(Measur_All)
}
