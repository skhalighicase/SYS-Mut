## Lat Update 6/9/2020

Thresh_Outlier<-function(Res_temp){

  ############
  Thresho=10 #we exclude mutation rates greater than 10%
  
  Target_RES=Res_temp[which(Res_temp$Percentage_NonSy_Wild<Thresho),]  ## we exclude mutation rates greater than 10%
  Rn_TG=dplyr::select(Target_RES,Percentage_NonSy_Wild,Perce_Impact_Final)
  Rn_TG_IDX=which(Rn_TG$Perce_Impact_Final>0)
  Rn_TG_U=Rn_TG[Rn_TG_IDX,]
  
  # Fit a lognormal distribution
  fit_params <- fitdistr(Rn_TG_U[,2],"lognormal")
  
  ##############################
  ## Q-Q plot
  par(mar=c(1,1,1,1))
  # create a vector of quantiles
  quants <-seq(0,1,length=81)[2:80]
  
  # find quantiles for the fitted distribution
  fit_quants <- qlnorm(quants,fit_params$estimate['meanlog'], fit_params$estimate['sdlog'])
  
  # find quantiles of the original data
  data_quants <- quantile(Rn_TG_U[,2],quants)
  
  # fit and data quantiles side by side
  Fited_and_Data_quant<-data.frame(fit_quants,data_quants)
  
  # create Q-Q plot
  plot(fit_quants, data_quants, xlab="Theoretical Quantiles", ylab="Sample Quantiles")
  title(main = "Q-Q plot of lognormal fit against data")
  abline(0,1)
  
  ###############################
  return(Fited_and_Data_quant)
}