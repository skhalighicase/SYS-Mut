
TrasnformCNA_Meth <-function(){
  File_Name=paste0("./Data/CommonGenesamples_",CANCER_TYPE,".RData")
  load (File_Name)
  
  ## Apply Sigmoid function on CNV and methylation data
  
  source("Centeredsigmoid.R")
  SigmoidCNV<-list(
    SigmoidCNV0 = Centeredsigmoid(Comm_Sam_CNV0),
    SigmoidCNV1 = Centeredsigmoid(Comm_Sam_CNV1)
  )
  SigmoidMeth<-list(
    SigmoidMeth_Max0 = Centeredsigmoid(Comm_Sam_Meth_Max0),
    SigmoidMeth_Max1 = Centeredsigmoid(Comm_Sam_Meth_Max1),
    
    SigmoidMeth_Min0 = Centeredsigmoid(Comm_Sam_Meth_Min0),
    SigmoidMeth_Min1 = Centeredsigmoid(Comm_Sam_Meth_Min1),
    
    SigmoidMeth_Mean0 = Centeredsigmoid(Comm_Sam_Meth_Mean0),
    SigmoidMeth_Mean1 = Centeredsigmoid(Comm_Sam_Meth_Mean1)
  )
  ## Apply soft-margin function on CNV and methylation data
  
  source("SoftThresh.R")
  ThreshCNV<-list(
    ThreshCNV0 = SoftThresh(Comm_Sam_CNV0),
    ThreshCNV1 = SoftThresh(Comm_Sam_CNV1)
  )
  ThreshMeth<-list(
    ThreshMeth_Max0 = SoftThresh(Comm_Sam_Meth_Max0),
    ThreshMeth_Max1 = SoftThresh(Comm_Sam_Meth_Max1),
    
    ThreshMeth_Min0 = SoftThresh(Comm_Sam_Meth_Min0),
    ThreshMeth_Min1 = SoftThresh(Comm_Sam_Meth_Min1),
    
    ThreshMeth_Mean0 = SoftThresh(Comm_Sam_Meth_Mean0),
    ThreshMeth_Mean1 = SoftThresh(Comm_Sam_Meth_Mean1)
  )
  
  FileName_Save=paste0("./Data/Transformed_CNV_Meth_",CANCER_TYPE,".RData")
  save(SigmoidCNV,SigmoidMeth,ThreshCNV,ThreshMeth,file=FileName_Save)
  
}
