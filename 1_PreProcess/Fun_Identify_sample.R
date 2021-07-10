## Identify the samples (Tumor Samples and Normal Samples)
Fun_Identify_sample <-function(TCGA_Barcode){

  PhenoType = -1 #Unknown
  p = grepl("TCGA", TCGA_Barcode)
  if (p){
    TCGA_Barcoderef=substr(TCGA_Barcode,1,nchar(TCGA_Barcode))
    if (!identical(substr(TCGA_Barcoderef,1,16),substr(TCGA_Barcode,1,min(16,nchar(TCGA_Barcode))))){
      TCGA_Barcode = "UnKnown"
    }
  }
  if (!(identical(TCGA_Barcode, "UnKnown"))){
    TCGA_Barcode = TCGA_Barcoderef
  }
  ShortBarcode=substr(TCGA_Barcode, 1, 12)
  TcgaType = substr(TCGA_Barcode, 13, 15)
  
  if (identical(TcgaType, '-01')){  ###01A 02A 01B, ....
    PhenoType=1  # Tumor sample
  }else if (identical(TcgaType, '-11')){
    PhenoType=0  # Normal sample
  }
  return (PhenoType)
}