## Softtheshold function

SoftThresh<-function(x){
  DD=as.matrix(x)
  c= 1
  GG = sign(DD)*((sqrt(DD^2 + c^2))-c)
}
