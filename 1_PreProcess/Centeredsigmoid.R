## Center Sigmoid function

Centeredsigmoid <-function(x){
  DD=as.matrix(x)
  c=1
  DD=( 1.0 - exp(-DD/c))/ ( 1.0 + exp(-DD/c))
  
  return(DD)
}