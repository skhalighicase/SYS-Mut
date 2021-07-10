## ## Read Influence Network
Fun_InfluenceNet<-function(){

Influ_net<-read.delim("./InputData/SYS_Net_HGNC.txt")

save(Influ_net,file="./Data/SYS_Net.RData")
}
