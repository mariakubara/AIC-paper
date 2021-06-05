
gamma <- function(dataVect){
  n <- length(dataVect)
  gammaMat <- matrix(0, n, n)
  
  for(i in 2:n){
    for(j in 1:(i-1)){
      gamma <- 0.5 * (dataVect[j]-dataVect[i])*(dataVect[j]-dataVect[i])
      gammaMat[i,j] <- gamma
      gammaMat[j,i] <- gamma
    }
  }
  
  return(gammaMat)
}


semiVar <- function(allGamma, binary){
  gammaMatrix <- allGamma * binary
  
  nonEmpty <- sum(colSums(gammaMatrix != 0))
  nonEmptyPairs <- nonEmpty/2
  
  sumGamma <- sum(colSums(gammaMatrix))
  
  semi <- sumGamma/nonEmptyPairs
  
  return(semi)
}

points<-read.csv("C:\\Users\\maria\\Desktop\\teaching\\data\\support_centres.csv")
points
crds<-cbind(points$lat,points$long)
crds

i<-5
library(spdep)
pkt.knn<-knearneigh(crds, k=i, longlat = NULL) # układ planarny, knn=zalezne

pkt.k.nb<-knn2nb(pkt.knn)
pkt.knn.mat <-nb2mat(pkt.k.nb)
binary <- pkt.knn.mat*i
binary
pkt.k.sym.nb<-make.sym.nb(pkt.k.nb) # usymetrycznienie macierzy
pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)

nrow(crds)
z<- runif(37)*20

semiVarKnn <- function(crds,z,maxKnn){

  #do gamma for z
  allGamma <- gamma(z)
  
  semiKnn <- rep(0,maxKnn)
  
  #iterate through knn
  for (i in 2:maxKnn){
    
    pkt.knn<-knearneigh(crds, k=i, longlat = NULL) # układ planarny, knn=zalezne
    
    pkt.k.nb<-knn2nb(pkt.knn)
    pkt.knn.mat <-nb2mat(pkt.k.nb)
    binary <- pkt.knn.mat*i
    binary
    
    semiK <- semiVar(allGamma,binary)
    
    semiKnn[i]<-semiK
  
  }
  
  plot(2:maxKnn, semiKnn[2:maxKnn])
  
  return(semiKnn)
}

g<- semiVarKnn(crds, z, 10)