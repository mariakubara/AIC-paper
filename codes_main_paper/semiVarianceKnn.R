
# Set of functions which are necessary to generate 
# semiVariance results for incremental knn

# gamma and semiVar are support functions, while semiVarKnn
# is the main function which allows to generate semiVariance results for
# data structure prepared in a following form:
# coordinates (x,y), dependent variable (assigned to each coordinate pair),
# and maxKnn which sets the range for semiVariance testing (k from 2 to maxKnn)

# semiVarKnn returns a vector in which semiVariance values are provided for
# each knn, knn=1 is set to zero and actualy calculated values start with knn=2 
# up to the last element of the vector which shows semiVariance for knn=maxKnn

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


semiVarKnn <- function(crds,z,maxKnn){

  #do gamma for z values (dependent variable)
  allGamma <- gamma(z)
  
  semiKnn <- rep(0,maxKnn)
  
  #iterate through knn
  for (i in 2:maxKnn){
    
    pkt.knn<-knearneigh(crds, k=i, longlat = NULL) # planar coordinates, knn dependent
    
    pkt.k.nb<-knn2nb(pkt.knn)
    pkt.knn.mat <-nb2mat(pkt.k.nb)
    binary <- pkt.knn.mat*i
    #binary
    
    semiK <- semiVar(allGamma,binary)
    
    semiKnn[i]<-semiK
  
  }
  
  #plot(2:maxKnn, semiKnn[2:maxKnn])
  
  return(semiKnn)
}
