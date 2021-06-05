library(sp)
library(doBy)
library(rgdal)
library(spdep)
library(maptools)
library(spatialreg) 


###################### Preparing raw data for the loop ######################

# loading company data
firms<-read.csv("dane/geoloc data.csv", header=TRUE, dec=",", sep=";")

# adding parameters – individual sectors
param<-data.frame(SEK_PKD7=c("A", "B", "C", "D" ,"E", "F", "G", "H", "I", "J", "K" ,"L", "M", "N", "O", "P", "Q", "R", "S"), SEK_agg=c("agri", "prod", "prod", "prod" ,"prod", "constr", "serv", "serv", "serv", "serv", "serv" ,"serv", "serv", "serv", "serv", "serv", "serv", "serv", "serv"), roa_ind=c(2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11))
# merging parameters and base dataset
firms1<-merge(firms, param, by="SEK_PKD7")
# adding parameters - aggregated sectors and their ROA
param2<-data.frame(SEK_agg=c("agri", "prod", "constr", "serv"), roa_sec=c(2,3.5,5,8))

# merging parameters and datasets with companies
firms1<-merge(firms1, param2, by="SEK_agg")

# premium for central localization (within close range to region's capital city)
firms1$roa_geo<-ifelse(firms1$powiatowe=="powiat Lublin", 1.5,0)

# calculating euclidean distance between a point and Lublin's centroid
coords<-as.matrix(data.frame(x=firms1$coords.x1, y=firms1$coords.x2)) 

core<-c(22.5666700, 51.2500000) # Lublin's center coordinates
firms1$dist<-spDistsN1(coords, core, longlat=TRUE)

# specifying roa_dist as a squared relation to the distance from the region's center
firms1$roa_dist<-0
firms1$roa_dist<- (1.5-(firms1$dist/100))^2

# assigning ROA mean as a sum with sectoral, geographical and distance-based roa
firms1$roa_param<-firms1$roa_sec + firms1$roa_geo + firms1$roa_dist 

# generating ROA with individualized mean
for(i in 1:dim(firms1)){
  firms1$roa[i]<-rnorm(1, firms1$roa_param[i], 0.045)}

# binary variables for aggregated sectors
firms1$agri<-ifelse(firms1$SEK_agg=="agri",1,0)
firms1$prod<-ifelse(firms1$SEK_agg=="prod",1,0)
firms1$constr<-ifelse(firms1$SEK_agg=="constr",1,0)
firms1$serv<-ifelse(firms1$SEK_agg=="serv",1,0)

# calculating euclinean distance between a point and Lublin's centroid
coords<-as.matrix(data.frame(x=firms1$coords.x1, y=firms1$coords.x2)) 
core<-c(22.5666700, 51.2500000) # Lublin's center coordinates
firms1$dist<-spDistsN1(coords, core, longlat=TRUE)

# randomly spread the points by epsilon (to avoid duplicates and issues in weighting matrix)
epsilon.x<-rnorm(dim(firms1)[1], mean=0, sd=0.015)
epsilon.y<-rnorm(dim(firms1)[1], mean=0, sd=0.015)

firms1$xxe<-firms1[,24]+epsilon.x
firms1$yye<-firms1[,25]+epsilon.y


##################### Ordering the sample #################### 
# random ordering of the sample 

set.seed(777)

firms1$los<-runif(dim(firms1)[1], 0,1)
firms1<-orderBy(~los, data=firms1)

######################### Parameters of the simulation #########################

# creating the structure of results matrix 
otherMat<-matrix(0,nrow=0,ncol=8)
colnames(otherMat)<-c('AIC', 'BIC', 'time', 'rho', 'knn', 'obsSample', 'minInterval', 'moran')

n.col<-99 # n.col+1 - number of max knn in the SWM
n.rowVec<-c(100,200,500,1000,10000) # sizes of subsamples

for(g in 2:5){
  
  n.row<-n.rowVec[g] # number of observations in subsample
  
  selector<-c(1:n.row) # initiating first selector assignment (which rows will be considered)
  
  while(max(selector)<=10000){ #limiting number of subsamples considered
    
    # partial results table to be filled during the loop iterations
    other.sdm<-matrix(0, nrow=n.col, ncol=8)
    
    eq<-roa~prod+constr+serv+dist # model structure
    
    # choosing the subsample
    datax<-firms1[selector,]
    
    # coordinates of points in the subsample (for SWM)
    crds<-as.matrix(datax[,c('xxe', 'yye')])
    
    #######################################

    # loop for model estimation with different knn in SWM
    for(i in 1:n.col){ # n.col sets number of iteration and knn considered
      
      pkt.knn<-knearneigh(crds, k=(i+1), longlat = NULL) 
      print(i)
      pkt.k.nb<-knn2nb(pkt.knn)
      pkt.k.sym.nb<-make.sym.nb(pkt.k.nb) # making SWM symmetric
      pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb)
      
      # Spatial Durbin Model estimation
      start.time <- Sys.time() # measuring execution time
      model.sdm<-lagsarlm(eq, data=datax, pkt.k.sym.listw, method="LU", type="mixed")
      end.time <- Sys.time()
      time.sdm<- difftime(end.time, start.time, units="secs")

      other.sdm[i,1]<-AIC(model.sdm) # AIC result
      other.sdm[i,2]<-BIC(model.sdm) # BIC result
      other.sdm[i,3]<-time.sdm
      other.sdm[i,4]<-model.sdm$rho
      other.sdm[i,5]<-i+1 # how many neighbours used in knn
      other.sdm[i,6]<-n.row #how many observations in the sample
      other.sdm[i,7]<-min(selector) #ID of first subsample observation (ordering of samples)
      other.sdm[i,8]<-moran.test(datax$roa, pkt.k.sym.listw)$estimate[1] # Moran statistic
    } 
    
    print(min(selector))
    otherMat<-rbind(otherMat,other.sdm)
    
    selector<-selector+n.row # increasing IDs in the selector to get consecutive IDs for the next sample
  }
}

write.csv2(otherMat, 'resultsEmpiricalSimulation.csv')