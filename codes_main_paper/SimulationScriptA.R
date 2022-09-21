# Simulation of type A

# Clusters with 20 observations each, all from 'service' sector

numObs<-list(c(10,20),c(25,20),c(50,20))

############## Set with concentrations, 20 points in each cluster ######################

# creating the structure of results matrix
resultsMat<-matrix(0, nrow=1, ncol=14) 
colnames(resultsMat)<-c('AIC', 'BIC', 'time', 'rho', 'knn', 'obsSample', 'seed', 'moran', 
                        "direct", "indirect", "total", "AICsem", "BICsem", "lambda")

# vector with seeds that will be used to ensure replicability of a simulation
seed<-seq(1,1000,50) 

# three sizes of data (200, 500 and 1000) with 20 replication each
for(n in 1:3){
  for(s in 1:20){
    
    # setting the seed for reproducibility 
    set.seed(seed[s])
    
    ################# Generating random points (centroids for clusters) ###################
    uniformSet <- runifpoint(numObs[[n]][1], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
    
    data<-uniformSet
    data.cords <- coords(data)
    data.cords.sp<-SpatialPoints(data.cords)
    
    #setting the projection aligned with spatstat - plenar projection
    proj4string(data.cords.sp)<-CRS("+proj=merc +datum=NAD83") 
    
    centroid<- as.matrix(cbind(data.cords.sp@coords[,1], data.cords.sp@coords[,2]))
    
    
    #2 draw points around the centroid in a given radius 
    pointsMat<-matrix(NA, nrow=0, ncol=2)
    
    for(i in 1:nrow(centroid)){
      point<-centroid[i,]
      howMany<-numObs[[n]][2]
      
      rmaximum <- 10000
      angle <- runif(howMany)*3.14*2
      radius <- runif(howMany)*rmaximum
      
      xi <- point[1]+radius*cos(angle)
      yi <- point[2]+radius*sin(angle)
      
      pointsMat<-rbind(pointsMat, cbind(xi,yi))
      
    }
    
    pointsMat.sp<-SpatialPoints(pointsMat)
    proj4string(pointsMat.sp)<-CRS("+proj=merc +datum=NAD83") 
    #plot(pointsMat.sp)
    
    ########## Check if points are withing the lubelskie window #####################
    
    # check which points are within powiat Lublin 
    # (geographical premium for being near the region's capital city)

    powLub <- inside.owin(pointsMat.sp@coords[,1], pointsMat.sp@coords[,2], W.lub)
    
    # change projection and create Spatial Points dataset
    pointsMat.sp<-spTransform(pointsMat.sp, CRS("+proj=longlat +datum=NAD83"))
    
    # creating a data frame with simulated companies characteristics
    data1<-cbind(pointsMat.sp@coords, powLub)
    data1<-as.data.frame(data1)
    
    # geographical premium
    data1$roa_geo<-ifelse(data1$powLub==1, 1.5,0)
    
    fulldata<-as.data.frame(data1)
    fulldata$sector<-'serv'

    fulldata$roa_sec<-NA
    fulldata[fulldata$sector=='serv',]$roa_sec<-8
    
    library(spdep)
    # calculating euclinean distance between a point and Lublin's centroid
    coords<-as.matrix(data.frame(x=fulldata$x, y=fulldata$y)) 
    core<-c(22.5666700, 51.2500000) # Lublin's center coordinates
    fulldata$dist<-spDistsN1(coords, core, longlat=TRUE)
    
    # specifying roa_dist as a squared relation to the distance from the region's center
    fulldata$roa_dist<-0
    fulldata$roa_dist<- (1.5-(fulldata$dist/100))^2
    
    # adding all elements to create unique mean utilizable for generating the ROA value
    fulldata$roa_param<-fulldata$roa_sec+fulldata$roa_geo+fulldata$roa_dist 
    
    # generating ROA with individualized mean
    for(i in 1:nrow(fulldata)){
      fulldata$roa[i]<-rnorm(1, fulldata$roa_param[i], 0.045)}
    
    # binary values for sector recognition
    fulldata$agri<-ifelse(fulldata$sector=="agri",1,0)
    fulldata$prod<-ifelse(fulldata$sector=="prod",1,0)
    fulldata$constr<-ifelse(fulldata$sector=="constr",1,0)
    fulldata$serv<-ifelse(fulldata$sector=="serv",1,0)
    
    # MODELLING
    eq<-roa~prod+constr+serv+dist # model equation
    datax<-fulldata
    
    n.knn<-49 # range of knn [2:50], n.knn+1
    
    # partial results table to be filled during the loop iterations
    other.sdm<-matrix(0, nrow=n.knn, ncol=14) 
    colnames(other.sdm)<-c('AIC', 'BIC', 'time', 'rho', 'knn', 'obsSample', 'seed', 'moran', 
                           "direct", "indirect", "total", "AICsem", "BICsem", "lambda")
    
    
    for(i in 1:n.knn){ 
      
      # i+1 -> how many neighbours in knn
      # i -> in which repetition are we in
      
      crds<-as.matrix(cbind(datax$x, datax$y))
      
      pkt.knn<-knearneigh(crds, k=(i+1), longlat = NULL) 
      print(i) # progress check
      
      pkt.k.nb<-knn2nb(pkt.knn)
      pkt.k.sym.nb<-make.sym.nb(pkt.k.nb) # making the W matrix symmetric 
      pkt.k.sym.listw<-nb2listw(pkt.k.sym.nb) # changing the class to listw 
      
      # Spatial Durbin Model estimation
      start.time <- Sys.time() # measuring system time
      model.sdm<-lagsarlm(eq, data=datax, pkt.k.sym.listw, method="LU", type="mixed", zero.policy=TRUE)
      end.time <- Sys.time()
      time.sdm<- difftime(end.time, start.time, units="secs")
      
      #calculating impacts
      W.c <- as(as_dgRMatrix_listw(pkt.k.sym.listw), "CsparseMatrix")
      trMat <- trW(W.c, type = "mult")
      model.SDM.imp <- impacts(model.sdm, tr = trMat, R = 2000)
      summary(model.SDM.imp, zstats = TRUE, short = TRUE)
      other.sdm[i, 9]<-model.SDM.imp$res$direct # direct
      other.sdm[i,10]<-model.SDM.imp$res$indirect # indirect
      other.sdm[i,11]<-model.SDM.imp$res$total # total
      
      # saving partial results
      other.sdm[i,1]<-AIC(model.sdm) # AIC.sdm
      other.sdm[i,2]<-BIC(model.sdm) # BIC.sdm
      other.sdm[i,3]<-time.sdm
      other.sdm[i,4]<-model.sdm$rho
      other.sdm[i,5]<-i+1 # how many neighbours used in knn
      other.sdm[i,6]<-nrow(datax) # how many observations in a sample
      other.sdm[i,7]<-seed[s] # seed used for simulation
      other.sdm[i,8]<-moran.test(datax$roa, pkt.k.sym.listw)$estimate[1] # Moran statistic
      
      #making SEM model
      model.sem<-errorsarlm(eq, data=datax, pkt.k.sym.listw,zero.policy=TRUE)
      other.sdm[i,12] <- AIC(model.sem) 
      other.sdm[i,13] <- BIC(model.sem) 
      other.sdm[i,14] <- model.sem$lambda
    }
    
    # merging partial results with general results table
    resultsMat<-rbind(resultsMat, other.sdm)
    
    # saving the results
    write.table(resultsMat, file="results_simulations/simulationCluster.csv", row.names = F, col.names = T)
    
  }
}