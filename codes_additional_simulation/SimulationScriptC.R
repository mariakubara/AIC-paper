# Simulation of type C

# Uniform 'agriculture' 50%, clusters per 25 observations 'service' 50%

numObs<-list(c(100,4,25),c(250, 10,25),c(500, 20,25))


# creating the structure of results matrix
resultsMat<-matrix(0, nrow=1, ncol=8) 
colnames(resultsMat)<-c('AIC', 'BIC', 'time', 'rho', 'knn', 'obsSample', 'seed', 'moran')

# vector with seeds that will be used to ensure replicability of a simulation
seed<-seq(1,1000,50) 

# three sizes of data (200, 500 and 1000) with 20 replication each
for(n in 1:3){
  for(s in 1:20){
    
    # setting the seed for reproducibility 
    set.seed(seed[s])
    
    ################# Generating random points from uniform distribution ###################
    uniform6 <- runifpoint(numObs[[n]][1], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
    
    #extracting the points 
    unif.cords <- coords(uniform6)
    unif.cords.sp<-SpatialPoints(unif.cords)
    
    #setting the projection aligned with spatstat - plenar projection
    proj4string(unif.cords.sp)<-CRS("+proj=merc +datum=NAD83") # plenarne
    
    # check which points are within powiat Lublin 
    powLub <- inside.owin(unif.cords.sp@coords[,1],unif.cords.sp@coords[,2], W.lub)

    # change projection and create Spatial Points dataset
    unif.cords.sp<-spTransform(unif.cords.sp, CRS("+proj=longlat +datum=NAD83")) #spherical
    
    # creating a data frame with simulated companies characteristics
    data1<-cbind(unif.cords.sp@coords, powLub)
    data1<-as.data.frame(data1)
    
    data1$roa_geo<-ifelse(data1$powLub==1, 1.5,0)
    data1$subsample <- 1
    
    names(data1)<- c("xi", "yi", "powLub", "roa_geo", "subsample")
    
    ################# Generating clusterised subsample ###################
    data2 <- runifpoint(numObs[[n]][2], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)
    
    data2.cords <- coords(data2)
    data2.cords.sp<-SpatialPoints(data2.cords)
    
    #setting the projection aligned with spatstat - plenar projection
    proj4string(data2.cords.sp)<-CRS("+proj=merc +datum=NAD83") 
    
    centrum2<- as.matrix(cbind(data2.cords.sp@coords[,1], data2.cords.sp@coords[,2]))
    
    # draw points around the centroid in a given radius 
    pointsMat2<-matrix(NA, nrow=0, ncol=2)
    for(i in 1:nrow(centrum2)){
      punkt<-centrum2[i,]
      ileClust<-numObs[[n]][3]
      
      rmaximum <- 10000
      angle <- runif(ileClust)*3.14*2
      radius <- runif(ileClust)*rmaximum
      
      xi <- punkt[1]+radius*cos(angle)
      yi <- punkt[2]+radius*sin(angle)
      
      pointsMat2<-rbind(pointsMat2, cbind(xi,yi))
      
    }
    
    pointsMat2.sp<-SpatialPoints(pointsMat2)
    #setting the projection aligned with spatstat - plenar projection
    proj4string(pointsMat2.sp)<-CRS("+proj=merc +datum=NAD83") 
    
    # check which points are within powiat Lublin 
    powLub <- inside.owin(pointsMat2.sp@coords[,1], pointsMat2.sp@coords[,2], W.lub)
    
    # change projection and create Spatial Points dataset
    pointsMat2.sp<-spTransform(pointsMat2.sp, CRS("+proj=longlat +datum=NAD83")) #sferyczne
    
    # creating a data frame with simulated companies characteristics
    data2<-cbind(pointsMat2.sp@coords, powLub)
    data2<-as.data.frame(data2)
    
    data2$roa_geo<-ifelse(data2$powLub==1, 1.5,0)
    data2$subsample<-4
    
    #################################################################
    ############ Merging subsamples ############################
    
    fulldata<-rbind(data1, data2)
    fulldata<-as.data.frame(fulldata)
    fulldata$sector<-NA
    fulldata[fulldata$subsample==1,]$sector<-'agri'
    #fulldata[fulldata$subsample==2,]$sector<-'constr'
    fulldata[fulldata$subsample==4,]$sector<-'serv'
    
    fulldata$roa_sec<-NA
    fulldata[fulldata$sector=='agri',]$roa_sec<-2
    #fulldata[fulldata$sector=='prod',]$roa_sec<-3.5
    #fulldata[fulldata$sector=='constr',]$roa_sec<-5
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
    other.sdm<-matrix(0, nrow=n.knn, ncol=8) 
    colnames(other.sdm)<-c('AIC', 'BIC', 'time', 'rho', 'knn', 'obsSample', 'seed', 'moran')
    
    
    for(i in 1:n.knn){ 
      
      # i+1 -> how many neighbours in knn
      # i -> in which repetition are we in
      
      crds<-as.matrix(cbind(datax$xi, datax$yi))
      
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
      
      # saving partial results
      other.sdm[i,1]<-AIC(model.sdm) # AIC.sdm
      other.sdm[i,2]<-BIC(model.sdm) # BIC.sdm
      other.sdm[i,3]<-time.sdm
      other.sdm[i,4]<-model.sdm$rho
      other.sdm[i,5]<-i+1 # how many neighbours used in knn
      other.sdm[i,6]<-nrow(datax) # how many observations in a sample
      other.sdm[i,7]<-seed[s] # seed used for simulation
      other.sdm[i,8]<-moran.test(datax$roa, pkt.k.sym.listw)$estimate[1] # Moran statistic
    }
    
    # merging partial results with general results table
    resultsMat<-rbind(resultsMat, other.sdm)
    
    # saving the results
    write.table(resultsMat, file="results_additional/simulationUniformCluster25.csv", row.names = F, col.names = T)
    
  }
}