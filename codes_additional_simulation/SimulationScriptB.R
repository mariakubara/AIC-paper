# Simulation of type B

# Uniform distribution, all from 'agriculture' sector
numObs<-c(200,500,1000)

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
    
    ################# Generating random points (centroids for clusters) ###################
    uniformSet <- runifpoint(numObs[n], win=W, giveup=1000, nsim=1, drop=TRUE, ex=NULL)

    data<-uniformSet
    data.cords <- coords(data)

    pointsMat.sp<-SpatialPoints(data.cords)
    proj4string(pointsMat.sp)<-CRS("+proj=merc +datum=NAD83") 
    plot(pointsMat.sp)
    
    ########## Check if points are withing the lubelskie window #####################
    
    # check which points are within powiat Lublin 
    # (geographical premium for being near the region's capital city)
    
    powLub <- inside.owin(pointsMat.sp@coords[,1], pointsMat.sp@coords[,2], W.lub)
    
    # change projection and create Spatial Points dataset
    pointsMat.sp<-spTransform(pointsMat.sp, CRS("+proj=longlat +datum=NAD83")) #sferyczne
    
    # creating a data frame with simulated companies characteristics
    data1<-cbind(pointsMat.sp@coords, powLub)
    data1<-as.data.frame(data1)
    
    # geographical premium
    data1$roa_geo<-ifelse(data1$powLub==1, 1.5,0)
    
    fulldata<-as.data.frame(data1)
    fulldata$sector<-'agri'

    fulldata$roa_sec<-NA
    fulldata[fulldata$sector=='agri',]$roa_sec<-2
    
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
    write.table(resultsMat, file="results_additional/simulationUniform.csv", row.names = F, col.names = T)
    
  }
}