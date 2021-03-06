#######################################################################
#loading packages

library(spdep)
library(sp)
library(maptools)
library(rgdal)
library(spatialreg)
library(devtools)
library(spatstat)
library(sf)
library(tidyverse)
library(osmdata)

#######################################################################
################ Preparing window ########################

woj<-readOGR("data","wojewodztwa") 
region<-woj[woj$jpt_nazwa_=="lubelskie",] 
region.plot<-region
region<-spTransform(region, CRS("+proj=merc +datum=WGS84")) 

W<-as(region, "owin") # converting to owin class (necessary for point pattern operations)
class(W)


################ Window with Lublin poviat (central location premium) #######
pow<-readOGR("data", "powiaty") 
pow<- spTransform(pow, CRS("+proj=longlat +datum=NAD83"))

region.lub<-pow[pow$jpt_nazwa_=="powiat Lublin",] 
region.lub<-spTransform(region.lub, CRS("+proj=merc +datum=WGS84")) 

W.lub<-as(region.lub, "owin") # conversion to owin class 
class(W.lub)

################ Preparation of lubelskie voivodeship borders in sf class for visualisations ###
woj.sf <- st_read("data/wojewodztwa.shp")
woj.sf <- st_transform(woj.sf, crs = "+proj=longlat +datum=NAD83")

# limit to lubelskie
lub.woj.sf <- woj.sf %>% filter(jpt_nazwa_=='lubelskie') 
