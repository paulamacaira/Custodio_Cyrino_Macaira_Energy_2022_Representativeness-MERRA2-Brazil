if(!require(dplyr)) {install.packages("dplyr"); library(dplyr) }
if(!require(lubridate)) {install.packages("lubridate"); library(lubridate) }
if(!require(tibble)) {install.packages("tibble"); library(tibble) }
if(!require(feather)) {install.packages("feather"); library(feather) }
if(!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse) }
if(!require(Metrics)) {install.packages("Metrics"); library(Metrics) }
if(!require(reshape2)) {install.packages("reshape2"); library(reshape2) }
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2) }
if(!require(BBmisc)) {install.packages("BBmisc"); library(BBmisc) }
if(!require(readxl)) {install.packages("readxl"); library(readxl) }
if(!require(hash)) {install.packages("hash"); library(hash) }
if(!require(gtools)) {install.packages("gtools"); library(gtools) }
if(!require(plotly)) {install.packages("plotly"); library(plotly) }
if(!require(raster)) {install.packages("raster"); library(raster) }
if(!require(rgdal)) {install.packages("rgdal"); library(rgdal) }
if(!require(ncdf4)) {install.packages("ncdf4"); library(ncdf4) }
if(!require(httr)) {install.packages("httr"); library(httr) }
if(!require(parallel)) {install.packages("parallel"); library(parallel) }
if(!require(forecast)) {install.packages("forecast"); library(forecast) }
if(!require(tseries)) {install.packages("tseries"); library(tseries) }
if(!require(fitdistrplus)) {install.packages("fitdistrplus"); library(fitdistrplus) }
if(!require(zoo)) {install.packages("zoo"); library(zoo) }
if(!require(cowplot)) {install.packages("cowplot"); library(cowplot) }
if(!require(numbers)) {install.packages("numbers"); library(numbers) }
if(!require(openxlsx)) {install.packages("openxlsx"); library(openxlsx) }
if(!require(tsibble)) {install.packages("tsibble"); library(tsibble) }
if(!require(fpp3)) {install.packages("fpp3"); library(fpp3) }

# login for downloading merra data
## Logar no site do MERRA
name <- "XXXX"
password <- "XXXX"

# directory for source code
## local onde esta arquivo do R: codigo principal, funcoes gerais, funcoes MERRA2, ...
dirsource <- "C:/ ..."

# base MERRA data directory
## Local Onde salvar pasta com arquivos baixados direto do site do MERRA2
dirmerrabase <- "C:/ ..."

# directory where MERRA data per point are stored
## local onde esta dados do MERRA2 por ponto de longitude e latitude
dirmerra <- "C: ..."

# load script for handling merra data
## Carregando funcoes especificas para carregar dados do MERRA2
source(paste0(dirsource,"/MERRA_data.R"))

##########################################################################################
##### DOWNLOAD OF MERRA DATA #############################################################
##########################################################################################

# the boundary of the box to download
## Intervalos de longitude e latitude que baixaremos os dados MERRA2
lon1<- -55.0
lat1<- -25.0
lon2<- -30.0
lat2<- -5.0

# define time span for download
## Perido que baixaremos os dados
date_seq<-seq(as.POSIXct("2008-01-01",tz="UTC"),as.POSIXct("2017-12-31",tz="UTC"),by="d")

# download
# some files may not be downloaded correctly (they are smaller)
# if this happens, delete them and repeat the downlaod

getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("U10M"),
                name,
                password,
                FALSE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("U50M"),
                name,
                password,
                FALSE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("V10M"),
                name,
                password,
                FALSE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("V50M"),
                name,
                password,
                FALSE)
getMERRADataBox(lon1,lat1,lon2,lat2,
                date_seq,c("DISPH"),
                name,
                password,
                FALSE)

# split date sequence due to memory restrictions
date_seq<-list(seq(as.POSIXct("2008-01-01",tz="UTC"),as.POSIXct("2017-12-31",tz="UTC"),by="d"))
setwd(dirmerrabase)

# convert to feather format
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"U10M","U10m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"U50M","U50m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"V10M","V10m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"V50M","V50m")
lapply(date_seq,convertMerraFeather,lon1,lat1,lon2,lat2,"DISPH","disph")

lonlat<-read_feather(paste(paste("./feather/LonLat","U10M",lon1,lat1,lon2,lat2,format(date_seq[[1]][1],"%Y%m%d"),format(date_seq[[1]][length(date_seq[[1]])],"%Y%m%d"),sep="_"),"/lonlat.feather",sep=""))
names(lonlat) <- c("long","lat")
write_feather(lonlat,paste(dirmerra,"/lonlat.feather",sep=""))

MerraDate <- seq(date_seq[[1]][1],date_seq[[length(date_seq)]][length(date_seq[[length(date_seq)]])],by="h")
hours <- as.POSIXct(rep(MerraDate[length(MerraDate)],23),tz="UTC")
hours <- hours + (1:23)*3600
MerraDate <- c(MerraDate,hours)
write_feather(as.data.frame(MerraDate),paste(dirmerra,"/MerraDate.feather",sep=""))

MerraDate <- read_feather(paste(dirmerra,"/MerraDate.feather",sep=""))
LonLat <- read_feather(paste(dirmerra,"/lonlat.feather",sep=""))
lonlat <- as.data.frame(LonLat)

pnames <- c("U10M","U50M","V10M","V50M","DISPH")

# change format from daily to point-wise files
# split into chunks due to memory restrictions
divs <- divisors(length(lonlat[,1])) 
binsize <- max(divs[which(divs<400)])
lll <- split(lonlat,f=rep(1:(length(lonlat[,1])/binsize),each=binsize))
for(ll in lll){
  invisible(apply(ll,1,saveMerraPoint,pnames,lon1,lat1,lon2,lat2,date_seq))
}


##########################################################################################
##### interpolation (nearest neighbor) ###################################################
##########################################################################################

# Data input from the wind farm administrator
# geographic coordinates
# turbine rotor height

dataset_wf = read_excel("Turbines_Wind_farms.xlsx", sheet =1, col_names = TRUE, skip = 0)

#confidential file "Turbines_Wind_farms.xlsx"
#dataset_wf
# Turbine # Latitude # Longitude # height # distance # lat_MERRA2 # lon_MERRA2

#Identifying the nearest MERRA-2 data point to the turbine
## Identificacao do ponto de dado MERRA-2 mais proximo do parque eolico

num_turbine = 24 #nrow(dataset_wf)
for (ind in 1:num_turbine) {
  
  Melhor_lat = 0
  Melhor_lon = 0
  Menor_d = Inf
  lat <- dataset_wf$Latitude[ind]
  long <- dataset_wf$Longitude[ind]
  
  for (i in 1:length(lonlat$long)) {
    a = lonlat$lat[i]
    b = lonlat$long[i]
    #Artigo [4]:
    rad <- pi/180
    d <- 6378.388*acos(sin(rad*lat) * sin(rad*a) + cos(rad*lat) * cos(rad*a) * cos(rad*b-rad*long))
    #d = sqrt(((dataset_wf$Latitude[ind] - a)^2)+((dataset_wf$Longitude[ind] - b)^2))
    if(d <= Menor_d ){
      Melhor_lat = a
      Melhor_lon = b
      Menor_d = d
    }
  }
  
  dataset_wf$distance = Menor_d
  dataset_wf$lat_MERRA2 = Melhor_lat
  dataset_wf$lon_MERRA2 = Melhor_lon
  
}


##########################################################################################
##### extrapolation (Hellman power law) ##################################################
##########################################################################################
ind = 1
Melhor_lat = dataset_wf$lat_MERRA2[ind]
Melhor_lon = dataset_wf$lon_MERRA2[ind]

Endereco = paste0(dirmerra,"/",Melhor_lon,"_",Melhor_lat,".feather")
base_hora <- read_feather(Endereco)

veloc_50M_hora <- sqrt(((base_hora$U50M)^2)+((base_hora$V50M)^2))
veloc_10M_hora <- sqrt(((base_hora$U10M)^2)+((base_hora$V10M)^2))

base_h = data.frame(veloc_50M=veloc_50M_hora,
                    veloc_10M=veloc_10M_hora, 
                    DISPH= base_hora$DISPH)

#validation
analise_dados = data.frame(var=colnames(base_h),Menor_zero =NA, Maior_25 =NA, Faltante=NA)
for (i in 1:ncol(base_h)) {
  c = base_h[is.na(base_h[,i]),i]
  analise_dados[3,(1+i)]= length(c)
  a = base_h[which(base_h[,i]<0),i]
  base_h[which(base_h[,i]<0),i] = NA
  analise_dados[1,(1+i)]= length(a)
  b = base_h[which(base_h[,i]>25),i]
  base_h[which(base_h[,i]>25),i] = NA
  analise_dados[2,(1+i)]= length(b)
}

#alpha friction coefficient
alpha <- (log(veloc_50M_hora)-log(veloc_10M_hora))/(log(50)-log(10+base_hora$DISPH))
#wind speed with power law
veloc_ext_hora <- veloc_50M_hora*((dataset_wf$height[ind]/50)^alpha)


###### Dados Diarios #######
dias = seq(as.POSIXct("2008-01-01",tz="UTC"),as.POSIXct("2017-12-31",tz="UTC"),by="d")

tam = length(base_hora$U10M)/24
veloc_50M_dia = c(data = NA, dim = tam)
veloc_10M_dia = c(data = NA, dim = tam)
DISPH_dia = c(data = NA, dim = tam)
for (i in 1:tam) {
  inicio = 1+((i-1)*24)
  fim = i*24
  veloc_50M_dia[i] = mean(veloc_50M_hora[inicio:fim]) 
  veloc_10M_dia[i] = mean(veloc_10M_hora[inicio:fim]) 
  DISPH_dia[i] = mean(base_hora$DISPH[inicio:fim]) 
}

alpha <- (log(veloc_50M_dia)-log(veloc_10M_dia))/(log(50)-log(10+DISPH_dia))
veloc_ext_dia <- veloc_50M_dia*((dataset_wf$height[ind]/50)^alpha)

base_d = data.frame(veloc_50M=veloc_50M_dia,
                    veloc_10M=veloc_10M_dia, 
                    DISPH= DISPH_dia)


###### Dados MESEs #######
inicio_data=as.Date("2008-01-01")
fim_data=as.Date("2017-12-31")
meses =length(seq(inicio_data,fim_data , by = "months"))

veloc_50M_mes = c(data = NA, dim = meses)
veloc_10M_mes = c(data = NA, dim = meses)
DISPH_mes = c(data = NA, dim = meses)

teste = data.frame(data = MerraDate, veloc_50M_hora = veloc_50M_hora,
                   veloc_10M_hora = veloc_10M_hora, DISPH_hora = base_hora$DISPH,
                   mes = month(MerraDate$MerraDate),ano = year(MerraDate$MerraDate))

cont = 1
for (i in year(inicio_data):year(fim_data)) {
  teste2 = filter(teste, ano == i)
  for (j in 1:12) {
    dados = filter(teste2, mes == j)
    veloc_50M_mes[cont] = mean(dados$veloc_50M_hora)
    veloc_10M_mes[cont] = mean(dados$veloc_10M_hora)
    DISPH_mes[cont] = mean(dados$DISPH_hora)
    cont = cont +1
  }
}

alpha <- (log(veloc_50M_mes)-log(veloc_10M_mes))/(log(50)-log(10+DISPH_mes))
veloc_ext_mes <- veloc_50M_mes*((dataset_wf$height[ind]/50)^alpha)

base_m = data.frame(veloc_50M=veloc_50M_mes,
                    veloc_10M=veloc_10M_mes, 
                    DISPH= DISPH_mes)

