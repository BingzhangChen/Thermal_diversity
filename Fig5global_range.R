#File of SST temperature
#Downloaded from https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA+High-resolution+Blended+Analysis&Variable=Sea+Surface+Temperature&group=0&submit=Search on 22 Dec 2018
setwd('~/OneDrive/Kremer/ideal')
rm(list=ls())
source('~/OneDrive/FlexEFT1D/Rscripts/Interpolate_WOA.R')
library(plot3D)
library(MASS)
library(fields)
tfile   <- "sst.day.mean.ltm.1971-2000.nc"
MLDfile <- "mldpd.mnltm.nc"
lon     <- ncread(tfile, 'lon')
lat     <- ncread(tfile, 'lat')
sst     <- ncread(tfile, 'sst')
NT      <- dim(sst)[3]

#The following function extracts the daily SST from a given station
stn_sst <- function(LON, LAT){
  if (LON < 0){
    LON <- LON + 360
  }else if (LON > 360){
    stop("Longitude incorrect!")
  }
  latlon <- list(x=LON, y=LAT)
  Y      <- numeric(NT)
  for (i in 1:NT){
    Z       <- sst[,,i]
    obj     <- list(x=lon, y=lat, z=Z)
    SURFACE <- interp.surface.grid(obj, latlon)
    Y[i]    <- SURFACE$z
  }
  return(Y)
}

#The following function gives the analytic seasonal temperature
days   <- 1:NT
SSTFun <- function(avgT, rangeT, days, offset = 60){
  y    <- cos(2. * pi * (days-offset) / NT)
  temp <- avgT - rangeT/2*y
  return(temp)
}

#Calculate mean sst
sst_mean <- apply(sst, c(1,2), mean)

#Range of seasonal temperature
sst_min   <- apply(sst, c(1,2), min)
sst_max   <- apply(sst, c(1,2), max)
sst_range <- sst_max - sst_min

#Model mean temperature and range
SST_avg_mod   <- c(6, 24)
SST_range_mod <- c(6, 12) #Two times Am

SST_mod        <- expand.grid(SST_avg_mod, SST_range_mod)
names(SST_mod) <- c('SST_mean', 'SST_range')
rm(SST_avg_mod, SST_range_mod)

#Find the grid closest to the target mean temperature and temperature ranges
set.seed(1000)
SST_mod$Lon <- NA
SST_mod$Lat <- NA

for (k in 1:nrow(SST_mod)){
  dT <- 0.15
  
  ss <- which(abs(sst_mean -SST_mod[k,1]) <= dT &
              abs(sst_range-SST_mod[k,2]) <= dT,  arr.ind=T)
  n  <- sample(1:nrow(ss), 1)
  SST_mod[k, 'Lon'] <- lon[ss[n,1]]
  SST_mod[k, 'Lat'] <- lat[ss[n,2]]
}

#Extract the daily temperature
SST1 <- matrix(NA, nr = nrow(SST_mod), nc = NT)
for (i in 1:nrow(SST_mod)){
  SST1[i,] <- stn_sst(SST_mod[i,'Lon'], SST_mod[i,'Lat'])
}

#Plot mean
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color
pdf('Fig5range.pdf', width=6, height=10)
nf <- layout(matrix(c(1,1,2,2,3:6), 4, 2, byrow = T)) #Change plotting layout
op <- par(font.lab = 1, 
          family ="serif",
          mar    = c(3.5,4,2,2),
          mgp    = c(2.3,1,0),
          oma    = c(1,1.5,0,0),
          cex.axis=1.3,
          cex.lab =1.3) 

lon1 <- seq(0, 360, by = 50)
lon0 <- lon1
lon1[lon1 > 180] <- lon1[lon1 > 180] - 360
image2D( sst_mean, lon, lat, xaxt = 'n', 
        xlab = 'Longitude (ºE)', ylab = 'Latitude (ºN)',
        clab = '')
axis(1, at = lon0, labels = lon1)
mtext('a) Annual mean temperature (ºC)', adj=0)

image2D( sst_range, lon, lat, xaxt = 'n', 
        xlab = 'Longitude (ºE)', ylab = 'Latitude (ºN)',
        clab = '')
axis(1, at = lon0, labels = lon1)
mtext('b) Temperature range (ºC)', adj=0)

#Add the four stations
for (k in 1:nrow(SST_mod)){
  text(SST_mod$Lon[k], SST_mod$Lat[k], 
       LETTERS[k], 
       col=2, 
       cex=1.2)
}

#Add the seasonal variations of temperature at these four stations
SST_mod$offset <- c(60, 262.5, 59.2, 42)
for (i in 1:nrow(SST_mod)){
  range1 <- max(SST1[i,]) - min(SST1[i,])
  mean1  <- mean(SST1[i,])
  SSTsim <- SSTFun(mean1, range1, days, offset= SST_mod[i, 'offset'])
  sstrange <- range(c(SST1[i,],SSTsim))
  plot(days, SST1[i,], ylim=sstrange, type = 'l',
       xlab = '', ylab = '')
  points(days, SSTsim, type = 'l', lwd = .9, lty = 2)
  
  # #Try spline
  # df <- data.frame(t=as.numeric(days), sst=SST1[i,])
  # cl <- gam(sst ~ s(t, bs = 'cc', k = 4), data = df)
  # vt <- predict(cl, newdata = data.frame(t=days))
  # points(days, vt, type = 'l', lwd = .9, lty = 1, col = 2)
  
  #txt <- paste(LATLONS[i,1], 'ºE', LATLONS[i,2],'ºN')
  mtext(paste0(letters[i+2],') St. ', LETTERS[i]), adj=0)
  if (i==1) legend('topleft', c('Observed','Simulated'),lty=1:2, cex=1.2)
}
mtext('Days',adj=0.5, side=1, outer=T)
mtext('Temperature (ºC)', adj=.25, side=2, outer=T)
dev.off()