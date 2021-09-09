#File of SST temperature
#Downloaded from https://www.esrl.noaa.gov/psd/cgi-bin/db_search/DBSearch.pl?Dataset=NOAA+High-resolution+Blended+Analysis&Variable=Sea+Surface+Temperature&group=0&submit=Search on 22 Dec 2018

source('~/Working/FlexEFT1D/Rscripts/Interpolate_WOA.R')
library(plot3D)
library(MASS)
tfile <- "sst.day.mean.ltm.1971-2000.nc"

lon  <- ncread(tfile, 'lon')
lat  <- ncread(tfile, 'lat')
time <- ncread(tfile, 'time')
sst  <- ncread(tfile, 'sst')

#Calculate mean sst
sst_mean <- apply(sst, c(1,2), mean)

#Range
sst_min <- apply(sst, c(1,2), min)
sst_max <- apply(sst, c(1,2), max)
sst_range <- sst_max - sst_min

#Plot mean
jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))  #A good color
pdf('FigS2range.pdf', width=6, height=10)
op <- par(font.lab = 1, 
          family ="serif",
          mar    = c(3.5,4,2,2),
          mgp    = c(2.3,1,0),
          cex.axis=1.3,
          cex.lab =1.3,
          mfrow  = c(3,1)   ) 

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

f2   <- data.frame(x=as.vector(sst_mean), y=as.vector(sst_range))
f2   <- na.omit(f2)
f2   <- f2[f2$y > 0.99,]
f2   <- kde2d(f2$x, f2$y, n = 50)
ZLIM <- c(0, 0.03)
f2$z[f2$z > ZLIM[2]] <- ZLIM[2]
image2D(f2, col = jet.colors(20), ylim=c(0,25), zlim=ZLIM,
     xlab='Mean temperature (ºC)',ylab='Temperature range (ºC)')
mtext('c)',adj=0)
dev.off()

#Plot selected stations
library(fields)
NT      <- dim(sst)[3]
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

LATLONS <- matrix(NA, nr = 2, nc = 2)
LATLONS[,1] <- 150  #Longitude
LATLONS[,2] <- c(20, 40)  #Latitude

SST1 <- matrix(NA, nr = 2, nc = NT)

for (i in 1:2){
    SST1[i,] <- stn_sst(LATLONS[i,1], LATLONS[i,2])
}

days   <- 1:NT

SSTFun <- function(avgT, rangeT, offset = 60){
    
  y   = cos(2. * pi * (days-offset) / NT)
  temp= avgT - rangeT/2*y
  return(temp)
}
#Plot an example to show the fits of seasonal temperature changes
pdf('tempcos.pdf', width = 5, height = 7)
op <- par(font.lab = 1, 
          family ="serif",
          mar    = c(1.5,2,2,.2),
          oma    = c(4,4,0,0),
          mgp    = c(2.3,1,0),
          cex.axis=1.1,
          cex.lab =1.1,
          mfrow  = c(2,1)   ) 

for (i in 1:2){
  range1 <- max(SST1[i,]) - min(SST1[i,])
  mean1  <- mean(SST1[i,])
  SSTsim <- SSTFun(mean1, range1)
  sstrange <- range(c(SST1[i,],SSTsim))
  plot(days, SST1[i,],ylim=sstrange, type = 'l',
       xlab = '', ylab = '', pch = 16, cex=.7)
  points(days, SSTsim, type = 'l', lwd = .9, lty = 2)
  txt <- paste(LATLONS[i,1], 'ºE', LATLONS[i,2],'ºN')
  mtext(paste(letters[i],')',txt),adj=0)
  if (i==1) legend('topleft', c('Observed', 'Simulated'),lty=1:2 )
}
mtext('Days',adj=0.5, side=1, outer=T)
mtext('Temperature (ºC)', adj=.5, side=2, outer=T, line=1)
dev.off()
