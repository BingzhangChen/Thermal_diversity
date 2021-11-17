library(Rcpp)
sourceCpp('utils.cpp')
#alphaG=1.9,    high diversity, species can coexist
#alphaG=1.1,  low diversity, competitive exclusion

#Read data:
avgtemp <- c(6, 24)
alphaG  <- c(1.0, 1.9)
NaG     <- length(alphaG)
Ntemp   <- length(avgtemp)
NAm     <- 2
Am      <- c(3,6)
nyear   <- 5
d_per_y <- 365
newx    <- 1:d_per_y
get_data <- function(i,j,k, u){
  #i: average temp
  #j: Diversity (alphaG)
  #k: amplitude
  #u: Activation energy of Kn
  sst  <- avgtemp[i]
  file <- paste0('Mod',sst,'_',j,'_',k,'_', u, '.out')
  dat  <- read.table(file, header=T)
  dat$TN <- dat$NO3+dat$ZOO+dat$PHY
  return(dat)
}

#system("scp qsb18204@buchanan.mathstat.strath.ac.uk:~/OneDrive/Kremer/ideal/Mod24*.out ./")
#system("scp qsb18204@buchanan.mathstat.strath.ac.uk:~/OneDrive/Kremer/ideal/Mod6*.out ./")
#Obtain final year
Nout <- nrow(get_data(1,1,1,1))
NY   <- (Nout-1)/nyear

#Plot only final year 
ff     <- (Nout-NY+1):Nout
months <- c("J","M","M","J","S","N")
dates  <- seq(15,365, by = 60)
U      <- 2 #Change U = 2 to plot results from temperature-dependent Kn
fname  <- paste0('FigS1NPZ',U,'.pdf')
pdf(fname,width=14,height=8)
par(font.lab  = 1,
    family    = "serif",
    cex.axis  = 1.7,
    cex.lab   = 1.7,
    mgp       = c(2.2,1,0),
    mfcol     = c(3,4),
    oma       = c(2,2,0,0),
    mar       = c(2,4,2,.2))
for (n in 1:2){  #Am
   for (i in 1:2){  #Mean temp
      #add environmental temperature
      temper <- get_data(i,1,n, U)$Temp[ff]
   
      YLIM <- range(c(get_data(i,1,n,U)$NO3[ff], get_data(i,2,n,U)$NO3[ff]))
      let  <- letters[(1+(i-1)*3 + (n-1)*6):(3+(i-1)*3 + (n-1)*6)]
      if (i == 1 && n == 1){
         YLAB <- c('DIN (µM)', 
                   expression(paste(P[t] *' (µM)')), 'ZOO (µM)')   
         Lab  <- 'Cold, low variability'
      }else{
         YLAB <- ''
         if (i == 2 && n == 1){
            Lab  <- 'Warm, low variability'
         }else if (i == 1 && n == 2){
            Lab  <- 'Cold, high variability'
         }else if (i == 2 && n == 2){
            Lab  <- 'Warm, high variability'
         }else{
            stop ("Index is wrong!")
         }
      }
   
      plot(newx, get_data(i,1,n,U)$NO3[ff], xaxt='n',
          xlab='',ylab=YLAB[1],
          ylim=YLIM,type = 'l')
      points(newx, get_data(i,2,n,U)$NO3[ff], type='l', col=2)
      axis(side=1, at= dates, label=months)
      mtext(let[1], adj=0, cex=1.5)
      mtext(Lab, cex=1.5)
      legendtxt = paste0(c('Low','High'), ' diversity')
      if (i==1&&n==1) legend("topleft", legendtxt, col=1:2, lty=1, cex=1.5)
   
      #Phytoplankton biomass plot
      YLIM <- range(c(get_data(i,1,n,U)$PHY[ff], get_data(i,2,n,U)$PHY[ff]))
      plot(newx, get_data(i,1,n,U)$PHY[ff], xlab='', xaxt='n',
          ylab=YLAB[2],
          ylim=YLIM,type = 'l')
      mtext(let[2],adj=0,cex=1.4)
      points(newx, get_data(i,2,n,U)$PHY[ff], type='l', col=2)
      axis(side=1, at=dates, label=months)
   
      #Zooplankton biomass plot
      YLIM <- range(c(get_data(i,1,n,U)$ZOO[ff], 
                      get_data(i,2,n,U)$ZOO[ff]))
      plot(newx, get_data(i,1,n,U)$ZOO[ff], 
          xlab='', 
          xaxt='n',
          ylab=YLAB[3],
          ylim=YLIM,type = 'l')
      mtext(let[3],adj=0,cex=1.4)
      points(newx, get_data(i,2,n,U)$ZOO[ff], type='l', col=2)
      axis(side=1, at=dates, label=months)
   }
}
mtext(side=1,'Month', outer=T, line=.8,cex=1.4)
dev.off()

#Fig. 6 seasonal variations of mean Topt, sigmaZ, mumax, and Ew (No temperature dependence on Kn)
U    <- 2  #Change U = 2 to plot results from temperature-dependent Kn
NROW <- 4  #NROW figures

fname <- paste0('Fig6avgTopt',U,'.pdf')
pdf(fname,width=14,height=3*NROW)
par(font.lab  = 1,
    family    = "serif",
    cex.axis  = 1.5,
    cex.lab   = 1.5,
    mgp       = c(2.2,1,0),
    mfcol     = c(NROW,4),
    oma       = c(2,2,0,0),
    mar       = c(2,4,2,.2))
for (n in 1:2){
  for (i in 1:2){
     #add environmental temperature
     temper <- get_data(i,1,n,U)$Temp[ff]
  
     #Mean Topt
     YLIM <- range(c(temper, TKrev(get_data(i,1,n, U)$avgZ[ff]), 
                             TKrev(get_data(i,2,n, U)$avgZ[ff])))
  
     let  <- letters[(NROW*(i-1)+1+NROW*2*(n-1)):(i*NROW+NROW*2*(n-1))]
     if (i == 1 && n == 1){
        YLAB <- c(expression(paste(bar(T[opt]) * " (ºC)")),
                  expression(paste(sigma[theta])),
                  expression(paste(µ[w] *' ('*d^-1* ')')),
                  expression(paste(E[w] *' (eV)'))
  		           )
        Lab  <- 'Cold, low variability'
     }else{
        YLAB <- ''
        if (i == 1 && n == 2){
           Lab  <- 'Cold, high variability'
        }else if (i == 2 && n == 1){
           Lab  <- 'Warm, low variability'
        }else if (i == 2 && n == 2){
           Lab  <- 'Warm, high variability'
        }else{
           stop ("Index is wrong!")
        }
     }
  
     k <- 1
     plot(newx, TKrev(get_data(i,1,n, U)$avgZ[ff]), 
         xaxt='n',
         xaxs='i',
         xaxs='i',
         xlab='', 
         ylab=YLAB[k],
         ylim=YLIM,type = 'l')
     mtext(let[k],adj=0,cex=1.5)
     points(newx, TKrev(get_data(i,2,n, U)$avgZ[ff]), type='l', col=2)
     axis(side=1, at=dates, label=months)

     lines(newx, temper, lty=2)
     mtext(Lab, cex=1.5)
     legendtxt <- paste0(c('Low','High'), ' diversity')
     if (i==1&&n==1) legend("topleft", legendtxt, col=1:2, lty=1, cex=1.5)
  
    #SigmaZ
     k <- k + 1
     YLIM <- c(0, 3)
     plot(newx, get_data(i,1,n, U)$sigmaZ[ff], xlab='', xaxt='n',
         ylab=YLAB[k],
         ylim=YLIM,type = 'l')
  
     mtext(let[k],adj=0,cex=1.5)
     points(newx, get_data(i,2,n,U)$sigmaZ[ff], type='l', col=2)
     axis(side=1, at=dates, label=months)
  
     #mumax
     k    <- k + 1
     YLIM <- c(.1, 1.5)
     plot(newx, get_data(i,1,n,U)$mumax[ff], xlab='', xaxt='n',
         ylab=YLAB[k],
         ylim=YLIM,type = 'l')
  
     mtext(let[k],adj=0,cex=1.5)
     points(newx, get_data(i,2,n,U)$mumax[ff], type='l', col=2)
     axis(side=1, at=dates, label=months)
  
     #Ew
     k <- k + 1
     YLIM <- c(-0.6, 1.1)
     plot(newx, get_data(i,1,n,U)$Ew[ff], xlab='', xaxt='n',
         ylab=YLAB[k],
         ylim=YLIM,
         type='l')
  
     mtext(let[k],adj=0,cex=1.5)
     points(newx, get_data(i,2,n, U)$Ew[ff], type='l', col=2)
     abline(h=0.32, lty=2)
     abline(h=0.65, lty=2)
     if (i == 1){
        text(200, 0.32, "0.32 eV", cex = 1.2)
        text(200, 0.64, "0.65 eV", cex = 1.2)
     }
     axis(side=1, at=dates, label=months)
  }
}
mtext(side=1,'Month', outer=T, line=.8,cex=1.5)
dev.off()