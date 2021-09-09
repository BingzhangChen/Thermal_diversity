#Plot distributions of relative abundances
NPHY <- 41
Topt <- -5 + 1:NPHY

pdf('Fig6RelBio.pdf',width=7,height=7)

par(font.lab  = 1,
    family    = "serif",
    cex.axis  = 1.2,
    cex.lab   = 1.2,
    mgp       = c(2.2,1,0),
    mfcol     = c(2,2),
    oma       = c(2,2,0,0),
    mar       = c(2,4,1,.2))

L    <- 0
envr <- c('Cold','Warm')
varb <- c('low', 'high')
for (n in 1:2){
  for (i in 1:2){
    L     <- L + 1
    for (j in 1:2){
      cff   <- get_data(i,j,n)[ff, ]
      iname <- which(names(cff) == 'PHY1')
      cff   <- cff[, iname:(iname+NPHY-1)]
  
      #Relative abundance
      cffs  <- get_data(i,j,n)[ff,'PHY' ] 
      cffr  <- cff/cffs 
  
      #Obtain average
      cffm  <- apply(cffr, 2,mean)
  
      if (j == 1){
         plot(Topt, cffm, type = 'h', lwd=2, 
              ylim=c(0,0.18),
              xlab= '',
              ylab= '')
         mtxt <- paste0(letters[L], ') ', envr[i], ', ',varb[n], ' variability')
         mtext(mtxt, adj=0, cex=1.1)
      }else{
         points(Topt-0.5, cffm, type="h", lwd=2, col=2)
      }
      legendtxt <- paste0(c('Low','High'), ' diversity')
      if (n==1 && i==1) legend("topright", legendtxt, col=1:2, lwd=3)
    }
  }
}
XLAB <- expression(paste(T[opt] * " (ºC)"))
YLAB <- "Relative biomass"
mtext(XLAB, side = 1,outer=T, line=.2,adj=.5)
mtext(YLAB, side = 2,outer=T, line=0,adj=.5)
dev.off()
