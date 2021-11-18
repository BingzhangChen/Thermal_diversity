setwd('~/OneDrive/Kremer/ideal')
temp    <- seq(0, 30, 0.1)
T.K     <- function(x, Tref = 15,kb=8.62E-5) 1/((Tref+273)*kb) - 1/((x + 273)*kb)
temp    <- T.K(temp)
numtemp <- length(temp)


###Growth rate function 
MU0 <- function(x, Z = T.K(20), 
                u0 = 0.59, 
                Ei = 0.22, 
                Ea0= 0.98, 
                beta= -0.2, 
                Eh0 = 4.43,
                phi = 0.21){
   #x is transfored temperature
   
   umax  <- u0  * exp(Ei  *Z)
   Ea    <- Ea0 * exp(beta*Z)
   Eh    <- Eh0 * exp(phi *Z)
   Ed    <- Eh - Ea
   B     <- 1 + Ea/Ed*exp(Eh*(x-Z))
   mu    <- umax*Eh/Ed * exp(Ea*(x-Z))/B
   return(mu)
}

#Check function
# y <- MU0(temp, Z=T.K(-5))
# plot(temp, y, type = 'l')
#Plot three typical TPCs and the fitness landscape under a given environmental temperature
pdf('Fig2example.pdf', height=8, width=4)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           mar     = c(4,3.5,3,1),
           mfrow   = c(2,1))


Topto <- c(5, 15, 25)
Topt_ <- T.K(Topto)
x1    <- seq(0,30,10)
y1    <- T.K(x1)
x1[length(x1)] <- paste0(x1[length(x1)],' ºC')

plot(temp,temp**2, ylim=c(0, 1), 
     xlab = 'Environmental temperature',
     ylab = bquote('Growth rate ('~d^-1~')'),
     type = 'n')
axis(3, at= y1, label = x1)

for (j in 1:length(Topt_)){ 
   mu_ <- MU0(temp, Z=Topt_[j])
   points(temp, mu_, type = 'l', col=j)
}
###

txt <- sapply(1:3, function(i) {
   as.expression(bquote(T[opt]~' = '~.(Topto[i])~'ºC'))
})
legend('topleft', txt, lty=1, col=1:length(Topt_) )

mtext('a', adj=0, cex=1.4)
#Add the envelop line
#curve(mumax, from = -3, to = 3, lty = 2, add = T)

par(mar = c(4,3.5,1,1))
#Plot fitness curve under constant environmental temperatures
Topto <- seq(-10, 50, 1)
Topt_ <- T.K(Topto)
envT0 <- c(5, 15, 25)  #Environmental temperature
envT  <- T.K(envT0)
mu_   <- numeric(length(Topt_))

plot(temp,temp**2, ylim=c(0, 1.5), 
     xlab = bquote(italic(theta)),
     ylab = bquote('Growth rate ('~d^-1~')'),
     type = 'n')

for (j in 1:length(envT)){ 
    mu_  <- MU0(envT[j], Z=Topt_)
    points(Topt_, mu_, type = 'l', col=j)
}
txt = paste('T =', envT0,'ºC')
legend('topright', txt, lty=1, col=1:length(envT) )
mtext('b', adj=0, cex=1.4)
dev.off()

#Below for Fig. 3: community activation energy for different thermal diversity
#Discretize (No approximation)
pdf_norm <- function(x, mean=T.K(20), V) exp(-(x-mean)**2/(2*V))/sqrt(2*pi*V)

#Topt of each species
Topt <- seq(T.K(-1), T.K(50), length.out=100)
Nsp  <- length(Topt)

#Compute the fraction of each species at different variance.
biof <- function(V){
   pfL <- vapply(1:Nsp, function(x)pdf_norm(x=Topt[x], V=V), 0)
   
   #dx of each species
   dx <- diff(Topt)
   dx <- c(dx, dx[length(dx)])
   
   #Compute fraction of each species
   PHYt<- sum(pfL*dx)
   biof<- pfL*dx/PHYt
   return(biof)
}
#Growth rate of each species at each temperature
VARS  <- c(0.1,1,5)
gr    <- array(NA, c(Nsp, length(temp),length(VARS)))
mucom <- array(NA, c(length(temp),length(VARS)))
AE    <- array(NA, c(length(temp)-1,length(VARS)))
for (k in 1:length(VARS)) {
   biomass <- biof(VARS[k])
   for (i in 1:length(temp)) {
      for (j in 1:Nsp) {
         gr[j, i, k] <- MU0(x = temp[i], Z = Topt[j])
      }
      mucom[i, k] <- sum(gr[, i, k] * biomass)
      
      #Calculate apparent E (derivative)
      if (i > 1) AE[i - 1, k] <-
         (log(mucom[i, k]) - log(mucom[i-1, k]))/(temp[i] - temp[i-1])
   }
   
}

pdf('Fig3var_mu.pdf', height=8, width=4)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           mar     = c(4,3.5,3,1),
           mfrow   = c(3,1))

#Plot probability distribution:
x1   <- seq(0,30,10)
y1   <- T.K(x1)
plot(Topt, biof(VARS[1]), type = 'l', 
     xlab = bquote(italic(theta)),
     ylab = 'Biomass density')
x1[length(x1)] <- paste0(x1[length(x1)],' ºC')
axis(3, at= y1, label = x1)

txt <- sapply(1:3, function(i) {
       as.expression(bquote(sigma[theta]^2~' = '~.(VARS[i])))
    })
legend('topright', legend=txt,
       lty=1, col=1:3 )
mtext('a', adj=0, cex=1.4)

for (i in 2:length(VARS)){
     points(Topt, biof(VARS[i]), type = 'l', col=i)
}

#Plot µ ~ temperature
par(mar = c(3.5,3.5,1,1))
plot(temp, log(mucom[,1]), xaxt='n', 
     ylim = c(-3,0.2),
     xlab = 'Environmental temperature (ºC)',
     ylab = expression(paste('Ln '*italic(µ)[w])),
     type = 'l')
axis(1, at=y1, label = x1)
mtext('b', adj=0, cex=1.4)

for (i in 2:length(VARS)) points(temp, log(mucom[,i]), type = 'l', col=i)

#Plot AE ~ temperature
plot(temp[-length(temp)], AE[,1], xaxt='n', 
     ylim = c(-1,1),
     xlab = 'Environmental temperature (ºC)',
     ylab = bquote(E[w] ~ '(eV)'),
     type = 'l')
axis(1, at= y1, label = x1)
abline(h=0.65, lty=2, lwd=.5)
text(1, 0.65, '0.65 eV')
abline(h=0.32, lty=2, lwd=.5)
text(1, 0.32, '0.32 eV')

mtext('c', adj=0, cex=1.4)
for (i in 2:length(VARS)){
   points(temp[-length(temp)], AE[,i], type='l',col=i)
}
dev.off()
