################################
# Code for analyzing intra- vs. interspecific growth rate ~ temperature relationships
# in marine phytoplankton
# Modified by B. Chen on 25/8/2021

####################################
T.K  <- function(x, Tref = 15,kb=8.62E-5) 1/((Tref+273)*kb) - 1/((x + 273)*kb)
load('PHYnls.Rdata')

PHYnls$ivTopt <- T.K(PHYnls$Topt)
Lm_um <- lm(Lnum ~ ivTopt, PHYnls)
Lm_Ea <- lm(LnEa ~ ivTopt, PHYnls)
Lm_Eh <- lm(LnEh ~ ivTopt, PHYnls)

pdf('Fig1PHYTrait_cor.pdf', height = 3, width = 8)
op <- par(font.lab = 1,
           family  = "serif",
           cex.lab = 1.2,
           cex.axis= 1.2,
           mgp     = c(2.2,1,0),
           oma     = c(2,2,0,0),
           mar     = c(2.5,4.3,2.3,.5),
           mfrow   = c(1,3))
x1   <- seq(0,30,10)
y1   <- T.K(x1)
x1[length(x1)] <- paste0(x1[length(x1)],'ºC')

plot(PHYnls$ivTopt, PHYnls$Lnum, 
     cex  = 0.5,
     xlab = '',
     ylab = expression(paste('Ln '*italic(µ)[m])))
mtext('a', adj=0, cex=1.4)
axis(3, at= y1, label = x1)
df  <- coef(summary(Lm_um))[, "Estimate"]
abline(df[1], df[2], lwd=2, lty = 1)

plot(PHYnls$ivTopt, PHYnls$LnEa, 
     cex  = 0.5,
     xlab = '',
     ylab = bquote('Ln' ~italic(Ea) ))
axis(3, at= y1, label = x1)
mtext('b', adj=0, cex=1.4)
df  <- coef(summary(Lm_Ea))[ , "Estimate"]
abline(df[1], df[2], lwd=2)

plot(PHYnls$ivTopt, PHYnls$LnEh,
     cex  = 0.5,
     xlab = '',
     ylab = bquote('Ln '~italic(E[h]) ))
mtext('c', adj=0, cex=1.4)
axis(3, at= y1, label = x1)

df  <- coef(summary(Lm_Eh))[, "Estimate"]
abline(df[1], df[2], lwd=2)

mtext(expression(paste('Optimal temperature (',italic(theta),')')),
      side=1, outer=T, line=.2, cex=1.3)
dev.off()

