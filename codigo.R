

# ------------------------------------------
# Algoritmo de optimizacion MV 
# Modelo MV-ESG 
# Autor: Carlos Zapata
# ------------------------------------------

rm(list=ls())
library(quantmod)
library(xts)
library(zoo)

# Parametros 
mu <- colMeans(retornos)
cov <- cov(retornos)
var <- diag(cov)
sigma <- sqrt(var)
corr <- cor(retornos)
n <- length(mu)

# ESG scores - Source: Refinitiv
iscore <- c(32,60,50,51)

# ------------------------------------------------------
# ------------------------------------------------------

# Optimizacion de los portafolios
# Modelo de markowitz

unos <- rep(1,n)
xMV <- t(mu)%*%solve(cov)%*%mu
yMV <- t(mu)%*%solve(cov)%*%unos
zMV <- t(unos)%*%solve(cov)%*%unos
dMV <- xMV*zMV-yMV^2
gMV <- (solve(cov,unos)%*%xMV-solve(cov,mu)%*%yMV)%*%(1/dMV)
hMV <- (solve(cov,mu)%*%zMV-solve(cov,unos)%*%yMV)%*%(1/dMV)

# Frontera eficiente de Markowitz
nport <- 100
rp <- seq(min(mu),max(mu),length=nport)
wpoMV <- matrix(0,ncol=n,nrow=nport)
sigmapoMV <- matrix(0,nrow=nport)
rpoMV <- matrix(0,nrow=nport)

for(i in 1:nport){
    w <- gMV+hMV*rp[i]
    wpoMV[i,] <- t(w)
    sigmapoMV[i,] <- sqrt(t(w)%*%cov%*%w)
    rpoMV[i,] <- t(w)%*%mu
}
# PMVG
wpmvg <- solve(cov,unos)%*%(1/zMV)
rownames(wpmvg) <- activos
rpmvg <- t(wpmvg)%*%mu
sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)

windows()
plot(sigma,mu,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.1),ylim=c(0,max(mu)),
     xlab="Riesgo",ylab="Retorno esperado")
lines(sigmapoMV,rpoMV,type="l")
points(sigmapmvg,rpmvg,pch=19)
text(sigma,mu,labels = activos,pos=4,cex=0.6)


# ------------------------------------------------------
# ------------------------------------------------------

# Modelo MV-ESG

unos <- rep(1,n)
x <- t(mu)%*%solve(cov)%*%mu
y <- t(mu)%*%solve(cov)%*%unos
z <- t(unos)%*%solve(cov)%*%unos

# Extension ESG
p <- t(mu)%*%solve(cov)%*%iscore
q <- t(unos)%*%solve(cov)%*%iscore
r <- t(iscore)%*%solve(cov)%*%iscore

# Calculo d ajustado
d <- (x*z*r+y*p*q+y*p*q-z*p^2-x*q^2-r*y^2)

# Vectores
g <- (solve(cov,mu)%*%(p*q)-solve(cov,mu)%*%(y*r)+solve(cov,unos)%*%(x*r)-solve(cov,unos)%*%(p^2)+
          solve(cov,iscore)%*%(p*y)-solve(cov,iscore)%*%(x*q))%*%(1/d)
h <- (solve(cov,mu)%*%(r*z)-solve(cov,mu)%*%(q^2)+solve(cov,unos)%*%(p*q)-solve(cov,unos)%*%(y*r)+
          solve(cov,iscore)%*%(y*q)-solve(cov,iscore)%*%(p*z))%*%(1/d)
f <- (solve(cov,mu)%*%(y*q)-solve(cov,mu)%*%(z*p)+solve(cov,unos)%*%(y*p)-solve(cov,unos)%*%(q*x)+
          solve(cov,iscore)%*%(z*x)-solve(cov,iscore)%*%(y^2))%*%(1/d)

# Ejemplo verificion
# rp <- 0.015
ip <- 0.4
rp <- 0.01
w <- g+h*rp+f*ip

# Construccion Frontera eficiente - (Superficie)
nport <- 100
rp <- seq(min(mu),max(mu),length=nport)
ip <- seq(min(iscore),max(iscore),length=nport)

wpo <- matrix(0,ncol=n,nrow=nport)
sigmapo <- matrix(0,ncol=nport,nrow=nport)
rpo <- matrix(0,ncol=nport,nrow=nport)
ipo <- matrix(0,ncol=nport,nrow=nport) 

for(i in 1:nport){
    for(j in 1:nport){
        w <- g+h*rp[i]+f*ip[j]
        wpo[i,] <- t(w)
        sigmapo[i,j] <- sqrt(t(w)%*%cov%*%w)
        rpo[i,j] <- t(w)%*%mu
        ipo[i,j] <- t(w)%*%iscore
    }
}

# PMVG
wpmvg <- solve(cov,unos)%*%(1/zMV)
rownames(wpmvg) <- activos
rpmvg <- t(wpmvg)%*%mu
sigmapmvg <- sqrt(t(wpmvg)%*%cov%*%wpmvg)
ipmvg <- t(wpmvg)%*%iscore

# Graficos
windows()
plot(sigma,mu,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.2),ylim=c(0,max(mu)),
     xlab="Riesgo",ylab="Retorno esperado",pch=19,bty="L")
points(sigmapo,rpo,col="gray")
lines(sigmapoMV,rpoMV,lwd=3)
text(sigma,mu,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,rpmvg,pch=19)
text(sigmapmvg,rpmvg,labels = "MVG",cex=0.6,pos=2)
legend("topleft", legend = c("Markowitz", "MV-ESG"),
       lty = c(1,1), lwd=c(2,1),col=c("black","gray"))


windows()
plot(sigma,iscore,main="Plano Riesgo-Retorno",xlim=c(0.04,max(sigma)*1.2),ylim=c(30,max(iscore)),
     xlab="Riesgo",ylab="ESG score",pch=19)
lines(sigmapo,ipo,type="l",col="gray")
text(sigma,iscore,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,ipmvg,pch=19)
text(sigmapmvg,ipmvg,labels = "MVG",cex=0.7,pos=2)


# ------------------------------------------------------------
## Forma convexa
# ------------------------------------------------------------

library(matlab)
library(MASS)
library(CVXR)

# PMVG - Markowitz - ESG
portMarkowitz <- function(mu, cov){
    w <- Variable(n)
    prob <- Problem(Minimize( quad_form(w, cov) ),
                    constraints = list(ones(1,n)%*%w==1, t(w)%*%mu==rpCVX[i], t(w)%*%iscore==ipCVX[j] )) #w >= 0
    result <- solve(prob)
    return(as.vector(result$getValue(w)))
}

# Verifica PMVG - ok
#rpCVX <- rpmvg
#ipCVX <- ipmvg
#w_Markowitz <- round(portMarkowitz(mu, cov),4)
#names(w_Markowitz) <- colnames(cov)
#w_Markowitz <- cbind(w_Markowitz)

nport <- 100
rpCVX <- seq(min(mu),max(mu),length=nport)
ipCVX <- seq(min(iscore),max(iscore),length=nport)

wpoCVX <- matrix(0,ncol=n,nrow=nport)
sigmapoCVX <- matrix(0,ncol=nport,nrow=nport)
rpoCVX <- matrix(0,ncol=nport,nrow=nport)
ipoCVX <- matrix(0,ncol=nport,nrow=nport) 

for(i in 1:nport){
    for(j in 1:nport){
        w <- portMarkowitz(mu, cov)
        wpoCVX[i,] <- w
        sigmapoCVX[i,j] <- sqrt(t(w)%*%cov%*%w)
        rpoCVX[i,j] <- t(w)%*%mu
        ipoCVX[i,j] <- t(w)%*%iscore
    }
}

windows()
plot(sigma,iscore,main="Plano Riesgo-Retorno",xlim=c(0,max(sigma)*1.2),ylim=c(30,max(iscore)),
     xlab="Riesgo",ylab="ESG score",pch=19)
lines(sigmapoCVX,ipoCVX,type="l",col="gray")
text(sigma,iscore,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,ipmvg,pch=19)
text(sigmapmvg,ipmvg,labels = "MVG",cex=0.7,pos=2)

windows()
plot(sigma,iscore,main="Plano Riesgo-Retorno",xlim=c(0.03,max(sigma)*1.2),ylim=c(50,max(iscore)),
     xlab="Riesgo",ylab="ESG score",pch=19)
lines(sigmapoCVX,ipoCVX,type="l",col="gray")
text(sigma,iscore,labels = activos,cex=0.6,pos=4)
points(sigmapmvg,ipmvg,pch=19)
text(sigmapmvg,ipmvg,labels = "MVG",cex=0.7,pos=2)


# ----------------------------------------------------------
# ----------------------------------------------------------
