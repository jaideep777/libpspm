library(tidyverse)

setwd("~/codes/libpspm/demo/")

png("RED_Daphnia.png", width = 1200*3, height=680*3, res=300)

par(mfrow = c(2,6), mar=c(4,4,4,1), oma=c(1,1,1,1), cex.lab=1.2, cex.axis=1.2)

### RED ####
# Analytical calculation of equilibrium distribution, if available
Ueq = function(x){
  a0<- 0.396
  phiA<- 0.749
  g0<-0.0838
  phiG<-0.7134
  m0<-1
  mort<-0.035
  mu0<-mort*m0/g0
  alpha<-0.10
  temp<- mu0/(1-phiG)
  coverage<-1-(1-alpha)/alpha*mu0/((mu0/(1-phiG))^(phiG/(phiG-1))*exp(mu0/(1-phiG))*as.numeric(pracma::gammainc(mu0/(1-phiG),phiG/(1-phiG)+1)[2]))
  Neq<-coverage/a0/(temp^(phiA/(phiG-1))*exp(temp)*as.numeric(pracma::gammainc(temp,phiA/(1-phiG)+1))[2])
  n0<-Neq*mort/g0
  return(n0*(x/m0)^(-phiG)*exp(mu0/(1-phiG)*(1-(x/m0)^(1-phiG)))*10000)
}

xeq = exp(seq(log(1), log(1e6), length.out=1000))
ueq = Ueq(xeq)

N = integrate(Ueq, 1, 1e6, abs.tol = 1e-6, rel.tol = 1e-6)
B = integrate(function(x){x*Ueq(x)}, 1, 1e6, abs.tol = 1e-6, rel.tol = 1e-6)

plot1 = function(file, N, title){
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  x = exp(seq(log(1), log(1e6), length.out=N))
  
  cols = rainbow(nrow(dat)/1, start=0, end=0.9, alpha=10/nrow(dat))
  
  matplot(x=x, y = t(dat[seq(1,nrow(dat),by=1),-c(1,2)]), col=cols, 
          type="l", lty=1, log="xy", ylim=c(1e-20, 1e4), xlab="Size", ylab="Density", main=title)
  
  lines(x=xeq, y=ueq, col=scales::alpha("orange", 0.7), lwd=4)
  points(x=x, y = t(dat[nrow(dat),-c(1,2)]), col="black", 
         type="l", lty=1, lwd=1, cex = 0.7)
}

plot1("RED_model/fmu_Redmodel.txt", 150, "FMU")
plot1("RED_model/ifmu_Redmodel.txt", 150, "IFMU")
plot1("RED_model/ifmu2_Redmodel.txt", 150, "IFMU(O2)")
plot1("RED_model/ebt_Redmodel.txt", 150, "EBT")
plot1("RED_model/cm_Redmodel.txt", 150, "CM")
plot1("RED_model/abm_Redmodel.txt", 150, "ABM")



### DAPHNIA ####

a = 0.75
mu = 0.1
r = 0.5
K = 3

xstar = (mu*(1+mu)*(2+mu)/2/a)^(1/3)
sstar = xstar/(1-xstar)
xeq = seq(0,1,length.out = 1000)
u_equil = function(x){
  a*r*sstar*(1-sstar/K)*(xstar-x)^(mu-1)/xstar^mu
} 
ueq = u_equil(xeq)

N = integrate(u_equil, 0, xstar*0.99999999999, abs.tol = 1e-6, rel.tol = 1e-6)
B = integrate(function(x){x*u_equil(x)}, 0, xstar*0.99999999999, abs.tol = 1e-6, rel.tol = 1e-6)

plot2 = function(file, N, title){
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  x=seq(0,1,length.out = N)
  # plot(x=x, y=exp(-8*x^3), type="l")
  
  cols = rainbow(nrow(dat)/1, start=0, end=0.9, alpha=10/nrow(dat))
  # cols = scales::viridis_pal(direction = -1, alpha = 10/nrow(dat))(nrow(dat))
  matplot(x=x, y = t(dat[seq(1,nrow(dat),by=1),-c(1,2,3)]), col=cols, 
          type="l", lty=1, log="y", ylim=c(0.1, 1e4), xlab="Size", ylab="Density", main=title)
  abline(v=xstar, col="grey")
  
  lines(x=xeq, y=ueq, col=scales::alpha("orange", 0.7), lwd=4)
  points(x=x, y = t(dat[nrow(dat),-c(1,2,3)]), col="black", 
         type="l", lty=1, lwd=1, cex = 0.7)
  # plot(dat$V2~dat$V1, type="l")
}

plot2("Daphnia_model/fmu_Daphnia.txt", 300, "FMU")
plot2("Daphnia_model/ifmu_Daphnia.txt", 300, "IFMU")
plot2("Daphnia_model/ifmu2_Daphnia.txt", 300, "IFMU(O2)")
plot2("Daphnia_model/ebt_Daphnia.txt", 300, "EBT")
plot2("Daphnia_model/cm_Daphnia.txt", 300, "CM")
plot2("Daphnia_model/abm_Daphnia.txt", 300, "ABM")

dev.off()
