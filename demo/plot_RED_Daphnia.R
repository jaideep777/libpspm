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

dat = read.delim("/home/jaideep/codes/libpspm/demo/RED_model/fmu_Redmodel.txt", header=F)
dat = dat[,-ncol(dat)]
N = 100
x = exp(seq(log(1), log(1e6), length.out=N))

matplot(x=x, y = t(dat[seq(1,nrow(dat),by=1),-c(1,2)]), col=rainbow(nrow(dat)/1, start=0, end=0.9, alpha=10/nrow(dat)), type="l", lty=1, log="xy", ylim=c(1e-20, 1e4))
xeq = exp(seq(log(1), log(1e6), length.out=30))
points(x=xeq, y= Ueq(xeq))
# plot(dat$V2~dat$V1, type="l")


### DAPHNIA ####

dat = read.delim("/home/jaideep/codes/libpspm/demo/Daphina_model/ebt_Daphnia.txt", header=F)
dat = dat[,-ncol(dat)]
a = 0.75
mu = 0.1
r = 0.5
K = 3
xstar = (mu*(1+mu)*(2+mu)/2/a)^(1/3)
sstar = xstar/(1-xstar)
N = 30
x=seq(0,1,length.out = N)
# plot(x=x, y=exp(-8*x^3), type="l")

matplot(x=seq(0,1,length.out = N), y = t(dat[seq(1,nrow(dat),by=1),-c(1,2)]), col=rainbow(nrow(dat)/1, start=0, end=0.9, alpha=10/nrow(dat)), type="l", lty=1, log="y", ylim=c(0.1, 1e4))
abline(v=xstar, col="black")
xeq = seq(0,1,length.out = 30)
points(x=xeq, y= a*r*sstar*(1-sstar/K)*(xstar-xeq)^(mu-1)/xstar^mu)
# plot(dat$V2~dat$V1, type="l")
