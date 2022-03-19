Ueq = function(x){
  (1-x)^2/(1+x)^4
}

plot1 = function(file, title){
  x = seq(0,1,length.out=26)
  dat = read.delim(file, header=F)
  dat = dat[,-ncol(dat)]
  Y = dat
  matplot(x=x, y=t(Y[,-c(1,2)]), type = "l", lty=1, col=rainbow(nrow(Y)/1, start=0, end=0.9, alpha=10/nrow(Y)), xlab="Size", ylab="Density", main=title, log="")
  points(x=x, y=Ueq(x))
}

par(mfrow = c(2,2), mar=c(4,4,4,1), oma=c(1,1,1,1), cex.lab=1.2, cex.axis=1.2)
plot1("~/codes/libpspm/cm_testmodel_equil.txt", "CM")
plot1("~/codes/libpspm/ebt_testmodel_equil.txt", "EBT")
plot1("~/codes/libpspm/fmu_testmodel_equil.txt", "FMU")
plot1("~/codes/libpspm/ifmu_testmodel_equil.txt", "IFMU")

