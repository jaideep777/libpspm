
plot_dists = function(folder, title, mtext = F){
  
  up1 = read.delim(paste0(dir, "/", folder, "/species_0_u.txt"), header=F)
  up2 = read.delim(paste0(dir, "/", folder, "/species_1_u.txt"), header=F)
  up3 = read.delim(paste0(dir, "/", folder, "/species_2_u.txt"), header=F)
  
  up1 = as.matrix(up1[,-ncol(up1)])
  up2 = as.matrix(up2[,-ncol(up2)])
  up3 = as.matrix(up3[,-ncol(up3)])
  up1[up1<0]=0
  up2[up2<0]=0
  up3[up3<0]=0
  
  hp1 = read.delim(paste0(dir, "/", folder, "/species_0_X.txt"), header=F)
  hp2 = read.delim(paste0(dir, "/", folder, "/species_1_X.txt"), header=F)
  hp3 = read.delim(paste0(dir, "/", folder, "/species_2_X.txt"), header=F)
  
  hp1 = as.matrix(hp1[,-ncol(hp1)])
  hp2 = as.matrix(hp2[,-ncol(hp2)])
  hp3 = as.matrix(hp3[,-ncol(hp3)])
  
  times = as.numeric(hp1[,1])
  
  image(x=times, y=hp1[1,-1],  z=log(1+100*up1[,-1]), col=scales::viridis_pal()(100), ylab="Height", xlab="Time (years)")
  mtext(title, side=3, line=1)
  if(mtext) mtext("lma = 0.0825", side=2, line=5)
  image(x=times, y=hp2[1,-1],  z=log(1+100*up2[,-1]), col=scales::viridis_pal()(100), ylab="Height", xlab="Time (years)")
  if(mtext) mtext("lma = 0.2625", side=2, line=5)
  image(x=times, y=hp3[1,-1],  z=log(1+100*up3[,-1]), col=scales::viridis_pal()(100), ylab="Height", xlab="Time (years)")
  if(mtext) mtext("lma = 0.4625", side=2, line=5)
}

##### Fixed input seed rain mode

dir = "~/codes/libpspm/demo/Plant_model"
setwd(dir)

# png("../size_dists.png", width = 1000*3, height=750*3, res=300)
par(mar=c(4,4,1,1), oma = c(1,4,2,1), cex.lab=1.2, cex.axis=1.2)
layout(mat = matrix(1:15, nrow=3, byrow=F))

plot_dists("outputs/fmu_120pts", "FMU", T)
plot_dists("outputs/ifmu_100pts", "IFMU (100 pts)")
plot_dists("outputs/ifmu_1000pts", "IFMU (1000 pts)")
plot_dists("outputs/ebt", "EBT")
plot_dists("outputs/abm_n1000_dt0.05_2", "ABM")


setwd(dir)
# dev.off()


setwd(dir)
seeds_fmu = read.delim("outputs/fmu_120pts/seed_rains.txt", header = F)
seeds_ifmu = read.delim("outputs/ifmu_100pts/seed_rains.txt", header = F)
seeds_ifmu1000 = read.delim("outputs/ifmu_1000pts/seed_rains.txt", header = F)
seeds_ebt = read.delim("outputs/ebt/seed_rains.txt", header = F)
seeds_abm = read.delim("outputs/abm_n1000_dt0.05_2/seed_rains.txt", header = F)

plot_seeds = function(y, title, ...){
  matplot(y = y, x=seeds_fmu$V1, type="l", lty=1, lwd=2, col=scales::alpha(c("purple", "green3", "mediumspringgreen", "cyan3", "darkgoldenrod2"), alpha=0.7), ylab="Seed rain", ...)
  mtext(title, line=1)
}

# png("../seed_rains.png", width = 660*3, height=766*3, res=300)
par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2, seeds_ifmu$V2, seeds_ifmu1000$V2, seeds_ebt$V2, seeds_abm$V2), "Species 1", xlab="")
plot_seeds(cbind(seeds_fmu$V3, seeds_ifmu$V3, seeds_ifmu1000$V3, seeds_ebt$V3, seeds_abm$V3), "Species 2", xlab="")
plot_seeds(cbind(seeds_fmu$V4, seeds_ifmu$V4, seeds_ifmu1000$V4, seeds_ebt$V4, seeds_abm$V4), "Species 3", xlab="Time (years)")
# dev.off()




##### Feedback mode


dir = "~/codes/libpspm/demo/Plant_model"
setwd(dir)

png("../size_dists_withFeedback.png", width = 1000*3, height=750*3, res=300)
par(mar=c(4,4,1,1), oma = c(1,4,2,1), cex.lab=1.2, cex.axis=1.2)
layout(mat = matrix(1:15, nrow=3, byrow=F))
  
plot_dists("outputs/fmu_new_f/", "FMU", T)
plot_dists("outputs/ifmu_new_f", "IFMU")
plot_dists("outputs/ifmu2_new_f", "IFMU(O2)")
plot_dists("outputs/ebt_new_f/", "EBT")
plot_dists("outputs/cm_new_f/", "CM")
dev.off()



  
setwd(dir)
seeds_fmu = read.delim("outputs/fmu_new_f/seed_rains.txt", header = F)
seeds_ebt = read.delim("outputs/ebt_new_f/seed_rains.txt", header = F)
seeds_ifmu = read.delim("outputs/ifmu_new_f/seed_rains.txt", header = F)
seeds_ifmu2 = read.delim("outputs/ifmu2_new_f/seed_rains.txt", header = F)
seeds_cm = read.delim("outputs/cm_new_f/seed_rains.txt", header = F)

plot_seeds = function(y, title, ...){
  matplot(y = y, x=seeds_fmu$V1, type="l", lty=1, lwd=2, col=scales::alpha(c("purple", "green3", "mediumspringgreen", "darkgoldenrod2", "magenta"), alpha=0.7), ylab="Seed rain", ...)
  mtext(title, line=1)
}

png("../seed_rains_withFeedback.png", width = 660*3, height=766*3, res=300)
par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2, seeds_ifmu$V2, seeds_ifmu2$V2, seeds_ebt$V2, seeds_cm$V2), "Species 1", xlab="")
plot_seeds(cbind(seeds_fmu$V3, seeds_ifmu$V3, seeds_ifmu2$V3, seeds_ebt$V3, seeds_cm$V3), "Species 2", xlab="")
plot_seeds(cbind(seeds_fmu$V4, seeds_ifmu$V4, seeds_ifmu2$V4, seeds_ebt$V4, seeds_cm$V4), "Species 3", xlab="Time (years)")
dev.off()


par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2, seeds_ebt$V2), "Species 1", xlab="", log="")
plot_seeds(cbind(seeds_fmu$V3, seeds_ebt$V3), "Species 2", xlab="", log="")
plot_seeds(cbind(seeds_fmu$V4, seeds_ebt$V4), "Species 3", xlab="Time (years)", log="")


par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2), "Species 1", xlab="")
plot_seeds(cbind(seeds_fmu$V3), "Species 2", xlab="")
plot_seeds(cbind(seeds_fmu$V4), "Species 3", xlab="Time (years)")

# png("../seed_rains.png", width = 660*3, height=766*3, res=300)
par(mfrow=c(3,1), mar = c(4,4,1,1), oma = c(1,1,4,1), cex.lab=1.2, cex.axis=1.2)
plot_seeds(cbind(seeds_fmu$V2, seeds_ebt$V2, seeds_cm$V2), "Species 1", xlab="")
plot_seeds(cbind(seeds_fmu$V3, seeds_ebt$V3, seeds_cm$V3), "Species 2", xlab="")
plot_seeds(cbind(seeds_fmu$V4, seeds_ebt$V4, seeds_cm$V4), "Species 3", xlab="Time (years)")
# dev.off()
