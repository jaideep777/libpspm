#### Multispecies test

library(plant)
p0 <- scm_base_parameters("FF16")
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p = p0, mutant = FALSE)
p1$seed_rain = 1
p1$control$environment_light_rescale_usually = F
p2 <- expand_parameters(trait_matrix(0.2625, "lma"), p = p1, mutant = FALSE)
p3 <- expand_parameters(trait_matrix(0.4625, "lma"), p = p2, mutant = FALSE)
p2$seed_rain = c(1, 1)*1
p3$seed_rain = c(1, 1,1)*1
# p3$cohort_schedule_times
# p2$control$environment_light_rescale_usually = F

# p0 <- scm_base_parameters("FF16")
# p2 <- expand_parameters(trait_matrix(0.0825, "lma"), p = p0, mutant = FALSE)
# p2$seed_rain = 1
# p2$control$environment_light_rescale_usually = F

ptm <- proc.time()
data1 <- run_scm_collect(p3)
proc.time() - ptm

t <- data1$time
h <- data1$species[[1]]["height", , ]
h2 <- data1$species[[2]]["height", , ]
h3 <- data1$species[[3]]["height", , ]
vs <- data1$species[[1]]["seeds_survival_weighted", , ]
vs2 <- data1$species[[2]]["seeds_survival_weighted", , ]
vs3 <- data1$species[[3]]["seeds_survival_weighted", , ]
ld <- data1$species[[1]]["log_density", , ]
ld2 <- data1$species[[2]]["log_density", , ]

le = data1$light_env

setwd("~/codes/pspm_package/demo/plant_model/output_3spp_eta12/")

n = 192 # 142
hp   = read.delim("species_0_X.txt", header=F, col.names = paste0("V", 1:n))
hp2  = read.delim("species_1_X.txt", header=F, col.names = paste0("V", 1:n))
hp3  = read.delim("species_2_X.txt", header=F, col.names = paste0("V", 1:n))
vsp  = read.delim("species_0_fec.txt", header=F, col.names = paste0("V", 1:n))
vsp2 = read.delim("species_1_fec.txt", header=F, col.names = paste0("V", 1:n))
vsp3 = read.delim("species_2_fec.txt", header=F, col.names = paste0("V", 1:n))
ldp  = read.delim("species_0_u.txt", header=F, col.names = paste0("V", 1:n))
ldp2 = read.delim("species_1_u.txt", header=F, col.names = paste0("V", 1:n))
ldp3 = read.delim("species_2_u.txt", header=F, col.names = paste0("V", 1:n))

lep = read.delim("light_profile_ind_plant.txt")
lep = lep[,-ncol(lep)]

par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))


matplot(x=t(lep), y=seq(0,20,length.out=200), type="l", lty=1, col=rainbow(200),
        ylab="Height", xlab="Light level")

seed_rain = read.delim("seed_rain.txt", header=F)
seed_rain = seed_rain[-ncol(seed_rain)]
matplot(seed_rain[,1], seed_rain[,-1], type="l", lty=1, lwd=2, ylim=c(0,2500),
        ylab="Seed rain", xlab="time")

hmean = function(h,l){
  rowSums(h*exp(l), na.rm = T)/rowSums(exp(l), na.rm=T)  
}

hmean1 = hmean(hp[,-1], ldp[,-1])
hmean2 = hmean(hp2[,-1], ldp2[,-1])
hmean3 = hmean(hp3[,-1], ldp3[,-1])
hmean_all = cbind(hp[,1], hmean1, hmean2, hmean3)
matplot(hmean_all[,1], hmean_all[,-1], type="l", lty=1, lwd=2, 
        ylab="Mean height", xlab="time")

matplot(y=cbind(as.numeric(ldp[191,-1]), 
                as.numeric(ldp2[191,-1]), 
                as.numeric(ldp3[191,-1])), x=200-ldp[,1],type="l", lty=1, lwd=2,
        xlab="age", ylab="Log density")

par(mfrow=c(2,2), mar=c(4,4,1,1), oma=c(1,1,1,1))

matplot(t, h, lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, h2, lty=1, col=scales::alpha("blue", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, h3, lty=1, col=scales::alpha("red", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))

matplot(hp$V1, hp[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", 0))
matlines(hp2$V1, hp2[,-1], lty=1, col=scales::alpha("cyan", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(hp3$V1, hp3[,-1], lty=1, col=scales::alpha("orange", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))

matplot(t, vs, lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, vs2, lty=1, col=scales::alpha("blue", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, vs3, lty=1, col=scales::alpha("red", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))

matplot(vsp$V1, vsp[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", 0))
matlines(vsp2$V1, vsp2[,-1], lty=1, col=scales::alpha("cyan", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(vsp3$V1, vsp3[,-1], lty=1, col=scales::alpha("orange", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))


matplot(ldp$V1, ldp[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", 0))
matlines(ldp2$V1, ldp2[,-1], lty=1, col=scales::alpha("cyan", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(ldp3$V1, ldp3[,-1], lty=1, col=scales::alpha("orange", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))

