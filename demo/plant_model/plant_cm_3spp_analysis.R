# 
# # code from https://traitecoevo.github.io/plant/articles/individuals.html
# 
# library(plant)
# 
#   
# # Individual plant
# ### TEST of individual plants ###
# env <- FF16_fixed_environment(1)
# y = grow_plant_to_time(FF16_Individual(s = FF16_Strategy()), seq(0,50,.2), env)
# dat=read.delim("~/codes/pspm_package/ind_plant.txt", header=F)
# dat_ref = read.delim("~/codes/plant_fate/tests/ind_plant.txt", header=F)
# 
# y = list(state=dat_ref, time=dat_ref$V1)
# colnames(y$state) = c("time", "height", "mortality", "fecundity", "area_heartwood")
# 
# par(mfrow=c(2,2), mar=c(4,4,1,1))
# plot(y$state[,"height"]~y$time, type="l")
# points(dat$V2~dat$V1, col="blue")
# plot(y$state[,"mortality"]~y$time, type="l")
# points(dat$V3~dat$V1, col="blue")
# plot(y$state[,"fecundity"]~y$time, type="l")
# points(dat$V4~dat$V1, col="blue")
# plot(y$state[,"area_heartwood"]~y$time, type="l")
# points(dat$V5~dat$V1, col="blue")
# 
# lp = read.delim("~/codes/pspm_package/light_profile_ind_plant.txt")
# 
# # points(mortality ~ time, yy, type="l", col="brown")
# # points(dat$V3~dat$V1, col="yellow3")
# 
# # ### TEST of single cohort ###
# # dat2=read.delim("/home/jaideep/codes/plant_fate/tests/cohort.txt", header=F)
# # points(dat2$V2~dat2$V1, col="cyan")
# # plot(dat2$V3~dat2$V1, col="cyan")
# # 
# # plot(dat2$V7~dat2$V1, col="green4")
# #   
# 
# # ### TEST of patch with single cohort ###
# # dat2=read.delim("/home/jaideep/codes/plant_fate_/tests/patch.txt", header=F)
# # plot(height ~ time, yy, type="l", ylim=c(0,18), main="1")
# # lines(dat2$V2~dat2$V1, col="cyan4")
# # plot(dat2$V3~dat2$V1, col="cyan")
# # 
# # plot(dat2$V7~dat2$V1, col="green4")
# 
# ### TEST of patch with multiple cohorts ###
# library(plant)
# p0 <- scm_base_parameters("FF16")
# p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p = p0, mutant = FALSE)
# p1$seed_rain = 1
# p1$control$environment_light_rescale_usually = F
# 
# ptm <- proc.time()
# data1 <- run_scm_collect(p2)
# proc.time() - ptm
#   
# t <- data1$time
# h <- data1$species[[1]]["height", , ]
# h2 <- data1$species[[2]]["height", , ]
# m <- data1$species[[1]]["mortality", , ]
# ld = data1$species[[1]]["log_density",,]
# fe = data1$species[[1]]["seeds_survival_weighted",,]
# le = data1$env
# 
# 
# matplot(t, h, lty=1, col=scales::alpha("black", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
# matlines(t, h2, lty=1, col=scales::alpha("blue", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
# 
# 
# heights=read.delim("~/codes/pspm_package_try_simple_plant/patch_full_hts.txt", header=F, col.names = paste0("V", 1:142))
# lds = read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_lds.txt", header=F, col.names = paste0("V", 1:142))
# vs = read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_vs.txt", header=F, col.names = paste0("V", 1:142))
# ms = read.delim("/home/jaideep/codes/plant_fate/tests/patch_full_ms.txt", header=F, col.names = paste0("V", 1:142))
#   
# # matplot(t, h, lty=1, col=scales::alpha("black", 0.25), type="l",
# #         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rai  n), xlim=c(0,110), ylim=c(0,18))
# 
# matplot(t, h, lty=1, col=scales::alpha("black", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain), xlim=c(0,105))
# matplot(heights$V1, heights[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Height (m)", xlim=c(0,105), add=T)
# 
# matplot(t, ld, lty=1, col=scales::alpha("black", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain), lwd=2)
# matplot(lds$V1, lds[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Log density", add=T)
# 
# 
# matplot(t, m, lty=1, col=scales::alpha("black", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Mortality (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
# matplot(ms$V1, ms[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Mortality", add=T)
# 
# matplot(t, fe, lty=1, col=scales::alpha("black", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Viable Seeds (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
# matplot(vs$V1, vs[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
#         las=1, xlab="Time (years)", ylab="Viable Seeds", add=T)
# 
# # image(z=as.matrix(heights[,-1]), x=heights[,1], xlim=c(0,105.32), zlim=c(0,20))
# # image(z=h, x=t, xlim=c(0,105.32), zlim = c(0,20))
# # 
# # image(z=log(1+as.matrix(vs[,-1])), x=vs[,1], xlim=c(0,105.32), zlim=c(0,log(1+25000)))
# # image(z=log(1+fe), x=t, xlim=c(0,105.32), zlim = c(0,log(1+25000)))
# 
# # 
# # lez=as.matrix(read.delim("/home/jaideep/codes/plant_fate/tests/le_z.txt", header=F)[-51])
# # leco=as.matrix(read.delim("/home/jaideep/codes/plant_fate/tests/le_co.txt", header=F)[-51])
# # 
# # plot(1e20, ylim=c(0,20), xlim=c(0,1))
# # for (i in 1:length(le)){
# #   lines(le[[i]][,"height"]~le[[i]][,"canopy_openness"], col=scales::alpha("black"))
# #   # lines(lez[i,]~leco[i,], col=scales::alpha("red"))
# # }
# # 
# 
# ### TEST INTERPOLATOR ###
# f1 = read.delim("/home/jaideep/codes/plant_fate/tests/interp.txt", sep=" ", header=F)
# f2 = read.delim("/home/jaideep/codes/plant_fate/tests/interp_orig.txt", sep=" ", header=F)
# 
# plot(f1$V2~f1$V1, type="l", col="red")
# points(f1$V3~f1$V1)
# points(f2$V2~f2$V1, pch=20, col="blue")
# 
# shape = 2
# mean_interval = 30
# scale = (gamma(1.0/shape)/shape/mean_interval)^shape
# p0 = shape*scale^(1.0 / shape) / gamma(1.0 / shape)

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
ld3 <- data1$species[[3]]["log_density", , ]

le = data1$env

setwd("~/codes/pspm_package/demo/plant_model/")

n = 142 # 192 # 
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
matplot(seed_rain[,1], seed_rain[,-1], type="l", lty=1, lwd=2, ylim=c(0,3500),
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
  # 
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

matlines(hp$V1, hp[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", 0))
matlines(hp2$V1, hp2[,-1], lty=1, col=scales::alpha("cyan", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(hp3$V1, hp3[,-1], lty=1, col=scales::alpha("orange", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Height (m)", main = paste0("With p1 | seed rain = ", p1$seed_rain))

matplot(t, vs, lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Viable seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, vs2, lty=1, col=scales::alpha("blue", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Viable seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, vs3, lty=1, col=scales::alpha("red", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Viable seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))

matlines(vsp$V1, vsp[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Viable seeds", main = paste0("With p1 | seed rain = ", 0))
matlines(vsp2$V1, vsp2[,-1], lty=1, col=scales::alpha("cyan", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Viable seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(vsp3$V1, vsp3[,-1], lty=1, col=scales::alpha("orange", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Viable seeds", main = paste0("With p1 | seed rain = ", p1$seed_rain))




matplot(t, exp(ld), lty=1, col=scales::alpha("black", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, exp(ld2), lty=1, col=scales::alpha("blue", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(t, exp(ld3), lty=1, col=scales::alpha("red", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain))

matlines(ldp$V1, ldp[,-1], lty=1, col=scales::alpha("red", 0.25), type="l",
        las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", 0))
matlines(ldp2$V1, ldp2[,-1], lty=1, col=scales::alpha("cyan", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain))
matlines(ldp3$V1, ldp3[,-1], lty=1, col=scales::alpha("orange", 0.25), type="l",
         las=1, xlab="Time (years)", ylab="Log density", main = paste0("With p1 | seed rain = ", p1$seed_rain))


plot(1e20, ylim=c(0,20), xlim=c(0,1))
for (i in 1:length(le)){
  lines(le[[i]][,"height"]~le[[i]][,"canopy_openness"], col=make_transparent("black"))
  # lines(lez[i,]~leco[i,], col=make_transparent("red"))
}
# x = exp(seq(log(0.01),log(18), length.out = 200))
x = seq(0,20,length.out = 200)
matplot(y = x, x= t(lep), type="l", lty=1, col= scales::alpha(rainbow(142, end = .8), 0.7), add=T)


##### FMU result 


#### Debig interpolator

a = read.delim("~/codes/pspm_package_try_simple_plant/interpolator_dump.txt", header=F, sep=" ")
plot(a$V3~a$V1, type="l", ylim=c(0.4652, 0.4656), xlim=c(0.09985, 0.0999))
points(a$V2~a$V1, type="l", col="red")
