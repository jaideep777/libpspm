library(tidyverse)

setwd("~/codes/libpspm/demo/")

# Error analysis
methods = c("FMU", "IFMU", "EBT", "ABM")
cols = c("purple", "green3", "cyan3", "darkgoldenrod2")

png("error_analysis_RED_Daphnia.png", height = 780*3, width = 800*3, res=300)

par(mfcol = c(2,2), mar=c(4,4,1,1), oma=c(1,1,3,1), cex.lab=1.2, cex.axis=1.2)

err_ebt = read.delim("Daphnia_model/ebt_error_analysis.txt")
err_fmu = read.delim("Daphnia_model/fmu_error_analysis.txt")
err_ifmu = read.delim("Daphnia_model/ifmu_error_analysis.txt")
err_abm = read.delim("Daphnia_model/abm_error_analysis.txt")

plot(x=1, y=NA, xlim=c(5,20000), ylim=c(1e-6, 1e-1), log="xy", ylab = "Biomass error", xlab = "Resolution")
err_fmu %>% with(points(Eb~Nf, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~Nf, type="o", col=cols[2], lwd=2))
err_ebt %>% with(points(Eb~Nf, type="o", col=cols[3], lwd=2))
err_abm %>% with(points(Eb~Nf, type="o", col=cols[4], lwd=2))
# legend(x = 6, y = 1e-4, legend = methods, col=cols, lty=1, lwd=2)
mtext("Daphnia model", line=1)

plot(x=1, y=NA, xlim=c(1,30000), ylim=c(1e-6, 1e-1), log="xy", ylab = "Biomass error", xlab = "Execution time (ms)")
err_fmu %>% with(points(Eb~tsys, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~tsys, type="o", col=cols[2], lwd=2))
err_ebt %>% with(points(Eb~tsys, type="o", col=cols[3], lwd=2))
err_abm %>% with(points(Eb~tsys, type="o", col=cols[4], lwd=2))
legend(x = 1.5, y = .5e-4, legend = methods, col=cols, lty=1, lwd=2)


## RED

err_ebt = read.delim("RED_model/ebt_error_analysis.txt")
err_fmu = read.delim("RED_model/fmu_error_analysis.txt")
err_ifmu = read.delim("RED_model/ifmu_error_analysis.txt")
err_abm = read.delim("RED_model/abm_error_analysis.txt")


plot(x=1, y=NA, xlim=c(1,50000), ylim=c(.5, 150000)/13000, log="xy", ylab = "Biomass relative error", xlab = "Resolution")
err_fmu %>% with(points(Eb~Nf, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~Nf, type="o", col=cols[2], lwd=2))
err_ebt %>% with(points(Eb~Nf, type="o", col=cols[3], lwd=2))
err_abm %>% with(points(Eb~Nf, type="o", col=cols[4], lwd=2))
# legend(x = 1, y = 1e-2, legend = methods, col=cols, lty=1, lwd=2)
mtext("RED model", line=1)

plot(x=1, y=NA, xlim=c(10,100000), ylim=c(.5, 150000)/13000, log="xy", ylab = "Biomass relative error", xlab = "Execution time (ms)")
err_fmu %>% with(points(Eb~tsys, type="o", col=cols[1], lwd=2))
err_ifmu %>% with(points(Eb~tsys, type="o", col=cols[2], lwd=2))
err_ebt %>% with(points(Eb~tsys, type="o", col=cols[3], lwd=2))
err_abm %>% with(points(Eb~tsys, type="o", col=cols[4], lwd=2))

dev.off()

