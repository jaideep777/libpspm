library(tidyverse)

source("demo/cod_model/parameters.cod.R")
source("demo/cod_model/simulator.8_matprocess.R")

res = simulate(0,45,display=T)

df <- read.delim("iebt_fishmodel_communityprops.txt", header=F)

p1 = df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "ssb", t > 0) |>
  pivot_wider(names_from = var, values_from=value) |> 
  ggplot() +
  geom_line(aes(x=t, y=ssb/1e9), col="cyan3") + 
  geom_point(data = tibble(y=res$summaries$SSB/1e9, x=0:249)[1:100,],
             aes(x=x,y=y), col="red") + 
  labs(y="SSB (MT)", x="Years")

p1.1 = df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "tsb") |>
  pivot_wider(names_from = var, values_from=value) |> 
  ggplot() +
  geom_line(aes(x=t, y=tsb/1e9), col="cyan3") + 
  geom_point(data = tibble(y=res$summaries$TSB/1e9, x=0:249)[1:100,],
             aes(x=x,y=y), col="red") + 
  labs(y="TSB (MT)", x="Years")


matobs<-c(0,0,0.0002,0.0029,0.0565,0.3240,0.6361,0.8451,0.9473,0.9838,0.9944,0.9964,1) # maturity-at-age
  
p2 = df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "mat") |>
  filter(t > 50) |>
  ggplot() + 
  geom_line(aes(x=a, y=value, group=t, col=t)) + 
  scale_color_viridis_c(direction = -1) + 
  geom_point(data = tibble(x=1:length(matobs), y=matobs), aes(x=x, y=y), col="red")+
  labs(y="maturity", x="age") + 
  xlim(c(0,15))

p3 = df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "U_tot") |> 
  filter(t > 50, a > 0) |>
  ggplot() + 
  geom_line(aes(x=a, y=-value/1e9, group=t, col=t)) + 
  scale_color_viridis_c(direction = -1) + 
  # scale_x_log10() + 
  # scale_y_log10() + 
  geom_point(data = tibble(y=res$N[250,]/1e9, x=1:30),
             aes(x=x,y=y), col="red") + 
  labs(y="Number (billions)", x="age") 

p4 <- df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "U_tot") |> 
  filter(a > 1) |>
  ggplot() + 
  geom_raster(aes(x=t, fill=log(1-value/1e6), y=a)) + 
  scale_color_viridis_c(direction = 1) + 
  labs(y="age", x="years", fill="log(number (mn))")

p5 <- df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "mat") |> 
  filter(a > 1) |>
  ggplot() + 
  geom_raster(aes(x=t, fill=value, y=a)) + 
  scale_color_viridis_c(direction = 1) + 
  labs(y="age", x="years", fill="mat")


res$mat = rowSums(res$N*res$MAT)/rowSums(res$N)

p6 = df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var %in% c("U_tot", "mat")) |> 
  pivot_wider(names_from=var, values_from=value) |> 
  group_by(t) |> 
  filter(t > 0) |>
  summarize(mat = sum(U_tot*mat)/sum(U_tot)) |>
  ggplot() + 
  geom_line(aes(x=t, y=mat), col="cyan3") + 
  geom_point(data = tibble(y=res$mat, x=0:249)[2:100,],
             aes(x=x,y=y), col="red") + 
  labs(y="maturity", x="Years") 


library(patchwork)  

print(
  (p1 + p1.1) / (p3 + (p3+scale_x_log10()+scale_y_log10()))  / (p2 + p6) / (p4 + p5)
)

