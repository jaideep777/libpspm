library(tidyverse)

df <- read.delim("iebt_fishmodel_communityprops.txt", header=F)

df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "ssb") |>
  pivot_wider(names_from = var, values_from=value) |> 
  ggplot(aes(x=t, y=ssb/1e9)) +
  geom_line()


mat<-c(0,0,0.0002,0.0029,0.0565,0.3240,0.6361,0.8451,0.9473,0.9838,0.9944,0.9964,1) # maturity-at-age
  
df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "mat") |>
  filter(t > 50) |>
  ggplot() + 
  geom_line(aes(x=a, y=value, group=t, col=t)) + 
  scale_color_viridis_c(direction = -1) + 
  geom_point(data = tibble(x=1:length(mat), y=mat), aes(x=x, y=y), col="red")+
  labs(y="maturity", x="age")

df |> setNames(c("t", "a", "var", "value")) |> 
  filter(var == "U_tot") |> 
  filter(t > 50, a > 1) |>
  ggplot() + 
  geom_line(aes(x=a, y=-value/1e9, group=t, col=t)) + 
  scale_color_viridis_c(direction = -1) + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(y="Number (billions)", x="age")

