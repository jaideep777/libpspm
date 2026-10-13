library(tidyverse)

u = read.delim("iebt_fish_x.txt", header=F)

x = seq(0, 10, length.out = 201)
t = seq(0,10, length.out = 101)


U <- (as.matrix(u[, -ncol(u)]))

image(U, x=t, y=x)

df <- data.frame(
  x = rep(x, times = length(t)),
  t = rep(t, each = length(x)),
  u = as.vector(t(U))
)

print(
ggplot(df, aes(x = t, y = x, fill = u)) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c() + 
  theme_minimal()
)


