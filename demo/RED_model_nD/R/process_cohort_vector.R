install.packages("ggplot2")
library(ggplot2)

data <- read.csv("cohort_vector.txt", header = FALSE, col.names = c("birthtime", "X", "Y", "Z", "spName"))
data <- data[1:100,]

ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_x_continuous(trans = "log")
  