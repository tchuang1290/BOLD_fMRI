library(ggplot2)

data(geyser, package = "MASS")
dim(geyser)

head(geyser,4)

ggplot(geyser) + 
  geom_histogram(aes(x = duration),
                 bins = 15, 
                 color = 'black', 
                 fill = 'grey')

ggplot(geyser) + 
  geom_histogram(aes(x = lag(duration)),
                 bins = 15, 
                 color = 'black', 
                 fill = 'grey')

ggplot(geyser) + 
  geom_histogram(aes(x = duration),
                 y = stat(density)
                 bins = 15, 
                 color = 'black', 
                 fill = 'grey')

d <- geyser$duration
d_short <- d[d<3]



