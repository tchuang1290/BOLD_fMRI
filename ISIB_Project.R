# BOLD data - ISIB Project 2021 
setwd("Documents/ISIB")
library(dplyr)
library(mosaic)
library(ggplot2)

# BOLD_subsequence.RData
# Data structure: 80:80:33:300 4D array
# Data type:      doubles
# Notes:          First 3 dimensions == voxel's relative position
#                 Fourth dimension == ea. voxel's 300-pt time series
#                 Data does not include washout points 
load("BOLD_subsequence.RData")

# mask.RData
# Data structure: 80:80:33 3D array
# Data type:      logical/boolean
# Notes:          ea. elt corresponds to the elt in BOLD_subsequence.RData with
#                   the same indices. 
#                 mask created by calculating mean of ea. original time series
#                   means > 70th percentile -> TRUE
#                   means < 70th percentile -> FALSE 
load("mask.RData")

# extract grey matter indices from mask; these will correspond to grey matter
# BOLD values in BOLD.subseq.
grey.ids <- which(mask, arr.ind = TRUE)

# extract white matter indices from mask; these will correspond to white matter
# BOLD values in BOLD.subseq
mask.inverse <- !mask
white.ids <- which(mask.inverse, arr.ind = TRUE)

# extract white matter BOLD values using indices from mask 
# TODO: Fix this 
# white.BOLD <- array(BOLD.subseq[as.logical(white.ids)], dim = c(147784, 300))
# mean(white.BOLD)
white.BOLD <- array(BOLD.subseq[mask.inverse], dim = c(147784, 300))
mean(white.BOLD)

# extract grey matter BOLD values using indices from mask
# grey.BOLD <- array(BOLD.subseq[as.logical(grey.ids)], dim = c(63416, 300))
# mean(grey.BOLD)
grey.BOLD <- array(BOLD.subseq[mask], dim = c(63416, 300))
high.mean <- mean(grey.BOLD)
grey.se <- sd(grey.BOLD)

# Step 1: Find overall mean 
overall.mean <- mean(BOLD.subseq)
overall.mean
overall.sd <- sd(BOLD.subseq)

z <- qnorm(0.025, lower.tail = FALSE)
z

(ci <- c(overall.mean - z*overall.sd, overall.mean + z*overall.sd))

# Step 2: Find white matter/low activity mean. This is our baseline.
low.mean <- mean(white.BOLD)
white.se <- sd(white.BOLD)

# Step 3: Look at each individual series in that region + evaluate their
#         distance from the mean. How likely are they to significantly 
#         deviate from the average, even though they don't display any
#         activity? This will allow us to find the extremes in the regions 
#         with no/low activity.
#           Construct a confidence interval for the range of likely 
#             "low-activity" values and then find the upper and lower bounds
#
#         3.1 Make an array of distances from mean 
#         3.2 Construct confidence interval from this array for likely distance
#              - do we actually need this? 

# initialize arrays to hold mean voxel data and distances from overall white
# matter means
low.mean.array <- c(0)
dist.array <- c(0)
dist.array2 <- c(0)

# take mean of each low-activity voxel's time sequence
# measure distance from overall low-activity mean
for(i in 1:147784){
  low.mean.array[i] <- mean(white.BOLD[i,])
  dist.array[i] <- abs((low.mean - low.mean.array[i]))
  # dist.array2[i] <- abs((overall.mean-low.mean.array[i]))
}
low.mean.array
dist.array
mean(low.mean.array)
mean(dist.array)
# mean(dist.array2)

z1 <- qnorm(.05, lower.tail = FALSE)

(ci <- low.mean + z1*white.se)

# (ci <- c(overall.mean - z*white.se, overall.mean + z*white.se))

# (ci <- c(mean(dist.array) - z*white.se, mean(dist.array) + z*white.se))

#orig.sample <- data.frame(y = low.mean.array)


#boot_data <- mosaic::do(100)*(
  #orig.sample %>%
    #sample_frac(replace = TRUE) %>%
    #count(y) %>%
    #mutate(prop = n/sum(n))
#)

#boot_data %>%
  #ggplot(aes(x = prop)) + 
  #geom_histogram()
  #theme_classic()


# Step 4: Areas of high activity should have greater distance from the mean.
#         Compare distances with respect to the mean, and our new threshold for
#         high activity.


# compare to grey matter dist from low mean benchmark 
high.mean.array <- c(0)
hldist.array <- c(0)
#hldist.array2 <- c(0)

for(i in 1:63416){
  high.mean.array[i] <- mean(grey.BOLD[i,])
  hldist.array[i] <- abs((low.mean - high.mean.array[i]))
  # hldist.array2[i] <- abs((overall.mean - high.mean.array[i]))
}
high.mean.array
hldist.array
mean(high.mean.array)
mean(hldist.array)
# mean(hldist.array2)

(ci <- c(high.mean - z*grey.se, high.mean + z*grey.se))

# Create histogram to visualize results 
# TODO: fine-tune visual
hist(hldist.array, 
     xlab = "Distance from Low-Activity Benchmark in BOLD AUs", 
     ylab = "Number of Voxels", 
     main = "Distribution of High-Activity Voxel Distance from Benchmark", 
     col = "red")

hist(high.mean.array,
     xlab = "Average Activity",
     ylab = "Number of Voxels",
     main = "Distribution of High-Activity Voxels",
     col = "green")

hist(low.mean.array,
     xlab = "Average Activity",
     ylab = "Number of Voxels",
     main = "Distribution of Low-Activity Voxels",
     col = "blue")

plot(low.mean.array)


# Step 5: Use Pearson's distance to see if we can retrieve the areas 
#         of activity based on the above algorithm. Should get correlation 
#         between expected ROIs and the high activity during stimuli.

#         Create an array of all mean voxel reads. Compare all mean voxel reads
#         to the described TRUE/FALSE vals. (Won't this recreate 
#         the data we already have?) 


# Step 6: Map areas of activity on the brain, compare to accepted literature.

require("grDevices") # for colours
x <- y <- seq(-4*pi, 4*pi, length.out = 27)
r <- sqrt(outer(x^2, y^2, "+"))
image(z = z <- cos(r^2)*exp(-r/6), col = gray.colors(33))
image(z, axes = FALSE, main = "Math can be beautiful ...",
      xlab = expression(cos(r^2) * e^{-r/6}))
contour(z, add = TRUE, drawlabels = FALSE)

# Volcano data visualized as matrix. Need to transpose and flip
# matrix horizontally.
image(t(volcano)[ncol(volcano):1,])

# A prettier display of the volcano
x <- 10*(1:nrow(volcano))
y <- 10*(1:ncol(volcano))
image(x, y, volcano, col = hcl.colors(100, "terrain"), axes = FALSE)
contour(x, y, volcano, levels = seq(90, 200, by = 5),
        add = TRUE, col = "brown")
axis(1, at = seq(100, 800, by = 100))
axis(2, at = seq(100, 600, by = 100))
box()
title(main = "Maunga Whau Volcano", font.main = 4)




