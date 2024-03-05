# BOLD data - ISIB Project 2021 

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

# audio_BOLD_subsequence.RData
# Data structure: 80:80:33:150 4-D array
# Data type:      doubles
# Notes:          Dimensions correspond to those of BOLD.subseq data. 
#                 Only 150 time points because this data is from audio-
#                   only stimuli
load("audio_BOLD_subsequence.RData")


# visual_BOLD_subsequence.RData
# Data structure: 80:80:33:150 4-D array
# Data type:      doubles
# Notes:          Dimensions correspond to those of BOLD.subseq data. 
#                 Only 150 time points because this data is from visual-
#                   only stimuli
load("visual_BOLD_subsequence.RData")


# extract grey matter indices from mask; these will correspond to grey matter
# BOLD values in BOLD.subseq.
grey.ids <- which(mask, arr.ind = TRUE)

# extract white matter indices from mask; these will correspond to white matter
# BOLD values in BOLD.subseq
mask.inverse <- !mask
white.ids <- which(mask.inverse, arr.ind = TRUE)

# extract white matter BOLD values using indices from mask 
white.BOLD <- array(BOLD.subseq[mask.inverse], dim = c(147784, 300))

# extract white matter BOLD values from audio-reliable data
aud.white.BOLD <- array(aud.BOLD.subseq[mask.inverse], dim = c(147784, 150))

# extract white matter BOLD values from visual-reliable data
vis.white.BOLD <- array(vis.BOLD.subseq[mask.inverse], dim = c(147784, 150))

# extract grey matter BOLD values using indices from mask
grey.BOLD <- array(BOLD.subseq[mask], dim = c(63416, 300))

# extract grey matter BOLD values from audio-reliable data
aud.grey.BOLD <- array(aud.BOLD.subseq[mask], dim = c(147784, 150))

# extract grey matter BOLD values from visual-reliable data
vis.grey.BOLD <- array(vis.BOLD.subseq[mask], dim = c(147784, 150))

# All Stimuli  ---------------------------------------------------------------- 

# Step 1: Find white matter/low activity time series means.
low.means <- colMeans(white.BOLD)

# Step 2: Look at each individual series in that region + evaluate their
#         distance from the mean. How likely are they to significantly 
#         deviate from the average, even though they don't display any
#         activity? This will allow us to find the extremes in the regions 
#         with no/low activity.
#

# initialize arrays to hold mean voxel data and distances from overall white
# matter means
low.dist.array <- c(0)

# take mean of each low-activity voxel's time sequence
# measure distance from overall low-activity mean
for(i in 1:147784){
  low.dist.array[i] <- 1 - cor(low.means, white.BOLD[i,])
}

# Plot the data 
hist(low.dist.array)

# Step 3: Areas of high activity should have greater distance from the mean.
#         Compare distances with respect to the mean, and our new threshold for
#         high activity.


# compare to grey matter dist from low mean benchmark
high.dist.array <- c(0)

high.means <- colMeans(grey.BOLD)

for(i in 1:63416){
  high.dist.array[i] <- 1 - cor(low.means, grey.BOLD[i,])
}

# Plot the data 
hist(high.dist.array)

# Get confidence interval 
low.mean <- mean(low.dist.array, na.rm = TRUE)
low.mean
high.mean <- mean(high.dist.array)
high.mean
low.sd <- sd(low.dist.array, na.rm = TRUE)
low.sd
high.sd <- sd(high.dist.array)
high.sd

z <- qnorm(0.05, lower.tail = FALSE)

# upper bound for low-activity 
upper.bound.low <- (low.mean + z * low.sd)

# lower bound for high-activity
lower.bound.high <- (high.mean - z * high.sd)



# Visual Reliable Stimuli ------------------------------------------------------
low.vis.array <- c(0)
low.vis.means <- colMeans(vis.white.BOLD)

for(i in 1:147784){
  low.vis.array[i] <- 1 - cor(low.vis.means, vis.white.BOLD[i,])
}


hist(low.vis.array)

# high-activity data
high.vis.array <- c(0)
high.vis.means <- colMeans(vis.grey.BOLD)

for(i in 1:63416){
  high.vis.array[i] <- 1 - cor(low.vis.means, vis.grey.BOLD[i,])
}

hist(high.vis.array)

# thresholds 
vis.low.mean <- mean(low.vis.array, na.rm = TRUE)
vis.low.mean
vis.high.mean <- mean(high.vis.array)
vis.high.mean
vis.low.sd <- sd(low.vis.array, na.rm = TRUE)
vis.low.sd
vis.high.sd <- sd(high.vis.array)
vis.high.sd

# upper bound for low-activity 
vis.upper.bound.low <- (vis.low.mean + z * vis.low.sd)

# lower bound for high-activity
vis.lower.bound.high <- (vis.high.mean - z * vis.high.sd)



# Audio-Reliable Stimuli -------------------------------------------------------

low.aud.array <- c(0)
low.aud.means <- colMeans(aud.white.BOLD)

for(i in 1:147784){
  low.aud.array[i] <- 1 - cor(low.aud.means, aud.white.BOLD[i,])
}


hist(low.aud.array)

# high-activity data
high.aud.array <- c(0)
high.aud.means <- colMeans(aud.grey.BOLD)

for(i in 1:63416){
  high.aud.array[i] <- 1 - cor(low.aud.means, aud.grey.BOLD[i,])
}

hist(high.aud.array)

# thresholds 
aud.low.mean <- mean(low.aud.array, na.rm = TRUE)
aud.low.mean
aud.high.mean <- mean(high.aud.array)
aud.high.mean
aud.low.sd <- sd(low.aud.array, na.rm = TRUE)
aud.low.sd
aud.high.sd <- sd(high.aud.array)
aud.high.sd

# upper bound for low-activity 
aud.upper.bound.low <- (aud.low.mean + z * aud.low.sd)

# lower bound for high-activity
aud.lower.bound.high <- (aud.high.mean - z * aud.high.sd)

# Step: Map areas of activity on the brain, compare to accepted literature.







