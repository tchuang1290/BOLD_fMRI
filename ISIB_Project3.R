# BOLD data - ISIB Project 2021 
# Written by: Megan Gelement and Ting Huang


# FUNCTION DEFINITIONS----------------------------------------------------------

# 2D ANALYSIS FUNCTIONS---------------------------------------------------------

# get_distances 
# Purpose:      Calculate Pearson's Distance (1 - correlation) between each 
#                 voxel read in a time series and that voxel's corresponding 
#                 benchmark
# Parameters:   benchmarks: list of benchmark values
#               data:       BOLD time series array containing each voxel's  
#                             entire time series
#               num.voxels: integer number of voxels to calculate distance for
# Returns:      dist.vector: vector containing all Pearson's distances.
# Notes:        TODO use apply function instead of for loop?
# 
get_distances <- function(benchmarks, data, num.voxels){
  dist.vector <- c(0)
  for(i in 1:num.voxels){
    dist.vector[i] <- 1 - cor(benchmarks, data[i,])
  }
  return(dist.vector)
}


# build_histogram
# Purpose:        Generate a histogram displaying Pearson's distance data   
# Parameters:     data: array of data to display
#                 xlabel: string label of x-axis
#                 ylabel: string label of y-axis
#                 title:  string title for histogram
#                 upper.bound: numerical upper bound of confidence interval
#                 lower.bound: numerical lower bound of confidence interval
# Returns:        None
# Notes:          TODO
build_histogram <- function(data, xlabel, ylabel, title, upper.bound, 
                            lower.bound){
  hist(data,
       xlab = xlabel, ylab = ylabel, main = title, col = "grey")
  abline(v=mean(data, na.rm = TRUE), col = "red", lwd = 2)
  abline(v=quantile(data, .95, na.rm = TRUE), col = "blue", lwd = 2)
  abline(v=quantile(data, .05, na.rm = TRUE), col = "blue", lwd = 2)
  abline(v=upper.bound, col = "green", lwd = 2)
  abline(v=lower.bound, col = "green", lwd = 2)
  legend("topright", 
         c("Mean", "95% Confidence Interval", "5th and 95th Percentiles"),
         col = c("red", "green", "blue"), lwd = 3, xpd = TRUE, inset = c(0,0),
         cex = 0.4)
}


remove_zeroes <- function(data, num.a = 80, num.b = 80, num.c = 33, num.d){
  for(i in 1:num.a){
    for(j in 1:num.b){
      for(k in 1:num.c){
        for(l in 1:num.d){
          if(data[i,j,k,l] == 0){
            data[i,j,k,l] = NA
          }
        }
      }
    }
  }
  return(data)
}



# 3D ANALYSIS FUNCTIONS---------------------------------------------------------

# is_all_zeroes
# Purpose:      Determine whether a given list contains all zeroes
# Parameters:   List to check for zeroes
# Returns:      Boolean (TRUE if all zeroes; FALSE otherwise)
# Notes:        Called by get_3D_distances; eliminates warnings stemming from
#                 calling cor() on list of zeroes
is_all_zeroes <- function(data){
  for(i in 1:length(data)){
    if(is.na(data[i])){
      next
    }
    else if(data[i] != 0){
      return(FALSE)
    }
  }
  return(TRUE)
}


# get_3D_distances
# Purpose:          Calculate Pearson's distance (1 - correlation) between
#                     each voxel's time series and the benchmark while
#                     retaining its 3D position.
# Parameters:       benchmarks: list of low-activity benchmarks
#                   data:       4D array of time series
#                   x, y, z:    first 3 dimensions (dimensions by time series)
# Returns:          3D vector of Pearson's distances
# Notes:            Not used yet. 
get_3D_distances <- function(benchmarks, data, x = 80, y = 80, z = 33){
  dist.vector <- array(rep(NA, x*y*z), dim=c(x, y, z))
  for(i in 1:x){
    for(j in 1:y){
      for(k in 1:z){
        if(!is_all_zeroes(data[i,j,k,])){
          dist.vector[i,j,k] <- 1 - cor(benchmarks, data[i,j,k,])
        }
      }
    }
  }
  return(dist.vector)
}


# apply_mask
# Purpose:     Retrieve data where mask value == TRUE
# Parameters: mask: a 3D boolean array
#             data: the 4D array on which to apply the mask
# Returns:    4D array containing only desired values, NAs elsewhere
# Notes:      Used to extract grey vs. white matter from overall dataset in 3D 
apply_mask <- function(mask, data){
  dims <- dim(data)
  new.mask <- array(rep(NA, dims[1]*dims[2]*dims[3]*dims[4]), dim = dims)
  for(i in 1:dims[1]){
    for(j in 1:dims[2]){
      for(k in 1:dims[3]){
        if(mask[i, j, k] == TRUE){
          new.mask[i,j,k,] = data[i,j,k,]
        }
      }
    }
  }
  
  return(new.mask)
}


# build_mask
# Purpose:    Create a 3D logical mask
# Parameters: TODO
# Returns:    3D boolean array
build_mask <- function(data, cutoff){
  dims <- dim(data)
  new.mask <- array(rep(NA, dims[1]*dims[2]*dims[3]), dim = dims)
  for(i in 1:dims[1]){
    for(j in 1:dims[2]){
      for(k in 1:dims[3]){
        if(!(is.na(data[i, j, k])) && (data[i,j,k] >= cutoff)){
          new.mask[i,j,k] = TRUE
        } else {
          new.mask[i,j,k] = FALSE
        }
      }
    }
  }
  
  return(new.mask)
}



# LOAD IN DATA -----------------------------------------------------------------


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




# PARTITION GREY AND WHITE MATTER-----------------------------------------------

# extract grey matter indices from mask
grey.ids <- which(mask, arr.ind = TRUE)

# extract white matter indices from mask
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

all.BOLD <- array(BOLD.subseq[na.rm = T], dim = c(211200, 300))



# CALCULATE AND DISPLAY DISTANCES: ALL STIMULI----------------------------------

# get benchmarks
low.benchmarks <- colMeans(white.BOLD, na.rm = TRUE)

# calculate white matter/low activity distances
low.dist.array <- get_distances(low.benchmarks, white.BOLD, 147784)

# calculate white matter/low activity confidence interval
low.mean <- mean(low.dist.array, na.rm = TRUE)
low.sd <- sd(low.dist.array, na.rm = TRUE)
z_alpha <- qnorm(0.025, lower.tail = FALSE)
ci.low.all <- c((low.mean - z_alpha * low.sd), (low.mean + z_alpha * low.sd))

# display white matter/low activity data as histogram
build_histogram(low.dist.array, "Distance from Benchmark (1 - correlation)",
                "No. Voxels", "Low-Activity Pearson's Distance, All Stimuli",
                ci.low.all[1], ci.low.all[2])

# calculate grey matter/high activity distances
high.dist.array <- get_distances(low.benchmarks, grey.BOLD, 63416)

# calculate grey matter/high activity confidence interval
high.mean <- mean(high.dist.array)
high.sd <- sd(high.dist.array)
ci.high.all <- c((high.mean - z_alpha * high.sd), 
                 (high.mean + z_alpha * high.sd))

# display grey matter/high activity data as histogram
build_histogram(high.dist.array, "Distance from Benchmark (1 - correlation)", 
                "No. Voxels", "High-Activity Pearson's Distance, All Stimuli",
                ci.high.all[1], ci.high.all[2])


# Graphs for all stimuli with benchmark

dist.array <- get_distances(low.benchmarks, all.BOLD, 211200)

# calculate grey matter/high activity confidence interval
all.mean <- mean(dist.array)
all.sd <- sd(dist.array)
ci.all <- c((all.mean - z_alpha * all.sd), (all.mean + z_alpha * all.sd))

# display grey matter/high activity data as histogram
build_histogram(dist.array, "Distance from Benchmark (1 - correlation)", 
                "No. Voxels", "All-Activity Pearson's Distance, All Stimuli",
                ci.low.all[1], ci.low.all[2])

BOLD.subseq1 <- remove_zeroes(BOLD.subseq, num.d = 300)
all.BOLD1 <- array(BOLD.subseq1, dim = c(211200, 300))

hist(all.BOLD1,
     xlab = "BOLD Values", 
     ylab = "No. Voxels", 
     main = "All-Activity, All Stimuli", 
     col = "grey")
abline(v=mean(all.BOLD, na.rm = TRUE), col = "red", lwd = 2)
abline(v=quantile(all.BOLD, .7, na.rm = TRUE), col = "green", lwd = 2)
legend("topright", 
       c("Mean", "70th Percentile"),
       col = c("red", "green"), lwd = 3, xpd = TRUE, inset = c(0,0),
       cex = 0.4)




# CALCULATE AND DISPLAY DISTANCES: VISUAL-RELIABLE STIMULI----------------------

# get benchmarks
low.vis.benchmarks <- colMeans(vis.white.BOLD, na.rm = TRUE)

# calculate white matter/low activity distances
low.vis.array <- get_distances(low.vis.benchmarks, vis.white.BOLD, 147784)

# calculate white matter/low activity confidence interval
low.vis.mean <- mean(low.vis.array, na.rm = TRUE)
low.vis.sd <- sd(low.vis.array, na.rm = TRUE)
ci.low.vis <- c((low.vis.mean - z_alpha * low.vis.sd), 
                (low.vis.mean + z_alpha * low.vis.sd))

# display white matter/low activity data as histogram
build_histogram(low.vis.array, "Distance from Benchmark (1 - correlation)", 
                "No. Voxels", 
                "Low-Activity Pearson's Distance, Visual-Reliable Stimuli",
                ci.low.vis[1], ci.low.vis[2])

# calculate grey matter/high activity distances
high.vis.array <- get_distances(low.vis.benchmarks, vis.grey.BOLD, 63416)

# calculate grey matter/high activity confidence interval
high.vis.mean <- mean(high.vis.array, na.rm = TRUE)
high.vis.sd <- sd(high.vis.array, na.rm = TRUE)
ci.high.vis <- c((high.vis.mean - z_alpha * high.vis.sd), 
                 (high.vis.mean + z_alpha * high.vis.sd))

# display grey matter/high activity data as histogram
build_histogram(high.vis.array, "Distance from Benchmark (1 - correlation)", 
                "No. Voxels", 
                "High-Activity Pearson's Distance, Visual-Reliable Stimuli",
                ci.high.vis[1], ci.high.vis[2])




# CALCULATE AND DISPLAY DISTANCES: AUDITORY-RELIABLE STIMULI--------------------

# get benchmarks
low.aud.benchmarks <- colMeans(aud.white.BOLD, na.rm = TRUE)

# calculate white matter/low activity distances
low.aud.array <- get_distances(low.aud.benchmarks, aud.white.BOLD, 147784)

# calculate white matter/low activity confidence interval
low.aud.mean <- mean(low.aud.array, na.rm = TRUE)
low.aud.sd <- sd(low.aud.array, na.rm = TRUE)
ci.low.aud <- c((low.aud.mean - z_alpha * low.aud.sd), 
                (low.aud.mean + z_alpha * low.aud.sd))

# display white matter/low activity data as histogram
build_histogram(low.aud.array, "Distance from Benchmark (1 - correlation)", 
                "No. Voxels", 
                "Low-Activity Pearson's Distance, Auditory-Reliable Stimuli",
                ci.low.aud[1], ci.low.aud[2])

# calculate grey matter/high activity distances
high.aud.array <- get_distances(low.aud.benchmarks, aud.grey.BOLD, 63416)

# calculate grey matter/high activity confidence interval
high.aud.mean <- mean(high.aud.array, na.rm = TRUE)
high.aud.sd <- sd(high.aud.array, na.rm = TRUE)
ci.high.aud <- c((high.aud.mean - z_alpha * high.aud.sd), 
                 (high.aud.mean + z_alpha * high.aud.sd))

# display grey matter/high activity data as histogram
build_histogram(high.aud.array, "Distance from Benchmark (1 - correlation)", 
                "No. Voxels", 
                "High-Activity Pearson's Distance, Auditory-Reliable Stimuli",
                ci.high.aud[1], ci.high.aud[2])



# 2D/3D VISUALIZATION-----------------------------------------------------------

# retrieve BOLD data by matter type (in 3D) 
high.aud.array.3D <- apply_mask(mask, aud.BOLD.subseq)
low.aud.array.3D <- apply_mask(!mask, aud.BOLD.subseq)

# plot slices
image(x = c(1:80), y = c(1:80), z = high.aud.array.3D[,,22, 22])
image(x = c(1:80), y = c(1:80), z = low.aud.array.3D[,,22, 22])
image(x = c(1:80), y = c(1:80), z = aud.BOLD.subseq[,,22, 22])

# get distances in 3D 
high.aud.3D.dist <- get_3D_distances(low.aud.benchmarks, high.aud.array.3D)

# set cutoff: TODO: check to make sure quantile is getting all the data
high.aud.cutoff <- quantile(high.aud.3D.dist, 0.95, na.rm = TRUE)

# partition out and plot high-activity data only
high.aud.mask <- build_mask(high.aud.3D.dist, high.aud.cutoff)

# retrieve BOLD data by cutoff
cutoff.aud.3D <- apply_mask(high.aud.mask, aud.BOLD.subseq)

# plot slice of cutoff data
image(x = c(1:80), y = c(1:80), z = cutoff.aud.3D[,,22, 22])
