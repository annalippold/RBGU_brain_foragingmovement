#****************************************************************************************************
# ALL BIRDS SPATIAL ECOLOGY -- Hidden Markov Models for foraging locations
#****************************************************************************************************

# Aim of script: find the number of feeding locations per feeding trip

#****************************************************************************************************
# Packages needed

#HMM
library(moveHMM)
library(data.table) # for rbindlist()
library(dplyr)
library(parallel) # for paralell computations (for-loop for maximum likelihood optimization)
library(adehabitatHR) # for utilisation distribution of feeding sites
library(sp) # transforming data into SpatialPoints

#* **************************************************************************************************
#* Data import
#* **************************************************************************************************

# Use dataHMM1 -> this is rediscreditized (myTraj_redis as dataframe), and colony is removed, 
# and trips have IDs, plus outlier and short trips (<=60min) are removed

#* **************************************************************************************************
#* Data preparation
#* **************************************************************************************************


# Certain columns have to have the right names for the functions in moveHMM to work. Rename.
# ID: an identifier for the track/individual; -> here: TRACK ID!!!!
# x: Easting or longitude coordinate;
# y: Northing or latitude coordinate.
# step: distance between 2 locations
# angle: angle between 2 steps

# Step and angle will be calculated by preData(), so remove them here to avoid duplicates
dataHMM2 <- dataHMM1 %>% 
  select(-c("dist", "rel.angle"))

dataHMM2 <- rename(dataHMM2, "ID"="birdID")

dataHMM2 <- prepData(dataHMM2) # this function prepares the data to work for moveHMM
summary(dataHMM2)

#**************************************************************************************************
# Fitting HMM *************************************************************************************
#**************************************************************************************************

#* 2 state model **********************************************************************************

# Generate random starting values based on a set of plausible starting values. Use histograms
# to find plausible starting values: 

# Plot the step distribution to make an "educated guess" on how to select the initial parameters
hist(dataHMM2$step)
# Are there steps of length 0?
whichzero <- which(dataHMM2$step == 0)
length(whichzero)/nrow(dataHMM) # 0.00002775388 % of steps are of 0 km length
min(dataHMM2$step, na.rm=T); max(dataHMM2$step, na.rm=T) # 0 - 16.81905
median(dataHMM2$step, na.rm=T) # 0.5231469
mean(dataHMM2$step, na.rm=T) # 1.505426
# Plot histogramm of turning angles
hist(dataHMM2$angle, breaks = seq(-pi, pi, length = 15))

#****************************************************************************************************
# Finding best model (best starting values)

# Create cluster of size ncores
ncores <- detectCores() - 1
cl <- makeCluster(getOption("cl.cores", ncores))

# Export objects needed in parallelised function to cluster
clusterExport(cl, list("dataHMM2", "fitHMM"))

# Number of tries with different starting values
niter <- 25

# Create list of starting values
allPar2 <- lapply(as.list(1:niter), function(x) {
  # Step length mean
  stepMean0 <- runif(2,
                     min = c(0.0001, 2),
                     max = c(0.007, 17))
  # Step length standard deviation
  stepSD0 <- runif(2,
                   min = c(0.0001, 1),
                   max = c(0.005, 5))
  # Proportion of steps with length = 0
  zeromass0 <- runif(2,
                     min = c(0, 0),
                     max = c(0.002, 0))
  
  # Turning angle mean
  angleMean0 <- c(0, 0)
  # Turning angle concentration
  angleCon0 <- runif(2,
                     min = c(0.01, 2),
                     max = c(0.1, 3))
  # Return vectors of starting values
  stepPar0 <- c(stepMean0, stepSD0, zeromass0)
  anglePar0 <- c(angleMean0, angleCon0)
  return(list(step = stepPar0, angle = anglePar0))
})

# Fit the niter models in parallel
allm2_parallel <- parLapply(cl = cl, X = allPar2, fun = function(par0) {
  m2 <- fitHMM(data = dataHMM2, nbStates = 2, stepPar0 = par0$step,
               anglePar0 = par0$angle)
  return(m2)
})


allm2_p_unlist <- unlist(lapply(allm2_parallel, function(m2) m2$mod$minimum))
allm2_p_unlist # rel. numerically stable


whichbest2 <- which.min(allm2_p_unlist)
# Best fitting model
mbest2 <- allm2_parallel[[whichbest2]]
mbest2

#Value of the maximum log-likelihood: -89886.76 
# Step length parameters:
#   ----------------------
#                state 1           state 2
# mean      0.19755194093 2.936646752891027
# sd        0.25635825319 2.243352187991005
# zero-mass 0.00005428633 0.000000009999658
# 
# Turning angle parameters:
#   ------------------------
#               state 1      state 2
# mean          0.05651358 -0.004791876
# concentration 0.10835964  1.553398264
# 
# Regression coeffs for the transition probabilities:
#   --------------------------------------------------
#            1 -> 2    2 -> 1
# intercept -1.654083 -1.552185
# 
# Transition probability matrix:
#   -----------------------------
#       [,1]      [,2]
# [1,] 0.8394422 0.1605578
# [2,] 0.1747710 0.8252290
# 
# Initial distribution:
#   --------------------
#   [1] 0.3325288 0.6674712

# Export tracks with 2 states for each bird
for (i in 1:length(unique(mbest2[[1]]$ID))) {
  tiff(paste("plots/allbirds_HMM tracks 2 states_per bird/", unique(mbest2[[1]]$ID)[i], ".tiff"))
  plot(mbest2, ask = F, animals = i)
  dev.off()
}


#****************************************************************************************************
#* Extract states for each observation
#* **************************************************************************************************

#* Add the state to each observation ***************************************************************

out <- mbest2$data
out$state <- viterbi(mbest2)

output <- out %>%
  left_join(dataHMM1, by = c('ID'='birdID', 'x', 'y')) #%>%
#dplyr::select(ID, step, angle, x, y, bird_loggerID, date,  state, time_diff, tripID,
#             tot_trip_time_min, tot_trip_time_h, trip_dist_km)


#save output

write.csv(output, "data_out/allbirds_HMM_2states_bestmodel_per bird.csv", row.names=FALSE)


#* ***************************************************************************************************
#* Calculate number of feeding sites
#* ***************************************************************************************************

#* Calculate nr of feeding sites per bird ****************************************************

# Remove all the trasitory locations, leaving only the feeding location points (state 1)
# State 1 characteristics: mean step length 0.33 km, mean turning angle 3.13

n_feeding <- output %>% 
  filter(state == 1) # state 1: feeding, state 2: searching

# Export that dataset to use for habitat presence script
write.csv(n_feeding, "data_out/allbirds_HMM_feeding state_bestmodel_perbird.csv")

# Remove foraging sites that have only one location

rmv_1loc_FS <- function(data, date_col) {
  # Calculate time differences (gaps) before and after each point
  gap_before <- c(Inf, difftime(data[[date_col]][-1], data[[date_col]][-nrow(data)], units = "mins"))
  gap_after <- c(difftime(data[[date_col]][-1], data[[date_col]][-nrow(data)], units = "mins"), Inf)
  
  # Identify points with large gaps both before and after
  is_single_point <- gap_before > 10 & gap_after > 10
  
  # Keep only points that are not single points
  return(data[!is_single_point, ])
}


n_feeding_filtered <- rmv_1loc_FS(n_feeding, "date")


# For kernelUD() to work, data needs to be of class SpatialPoints
n_feeding_WGS84 <- SpatialPointsDataFrame(coords=as.data.frame(cbind(n_feeding_filtered$x, 
                                                                     n_feeding_filtered$y)), 
                                          data= as.data.frame(n_feeding_filtered), 
                                          proj4string = CRS("+proj=longlat +datum=WGS84"))

n_feeding_mtm8 <- spTransform(n_feeding_WGS84, 
                              CRS("+proj=tmerc +lat_0=0 +lon_0=-73.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 
                            +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))

# Run utilisation distribution on all locations of state 1 per bird

# The function kernelUD of the package adehabitatHR estimates the UD. The UD is estimated in each pixel 
# of a grid superposed to the relocations. (adehabitatHR manual)
# FS = feeding sites
# h=LSCV -> least square cross validation method. Makes more sense for feeding locations than 
# h=href, which I used for home range 


FS_kern <- kernelUD(n_feeding_mtm8[1], h="LSCV", grid = 400, extent = 0.5) # indicate ID column and h (smoothing par). 
FS_kern
image(kernelUD(n_feeding_mtm8[1], grid=300, extent=0.1)) # all birds one image

# get area of FS
FSarea95_href <- getverticeshr(FS_kern, percent = 95)

for (i in 1:length(FSarea95_href)) {
  tiff(paste("plots/allbirds_ind feeding sites UD hLSCV 95_2states/UD95_", FSarea95_href[[1]][i], ".tiff"))
  plot(FSarea95_href[i,], col = "darkred", xlim = c(320000, 280000), ylim = c(5030000, 5102728))
  #points(n_feeding_mtm8@coords)
  title(main = paste("95% UD kernel", FSarea95_href[[1]][i], 
                     "\ntotal area:", round(FSarea95_href[[2]][i]/100, digits = 2), "km2"))
  points(x = 309449, y = 5063782, pch = 23, bg = "chartreuse", cex = 2) # Add marker for colony
  dev.off()
}

# Images without colony or any writing to make automatic counting easier
for (i in 1:length(FSarea95_href)) {
  # Use trimws() to remove any leading/trailing spaces from FSarea95_href[[1]][i]
  file_name <- paste0("plots/for_analysis/allbirds_indFS_UD_LSCV_95_2states/UD95_", trimws(FSarea95_href[[1]][i]), ".tiff")
  
  tiff(file_name)
  
  plot(FSarea95_href[i,], col = "black")
  # points(n_feeding_mtm8@coords)
  dev.off()
}


#* ******************************************************************************************************
#* Count nr of foraging sites
#* ******************************************************************************************************

# base code from chatGPT 4, with own improvements (Sept 17 2024)

# Install EBImage 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EBImage")
library(EBImage)

# Test on one image before making a loop, also display to make sure it worked --------------------------

testimage <- readImage("plots/for_analysis/allbirds_indFS_UD_LSCV_95_2states/UD95_21LD24.tiff")
gray_img <- channel(testimage, "gray")

# Enhance contrast (optional)
gray_img <- normalize(gray_img)

# Apply a threshold to binarize the image (black particles)
threshold_value <- 0.55  # Adjust this value if needed
binary_img <- gray_img < threshold_value  # Black particles will be white (1)

# Apply morphological operations to enhance particle detection
binary_img <- dilate(binary_img, makeBrush(5, shape = 'disc'))  # Adjust brush size if needed
binary_img <- erode(binary_img, makeBrush(3, shape = 'disc'))   # Adjust brush size if needed

# Label the binary image
labeled_img <- bwlabel(binary_img)

# Count particles
particle_count <- max(labeled_img) # Subtract 1 if the background label is 0
# Show the original and processed images
display(gray_img, method = "raster")
display(binary_img, method = "raster")
display(colorLabels(labeled_img), method = "raster")

# Output the number of particles
cat("Number of particles:", particle_count, "\n")

#----------------------------------------------------------------------------------------------------


# Define the folder containing the TIFF files
folder_path <- "plots/for_analysis/allbirds_indFS_UD_LSCV_95_2states"

# List all TIFF files in the folder
file_list <- list.files(folder_path, pattern = "\\.tiff$", full.names = TRUE)

# Initialize a dataframe to store results
FS_results <- data.frame(birdID = character(), n_FS_raw = numeric(), stringsAsFactors = FALSE)

# Loop through each file
for (file_path in file_list) {
  # Read and process each image
  testimage <- readImage(file_path)
  gray_img <- channel(testimage, "gray")
  
  # Enhance contrast (optional)
  gray_img <- normalize(gray_img)
  
  # Apply a threshold to binarize the image (black particles)
  threshold_value <- 0.55  # Adjust this value if needed
  binary_img <- gray_img < threshold_value  # Black particles will be white (1)
  
  # Apply morphological operations to enhance particle detection
  binary_img <- dilate(binary_img, makeBrush(5, shape = 'disc'))  # Adjust brush size if needed
  binary_img <- erode(binary_img, makeBrush(3, shape = 'disc'))   # Adjust brush size if needed
  
  # Label the binary image
  labeled_img <- bwlabel(binary_img)
  
   # Count particles
  particle_count <- max(labeled_img)  # Count distinct labels in the labeled image
  
  # Extract birdID from the file name
  file_name <- basename(file_path)
  birdID <- sub("UD95_(\\d{2}[A-Z]+\\d{2})\\.tiff", "\\1", file_name)
  
  # Append results to the dataframe
  FS_results <- rbind(FS_results, data.frame(birdID = birdID, n_FS_raw = particle_count, stringsAsFactors = FALSE))
}


write.csv(FS_results, "data_out/allbirds_nFS_results_raw.csv", row.names = F)


remove(background_label, birdID, binary_img, gray_img, labeled_img, orgimage, particle_count, 
       particle_labels, testimage, threshold_value)



