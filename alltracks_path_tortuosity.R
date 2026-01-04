#* ************************************************************************************************************
#* PATH TORTUOSITY
#* ************************************************************************************************************

# Aim of script: measure how much paths returning from colony deviate from straight line (direct path)
# Steps: 
# 1) Import GPS data that has trip IDs
# 2) For every trip, find the furthest transit (state 2?) location. Measure distance between that point
#    and nest. Save that value in new column
# 3) Measure the actual path length between furthest transit location and nest. Save in new column
# 4) Calculate the deviation of the actual path from the ideal direct path. 

#* ************************************************************************************************************
#* Packages
#* ************************************************************************************************************

library(dplyr)
library(seabiRds) # for getColDist()


#* ************************************************************************************************************
#* Data import
#* ************************************************************************************************************

# import data with HMM states & trip ID
df_GPS <- read.csv("data_in/allbirds_HMM_2states_bestmodel_per bird.csv", header = T, sep = ",")

df_bio <- read.csv("data_in/allbirds_bioldata.csv", header = T, sep = ",")

# merge datasets
df_full <- left_join(df_GPS, df_bio, by = c("ID" = "birdID"))

#* ************************************************************************************************************
#* 
#* ************************************************************************************************************

# Function parameters: df = dataset with GPS locations; id = name of ID column (here: bird_tripID)
# colonyLat: latitude of colony location, colonyLon = long of colony location. Individual nest location data
# could also be used if available. 
# Lat and lon columns need to be called x and y

return_path_tortuosity <- function(df, id, nestLat, nestLon, xlim, ylim, export_plots) {
  
  out_data <- data.frame(matrix(ncol = 5, nrow = length(unique(df[[id]]))))
  colnames(out_data) <- c("bird_tripID", "birdID", "dist_last_FS", "km_return_path", "tortuosity")
  
  # subset trip
  for (i in 1:length(unique(df[[id]]))) {
    in_data <- df[which(df[[id]] == unique(df[[id]])[i]),] 
    print(unique(in_data[[id]])) # see where the function is at
    
    # Reverse the entire dataset by row
    in_data_reversed <- in_data[nrow(in_data):1, ]
    
    # Find the index of the first foraging location (state 1) which is actually the last before return
    index_last_FS <- which(in_data_reversed$state == 1)[1]
    
    # Check if index_last_FS is NA and skip to the next iteration if it is
    if (is.na(index_last_FS)) {
      next
    }
    print("found foraging site")
    
    # Get distance to colony from the last loc of the last foraging site. 
    dist_last_FS <- getColDist(lat = in_data_reversed$y[index_last_FS], 
                               lon = in_data_reversed$x[index_last_FS], 
                               colonyLat = unique(in_data_reversed[[nestLat]]), 
                               colonyLon = unique(in_data_reversed[[nestLon]]))
    
    # if distance to colony is less than 1km, skip this trip
    if (dist_last_FS < 1.01) {
      next
    }
    
    print("found distance to colony")
    
    # get total distance actually traveled by bird on return path (add steps in km + dist last loc to col)
    #km_return_path <- sum(in_data_reversed$step[1:index_last_FS]) + 
    #  getColDist(lat = in_data_reversed$y[1], 
                 # lon = in_data_reversed$x[1], 
                 # colonyLat = unique(in_data_reversed[[nestLat]]), 
                 # colonyLon = unique(in_data_reversed[[nestLon]]))
    
    # get total path length
    for (j in 1:index_last_FS) {
      
      if (j == 1) {
        a <- getColDist(lat = in_data_reversed$y[j], 
                        lon = in_data_reversed$x[j], 
                        colonyLat = unique(in_data$nest_lat), 
                        colonyLon = unique(in_data$nest_lon))
        
        km_return_path <- a
        
        b <- 0
        
      }

      
      if (j > 1) {
        
        b <- getColDist(lat = in_data_reversed$y[j], 
                        lon = in_data_reversed$x[j], 
                        colonyLat = in_data_reversed$y[j-1], 
                        colonyLon = in_data_reversed$x[j-1])
      
      }
      
      km_return_path <- km_return_path + b 
      
    }
    
    # Calculate tortuosity of return path to colony
    tortuosity <- dist_last_FS / km_return_path
    
    # put results into out_data
    out_data[i, "bird_tripID"] <- unique(in_data_reversed[[id]])
    out_data[i, "birdID"] <- sub("_(.*)", "", unique(in_data_reversed[[id]]))
    out_data[i, "dist_last_FS"] <- dist_last_FS
    out_data[i, "km_return_path"] <- km_return_path
    out_data[i, "tortuosity"] <- tortuosity
    
    # User prompted plot export
    if (export_plots) {
    # export plots of foraging trips with return path highlighted
    tiff(paste("plots/trips_w_return_path/tortuosity_", unique(in_data_reversed$bird_tripID), ".tiff"))
    
    plot(in_data_reversed$x, in_data_reversed$y, col = in_data_reversed$state,
         xlim = xlim, ylim = ylim,
         main = unique(in_data_reversed$bird_tripID))
    
    # Add colony
    points(x = unique(in_data_reversed[[nestLon]]), 
           y = unique(in_data_reversed[[nestLat]]), 
           pch = 8, col= "gold", cex = 2)
    
    # Add a path between consecutive points
    lines(in_data_reversed$x, in_data_reversed$y, col = "grey20")
    
    # Add filled circles for points between the first point and index_last_FS
    points(in_data_reversed$x[1:index_last_FS], 
           in_data_reversed$y[1:index_last_FS], 
           col = in_data_reversed$state[1:index_last_FS], pch = 19)
    
    text_label <- paste("Distance to Colony:", round(dist_last_FS, 2), "km",
                        "\nTotal Distance Returned:", round(km_return_path, 2), "km",
                        "\nTortuosity:", round(tortuosity, 2))
    # Adjust y position for text placement
    usr <- par()$usr  # Get the current plotting region
    text(x = mean(usr[1:2]), 
         y = usr[4] - 0.12 * (usr[4] - usr[3]), 
         labels = text_label, 
         pos = 3, 
         cex = 0.8, 
         adj = 0.5, 
         xpd = TRUE) 
    
    dev.off()
    
    }
    
  }

  # Remove rows with NA values from out_data before returning
  out_data <- na.omit(out_data)
  
  return(out_data)
    
}

# indicate dataset, id column, and lat & lon of colony, lat & lon bounds of plot, export plot T or F
tort_res <- return_path_tortuosity(df_full, id = 19, 
                                   nestLat = 25, 
                                   nestLon = 26,
                                   xlim = c(-73.18, -74.10), 
                                   ylim = c(45.30, 46.10),
                                   export_plots = F)


# Summarize results per bird
tort_sum <- tort_res %>%
  group_by(birdID) %>%
  summarise(
    av_tort = mean(tortuosity, na.rm = TRUE),
    med_tort = median(tortuosity, na.rm = TRUE),
    sd_tort = sd(tortuosity, na.rm = TRUE),
    n = n()
  )

write.csv(tort_sum, "data_out/tortuosity_summ_per_bird.csv", row.names = F)

# Only Stefans birds:
write.csv(tort_sum[grep("GBC", tort_sum$birdID), ], "data_out/tortuosity_summ_per_bird_phytotron.csv", row.names = FALSE)


trip_measures_per_bird <- df_GPS %>% 
  group_by(ID) %>% 
  summarise(
    av_trip_dist_km = mean(trip_dist_km, na.rm = T),
    sd_trip_dist_km = sd(trip_dist_km, na.rm = T),
    av_trip_dur_h = mean(tot_trip_time_h),
    sd_trip_dur_h = sd(tot_trip_time_h), 
    n_trips = length(unique(bird_tripID))
  )

# Export
write.csv(trip_measures_per_bird, "data_out/trip_measures_per_bird.csv", row.names = F)

trip_measures_per_trip <- df_GPS %>% 
  group_by(bird_tripID) %>% 
  summarise(
    trip_dist_km = unique(trip_dist_km),
    trip_duration_h = unique(tot_trip_time_h),
    trip_date = unique(as.Date(date))
  )


# Export
write.csv(trip_measures_per_trip, "data_out/trip_measures_per_trip.csv", row.names = F)

# Visualize the probability distributions of step lengths for each bird

# Create a folder to save the plots if it doesn't already exist
if (!dir.exists("plots/step_histograms")) {
  dir.create("plots/step_histograms", recursive = TRUE)
}

for (i in 1:length(unique(df_GPS$ID))) {
  
  # Subset data for the current bird
  bird_data <- df_GPS[which(df_GPS$ID == unique(df_GPS$ID)[i]),]
 
  print(unique(bird_data$ID)) # see where the function is at
  
  # Create the histogram plot as a TIFF file
  tiff(filename = paste0("plots/step_histograms/step_histogram_", unique(bird_data$ID), ".tiff"), 
       width = 600, height = 500, res = 130)
  
  hist(bird_data$step, 
       main = paste("Step Length Distribution for ", unique(bird_data$ID)),
       xlab = "Step Length [km]",
       col = "grey",
       border = "black")
  
  dev.off()  # Close the TIFF device to save the plot
}



if (!dir.exists("plots/trip_histograms")) {
  dir.create("plots/trip_histograms", recursive = TRUE)
}


for (i in 1:length(unique(df_GPS$ID))) {
  
  # Subset data for the current bird
  bird_data <- df_GPS[which(df_GPS$ID == unique(df_GPS$ID)[i]),]
  
  print(unique(bird_data$ID)) # see where the function is at
  
  # Create the histogram plot as a TIFF file
  tiff(filename = paste0("plots/trip_histograms/trip_histogram_", unique(bird_data$ID), ".tiff"), 
       width = 600, height = 500, res = 130)
  
  hist(unique(bird_data$trip_dist_km), 
       main = paste("Trip Distance Distribution for ", unique(bird_data$ID)),
       xlab = "Trip Distance [km]",
       col = "grey",
       border = "black")
  
  dev.off()  # Close the TIFF device to save the plot
}







