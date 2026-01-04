# Remove foraging trips that have artificial sections due to GPS gaps (threshold can be defined in function call)


filter_robot_trips <- function(df, id, threshold) {
  
  library(dplyr)
  
  # Convert id to a symbol
  id_col <- ensym(id)
  print(paste0("n total trips: ", length(unique(df[[as.character(id_col)]]))))
  
  # Identify trips with "robotic" sections
  robot_trips <- df %>%
    filter(state == 2) %>% # filter out only transit locations
    group_by(!!id_col) %>%
    filter(
      any(rle(round(angle, 5))$lengths >= threshold) &  # Check for repeating rounded angles
        any(rle(round(step, 4))$lengths >= threshold)    # Check for repeating rounded steps
    ) %>% 
    ungroup()
  
  print(paste0("n trips with artificial sections: ", length(unique(robot_trips[[as.character(id_col)]]))))
  print("IDs of trips to remove:")
  print(unique(robot_trips[[as.character(id_col)]]))
  
  # Exclude flagged trips from the main dataset
  df_clean <- df %>%
    filter(!(.data[[as.character(id_col)]] %in% robot_trips[[as.character(id_col)]]))
  print(paste0("n total trips after cleaning: ", length(unique(df_clean[[as.character(id_col)]]))))
  
  return(df_clean)
}

# Example usage: 
df_GPS_HMM_clean <- filter_robot_trips(df = df_GPS_HMM, id = "bird_tripID", threshold = 5)
