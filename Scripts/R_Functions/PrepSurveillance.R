#######################
####### Purpose #######
#######################

# Processes sampling file that is user inputted indicating
# collection date, longitude, latitude, quantity and acreage of sample effort
  # Currently, only option is to have landscape prediction for county of interest
  # i.e. grid.opt="ras"
# Creates a sample design matrix, indicating which grid cells to sample 
# during specific weeks 

#######################
######## Function #####
#######################

# Inputs:
  # sampling- file with sampling scheme
  # inc - grid cell resolution (not sure if this will be true with new MakeGrid implementation)

#Outputs:
#1-sample.design- dataframe containing sampling info to be used in Surveilance.R
  #col1-date
  #col2-lat
  #col3-long
  #col4-quantity
  #col5-acres
  #col6-FY start date
  #col7-week #


LoadSurveillanceDesign<-function(inc){

  # Read sampling file
  # sampling file should be in the tile of interest 
  county_shapefile <- st_read("../Counties-of-Interest/SarasotaFL/partnership_shapefiles_24v2_12115/PVS_24_v2_county_12115.shp")
  sample.design <- read.csv("../Counties-of-Interest/SarasotaFL/sampling-Sarasota-FL.csv") # maybe make this more generic so user doesn't have to put in path?
  # need to ensure surveillance county is in tile that is being looked at 
   working_crs <- 3087
  min_sample_thresh = 10
  # Aggregate by week and store weeks to be sampled in variable
  # Convert "collection_date" to Date format
  names(sample.design)[1] <- "dates"  # Rename the first column to 'dates'
  sample.design$dates <- as.Date(sample.design$dates, format = "%m/%d/%Y") # change to standard date format
  
  # Apply mutate on the entire dataframe, not just on a vector
  sample.design <- sample.design %>%
    mutate(# if collection date is in october or later, fiscal year starts in that year
      # can handle out of order dates from spreadsheet
      fiscal_year_start = as.Date(paste0(ifelse(month(dates) >= 10, year(dates), year(dates) - 1), "-10-01")), 
      # Which fiscal week are we in
      fiscal_week = as.integer(difftime(dates, fiscal_year_start, units = "weeks")) + 1,
      # Make week number stays between 1-52 (handle any overflow)
      fiscal_week = ifelse(fiscal_week < 1, fiscal_week + 52, fiscal_week),  # Handle weeks that are before the fiscal year start
      fiscal_week = ifelse(fiscal_week > 52, fiscal_week %% 52, fiscal_week)  # Ensure the week number stays between 1 and 52
    )
  
  # Decide which cells are sampling locations and add "1" to 8th column of grid
  names(sample.design)[2] <- "latitude"
  names(sample.design)[3] <- "longitude"
  
  sample_coords <- st_as_sf(sample.design, coords = c("longitude", "latitude"), crs=4269) # pull out x and y coords from sample.design, convert to sf object
  sample_sf_transformed <- st_transform(sample_coords, crs = 26917) # transform coordinates to meter based CRS
  sample_coords_transformed <- st_coordinates(sample_sf_transformed) # extract the coordinates only
  sample_coords_transformed <- sample_coords_transformed / 1000 # convert to km
  
  
  # Add column to sample.design so the cell # where sampling occurs can be updated
  #sample.design$sampling_loc <- 0L
  sample.design$sampling_loc <- vector("list", nrow(sample.design))
  
  return(sample.design)
  # Loop over each sampling point and check proximity to grid centroids
  for (i in 1:nrow(sample_coords_transformed)) {
    # Extract the x and y coordinates of the current sample point
    sample_x <- sample_coords_transformed[i, 1]
    sample_y <- sample_coords_transformed[i, 2]
    
    # Get the acreage for the current sample point (assuming you have an 'acres' column in the dataframe)
    names(sample.design)[5] <- "acres"  # Rename the first column to 'dates'
    acres <- sample.design$acres[i]
    
    # Convert acres to square kilometers
    area_km2 <- acres * 0.00404686
    
    # Resolution of a grid cell in km2 calculation
    # VR: need to check this conversion more? is this correct??
    #numerator = inc * 1000 * inc * 1000
    #denominator = 1000000
    #grid_cell_area = numerator/denominator
    grid_cell_area = inc * inc
    
    # Calculate the number of grid cells to sample based on the area (rounding up to ensure entire area is covered)
    num_cells_to_sample <- ceiling(area_km2 / grid_cell_area)
    
    # Calculate the Euclidean distance (dist between 2 points) from this sample point to each centroid in the grid
    # square root [(xf-xi)^2 + (yf-yi)^2]
    distances <- sqrt((grid[, 6] - sample_x)^2 + (grid[, 7] - sample_y)^2)
    
    # Check if the minimum distance is within the threshold
    # Using threshold of 10 meters
    if (min(distances) <= min_sample_thresh) {
      # Get indices of the closest `num_cells_to_sample` grid cells
      sorted_indices <- order(distances)
      sampled_cell_indices <- sorted_indices[1:num_cells_to_sample]
      
      # Store the sampled grid cell indices
      sample.design$sampling_loc[[i]] <- sampled_cell_indices
    } else {
      # If no nearby grid cell found, store NA or empty list
      sample.design$sampling_loc[[i]] <- NA  # or list() if you prefer
    }
  }
      
  
}

MatchGridtoCell <- function(sample.design,inc,grid){
  sample_coords <- st_as_sf(sample.design, coords = c("longitude", "latitude"), crs=4269) # pull out x and y coords from sample.design, convert to sf object
  sample_sf_transformed <- st_transform(sample_coords, crs = 26917) # transform coordinates to meter based CRS
  sample_coords_transformed <- st_coordinates(sample_sf_transformed) # extract the coordinates only
  sample_coords_transformed <- sample_coords_transformed / 1000 # convert to km
  
  
  # Add column to sample.design so the cell # where sampling occurs can be updated
  #sample.design$sampling_loc <- 0L
  sample.design$sampling_loc <- vector("list", nrow(sample.design))
  
  # Loop over each sampling point and check proximity to grid centroids
  for (i in 1:nrow(sample_coords_transformed)) {
    # Extract the x and y coordinates of the current sample point
    sample_x <- sample_coords_transformed[i, 1]
    sample_y <- sample_coords_transformed[i, 2]
    
    # Get the acreage for the current sample point (assuming you have an 'acres' column in the dataframe)
    names(sample.design)[5] <- "acres"  # Rename the first column to 'dates'
    acres <- sample.design$acres[i]
    
    # Convert acres to square kilometers
    area_km2 <- acres * 0.00404686
    
    # Resolution of a grid cell in km2 calculation
    # VR: need to check this conversion more? is this correct??
    #numerator = inc * 1000 * inc * 1000
    #denominator = 1000000
    #grid_cell_area = numerator/denominator
    grid_cell_area = inc * inc
    
    # Calculate the number of grid cells to sample based on the area (rounding up to ensure entire area is covered)
    num_cells_to_sample <- ceiling(area_km2 / grid_cell_area)
    
    # Calculate the Euclidean distance (dist between 2 points) from this sample point to each centroid in the grid
    # square root [(xf-xi)^2 + (yf-yi)^2]
    distances <- sqrt((grid[, 6] - sample_x)^2 + (grid[, 7] - sample_y)^2)
    
    # Check if the minimum distance is within the threshold
    # Using threshold of 10 meters
    if (min(distances) <= min_sample_thresh) {
      # Get indices of the closest `num_cells_to_sample` grid cells
      sorted_indices <- order(distances)
      sampled_cell_indices <- sorted_indices[1:num_cells_to_sample]
      
      # Store the sampled grid cell indices
      sample.design$sampling_loc[[i]] <- sampled_cell_indices
    } else {
      # If no nearby grid cell found, store NA or empty list
      sample.design$sampling_loc[[i]] <- NA  # or list() if you prefer
    }
  }
  
  return(sample.design)
}

#######################
####### Purpose #######
#######################

# Identifies which land raster tiles correspond to the 
# location data in the surveillance sample.design

#######################
######## Function #####
#######################

# Inputs:
# tile_path: file location to all landscape tiles
# sample.design: surveillance dataframe
# inc - grid cell resolution (not sure if this will be true with new MakeGrid implementation)

#Outputs:
# list containing filepaths to load in

FindSurveillanceTiles <- function(tile_path,sample.design){
  fs <- list.files(tile_path, full.names=TRUE)
  nm <- unlist(tstrsplit(fs, '/', keep=5))
  nm <- unlist(tstrsplit(nm, '[_.]', keep=2))
  plands_list <- vector(mode="list", length=length(fs))
  plands_names <- vector(mode="list")
  samp_xmin = sample.design
  for(fi in 1:length(fs)){
    curr_tile <- terra::rast(fs[fi])
    current_extent <- ext(curr_tile)
    rast_xmin <- xmin(current_extent)
    rast_xmax <- xmax(current_extent)
    rast_ymin <- ymin(current_extent)
    rast_ymax <- ymax(current_extent)
    
    ## is sample within xmin and xmax?
    # need to grab a counties x and y coordinates!!
    if((r_xmin > rast_xmin) & (samp_ymin > rast_ymin)){
      right_plands.append(fs[[fi]])
    }
    
    plands_names.append(nm[[fi]])
  }
  # stick lands into a sprc object
  plands_sprc <- terra::sprc(plands_list)
  names(plands_sprc) <- lapply(plands_list, names)
  
  return(plands_names)
}
  
  
