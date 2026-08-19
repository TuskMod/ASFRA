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

LoadSurveillanceDesign<-function(parameters){
  
  if(parameters$sample != 1){
    return(NULL)
  }

  # Read sampling file
  # sampling file should be in the tile of interest 
  county_shapefile <- "../Counties-of-Interest/SarasotaFL/partnership_shapefiles_24v2_12115/PVS_24_v2_county_12115.shp"
  sample.design <- read.csv("../Counties-of-Interest/SarasotaFL/sampling-Sarasota-FL.csv") # maybe make this more generic so user doesn't have to put in path?
#  county_shapefile <- "../Counties-of-Interest/MuscogeeGA/partnership_shapefiles_24v2_13215/PVS_24_v2_county_13215.shp"
#  sample.design <- read.csv("../Counties-of-Interest/MuscogeeGA/sampling-Muscogee-GA.csv") # maybe make this more generic so user doesn't have to put in path?
  
   # need to ensure surveillance county is in tile that is being looked at 
   working_crs <- 3087
 # min_sample_thresh = 10
  # Aggregate by week and store weeks to be sampled in variable
  # Convert "collection_date" to Date format
  names(sample.design)[1] <- "dates"  # Rename the first column to 'dates'
  sample.design$dates <- as.Date(sample.design$dates, format = "%m/%d/%Y") # change to standard date format
  sample.design$shp_file <- county_shapefile
  
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
  sample.design$sampling_loc <- vector("list", nrow(sample.design))
  
  return(sample.design)

}

##
# Purpose: Matching Surveillance Locations to specific grid cells

##

MatchGridstoCell <- function(sample.prep,parameters,grid,lands_data){
  if(parameters$sample != 1){
    return(NULL)
  }
  # get CRDS for tile
  curr_tile <- terra::rast(lands_data[[1]])
  custom_crs <- crs(curr_tile)
  # now - get CRS for sample.design
  # assuming it a CRS!
  sample_coords <- st_as_sf(sample.prep, coords = c("longitude", "latitude"))
  st_crs(sample_coords) <- 4269
  # pull out x and y coords from sample.prep, convert to sf object
  sample_sf_transformed <- st_transform(sample_coords,custom_crs)
  # transform coordinates to meter based CRS
  sample_coords_transformed <- st_coordinates(sample_sf_transformed) # extract the coordinates only
  #sample_coords_transformed <- sample_coords_transformed / 1000 # convert to km
#  sample_bbox <- st_bbox(st_read(sample.prep$shp_file[[1]]))
 # sample_xmin <- sample_bbox["xmin"][[1]]
  #sample_xmax <- sample_bbox["xmax"][[1]]
  #sample_ymin <- sample_bbox["ymin"][[1]]
  #sample_ymax <- sample_bbox["ymax"][[1]]

  current_extent <- ext(curr_tile)
  tile_centroid <- centroids(curr_tile)
  tile_x <- tile_centroid$x[[1]]
  tile_y <- tile_centroid$y[[1]]
  rast_xmin <- xmin(current_extent)[[1]]
  rast_xmax <- xmax(current_extent)[[1]]
  rast_ymin <- ymin(current_extent)[[1]]
  rast_ymax <- ymax(current_extent)[[1]]
  
  # Add column to sample.design so the cell # where sampling occurs can be updated
  #sample.design$sampling_loc <- 0L
  sample.prep$sampling_loc <- vector("list", nrow(sample.prep))
  # Loop over each sampling point and check proximity to grid centroids
  for (i in 1:nrow(sample_coords_transformed)) {
    # Extract the x and y coordinates of the current sample point
    sample_x <- sample_coords_transformed[i, 1]
    sample_y <- sample_coords_transformed[i, 2]
    
    
    
    # Get the acreage for the current sample point (assuming you have an 'acres' column in the dataframe)
    names(sample.prep)[5] <- "acres"  # Rename the first column to 'dates'
    acres <- sample.prep$acres[i]
    
    # Convert acres to square kilometers
    area_km2 <- acres * 0.00404686
    
    # Resolution of a grid cell in km2 calculation
    # VR: need to check this conversion more? is this correct??
    #numerator = inc * 1000 * inc * 1000
    #denominator = 1000000
    #grid_cell_area = numerator/denominator
    grid_cell_area = parameters$inc * parameters$inc
    
    # Calculate the number of grid cells to sample based on the area (rounding up to ensure entire area is covered)
    num_cells_to_sample <- ceiling(area_km2 / grid_cell_area)
    
    # Calculate the Euclidean distance (dist between 2 points) from this sample point to each centroid in the grid
    # square root [(xf-xi)^2 + (yf-yi)^2]
    # get centroid of grid! use this to find each cell's centroid lat/long 
    
    tile_coords <- crds(curr_tile)
    distances <- (sqrt((tile_coords[,1] - sample_x)^2 + (tile_coords[,2] - sample_y)^2))/1000
  
    # Check if the minimum distance is within the threshold
    # Using threshold of 100 meters
    min_sample_thresh <- parameters$inc/sqrt(2)
    #min_sample_thresh <- 100
    if (min(distances) <= min_sample_thresh) {
      # Get indices of the closest `num_cells_to_sample` grid cells
      sorted_indices <- order(distances)
      sampled_cell_indices <- sorted_indices[1:num_cells_to_sample]
      
      # Store the sampled grid cell indices
      sample.prep$sampling_loc[[i]] <- sampled_cell_indices
    } else {
      # If no nearby grid cell found, store NA or empty list
      sample.prep$sampling_loc[[i]] <- NA  # or list() if you prefer
    }
  }
  
  return(sample.prep)
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

FindSurveillanceTiles <- function(parameters,tile_path,sample.design){
  #print(parameters)
  if(parameters$sample != 1){
    return(NULL)
  }
  fs <- list.files(tile_path, full.names=TRUE)
  nm <- unlist(tstrsplit(fs, '/', keep=5))
  county_shapefile <- sample.design$shp_file[[1]]
  curr_tile <- terra::rast(fs[1])
  custom_crs <- crs(curr_tile)
  
  # need surveillance area - so there's some bounding box 
  # associated with the area , an extent 
  right_plands <- vector(mode="list")
  plands_names <- c()
  counter = 1
  
  sample_geom <-  st_read(county_shapefile)
  st_crs(sample_geom) <- 4269
  # pull out x and y coords from sample.prep, convert to sf object
  sample_sf_transformed <- st_transform(sample_geom,custom_crs)
  sample_bbox <- st_bbox(sample_sf_transformed)
  # transform coordinates to meter based CRS
  sample_xmin <- sample_bbox["xmin"][[1]]
  sample_xmax <- sample_bbox["xmax"][[1]]
  sample_ymin <- sample_bbox["ymin"][[1]]
  sample_ymax <- sample_bbox["ymax"][[1]]
  for(fi in 1:length(fs)){
    fi_tile <- terra::rast(fs[fi])
   # custom_crs <- crs(curr_tile)
    #curr_tile <- project(curr_tile,"NAD83")
    current_extent <- ext(fi_tile)
    rast_xmin <- xmin(current_extent)[[1]]
    rast_xmax <- xmax(current_extent)[[1]]
    rast_ymin <- ymin(current_extent)[[1]]
    rast_ymax <- ymax(current_extent)[[1]]
    
    ## is sample within xmin and xmax?
    # need to grab a counties x and y coordinates!!
    old_counter <- counter
    if((sample_ymin >= rast_ymin) & (sample_ymax <= rast_ymax)){ 
       if ((sample_xmin >= rast_xmin) & (sample_xmax <= rast_xmax)){
         right_plands[[counter]] <- terra::mean(fi_tile)
         plands_names[[counter]]<-fs[fi]
         if (length(nm[fi]) == 3){
            names(right_plands[[counter]]) <- nm[fi][3]
         }
         if (length(nm[fi]) == 2){
           names(right_plands[[counter]]) <- nm[fi][2]
         }
         counter <- counter + 1
       }
    }
    if(old_counter == counter){
      if(TileinCounty(sample_bbox,current_extent)){
        right_plands[[counter]] <- terra::mean(fi_tile)
        plands_names[[counter]]<-fs[fi]
        if (length(nm[fi]) == 3){
          names(right_plands[[counter]]) <- nm[fi][3]
        }
        if (length(nm[fi]) == 2){
          names(right_plands[[counter]]) <- nm[fi][2]
        }
        counter <- counter + 1
      }
    }
    rm(fi_tile)
  }
  
#  if(length(right_plands) >= 2){
#    new_tile <- JoinTogetherTiles(county_shapefile,right_plands)
#  }
    
  # stick lands into a sprc object
#  plands_sprc <- terra::sprc(right_plands)
#  names(plands_sprc) <- lapply(right_plands, names)
  
  return(plands_names)
}


JoinTogetherTiles <- function(parameters,tile_path,county.shp,sample.design){
  if(parameters$sample != 1){
    return(NULL)
  }
  pland_files <- list.files(tile_path, full.names=TRUE)
  nm <- unlist(tstrsplit(pland_files, '/', keep=5))
  curr_tile <- terra::rast(pland_files[1])
  custom_crs <- crs(curr_tile)

  county.dat <- st_read(county.shp[[1]])
  pland_dat <- c()
  county.centroid <- st_centroid(county.dat)
  merged_raster <- NULL
  
  coords <- st_transform(county.centroid, custom_crs) # transform coordinates to meter based CRS
  buffered_point <- st_bbox(st_buffer(coords, dist = 50000))
  sampled_extent <- ext(buffered_point)
  #ymax <- buffered_point$ymax
  #ymin <- buffered_point$ymin
  #xmax <- buffered_point$xmax
  #xmin <- buffered_point$xmin
  first_tile <- NULL
  all_rasters <- c()
  counter <- 0
  for(i in 1:length(pland_files)){
    pland <- terra::mean(terra::rast(pland_files[i]))
    custom_crs <- crs(pland)
    current_extent <- ext(pland)
    #pland_dat <- c(pland_dat,pland)
    if(TileinCounty(buffered_point,current_extent)){
       all_rasters <- c(all_rasters,pland)
       counter <- counter + 1
    }
  }
  
  rast_list <- sprc(all_rasters)
  merged_raster <- mosaic(rast_list,fun="min")
  
 # st_crs(county.centroid) <- "NAD83"
  #coords <- st_transform(county.centroid, custom_crs) # transform coordinates to meter based CRS
  #buffered_point <- st_bbox(st_buffer(coords, dist = 50000))
  #print(ext(buffered_point))
  #ymax <- buffered_point$ymax
  #ymin <- buffered_point$ymin
  #xmax <- buffered_point$xmax
  #xmin <- buffered_point$xmin
 
  # ok - this means that we need to add ANOTHER tile before we do any sort of cropping
  # using the buffered_points bbox - let us find any extra tiles! 
  ## do a check here seeing if buffered_points is truly within tile
  ## if not - search for next tiles that fit!!
  
  #crop_bbox <- ext(xmin,xmax,ymin,ymax)
  cropped_tile <- crop(merged_raster,sampled_extent,mask=TRUE)
  #print("Cropping tile?")
  #print(dim(crds(cropped_tile))[[1]])
  #print(cropped_tile)
  if((dim(crds(cropped_tile))[[1]] < 40000)){
   # print("returning normal tile!")
    pland_names <- FindSurveillanceTiles(parameters,tile_path,sample.design)
   # completed_tile <- terra::mean(terra::rast(pland_names[[1]]))
  #  names(completed_tile) <- 
    return(pland_names)
    #print(pland_names[[1]])
   # return(terra::mean(terra::rast(pland_names[[1]])))
  }
  
  rf <- writeRaster(cropped_tile,filename = file.path("Input","muscogee.tif"),overwrite=TRUE)
  pland_names <- c()
  pland_names <- c(pland_names,"Input/muscogee.tif")

  return(pland_names)
}
  



    
    


## Identifies if a tile is within a county limits
## Inputs: county.rast (county bounding box), tile.rast (tile bounding box)
## Returns: True/False based on county and tile boudning box

TileinCounty <- function(county.rast,tile.rast){
 
  county_xmin <- county.rast["xmin"][[1]]
  county_xmax <- county.rast["xmax"][[1]]
  county_ymin <- county.rast["ymin"][[1]]
  county_ymax <- county.rast["ymax"][[1]]

#  county_xmin <- xmin(county.rast)[[1]]
#  county_xmax <- xmax(county.rast)[[1]]
#  county_ymax <- ymax(county.rast)[[1]]
#  county_ymin <- ymin(county.rast)[[1]]
  
  
  tile_xmin <- xmin(tile.rast)[[1]]
  tile_xmax <- xmax(tile.rast)[[1]]
  tile_ymax <- ymax(tile.rast)[[1]]
  tile_ymin <- ymin(tile.rast)[[1]]
  
  within_x <- FALSE
  within_y <- FALSE
  
  
  # conditions when a county will NEVER fit in a tile 
  if((county_xmin > tile_xmax)|| (county_xmax < tile_xmin)){
    return(FALSE)
  }
  
  if((county_ymin > tile_ymax)  || (county_ymax < tile_ymin)){
    return(FALSE)
  }
  
  # these will find counties that will definitely be in the
  # county
  if((county_xmin <= tile_xmax) & (county_ymin <= county_ymax)){
    return(TRUE)
  }
  
  
  
  return(FALSE)
  
}




ReadTileFolders <- function(tile_fp){
  
  edge_path <- "Landscape_Setup/NND_Lands/all_tile_attribs_edge.csv"
  full_path <- "Landscape_Setup/NND_Lands/all_tile_attribs.csv"
  split_vals = strsplit(tile_fp,"/",fixed=TRUE)
  tile_fn <- split_vals[[1]][[5]]
  fn_values <- strsplit(tile_fn,"_",fixed=TRUE)
  tile_data <- NULL
  if(length(fn_values) == 2){
    tile_data <- read.csv(full_path)
  } else{
    tile_data <- read.csv(edge_path)
  }
  tile_data <- tile_data[as.numeric(fn_values[[1]][[2]])-1,]
  # all_cols format comes from the loaded in RDS file that already exists 
  all_cols <- c(
    "index", "nnd_med", "nnd_range", "nnd_mean", "nnd_sd",
    "x", "y", "state", "gamma.shape", "gamma.scale",
    "moranI", "gearyC", "tc", "mast", "rgd", "rds",
    "dayl", "prcp", "tmin", "tmax", "drt", "contag",
    "aggindex", "entropy", "simpindx", "nnd_cv",
    "mvmt.mean", "mvmt.cv"
  )
  
 # tile_data <- tile_data[, -1]
  
  colnames(tile_data) <- all_cols
  return(tile_data)
}
  
