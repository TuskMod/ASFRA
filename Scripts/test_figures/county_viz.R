  library(tidyverse)
  library(terra)
  library(sf)
  
  # create grid - same as in Make_Grid.R
  # create SpatRaster
  # Let us also visualize landscape
  
  FindAllCountyCSV <- function(county_path){
      counties <- list.dirs(county_path,recursive=FALSE)
      
      samp_data <- c()
      county_shp <- c()
   
      for (i in 1:length(counties)){
        county <- counties[[i]]
        csv_file <- list.files(county,pattern="\\.csv$")
        csv_path <- paste0(county,"/",csv_file)
        #csv_path <- paste0(csv_path,csv_file)
        shp_dir <- list.dirs(county,recursive=FALSE)
        shp_files <- list.files(shp_dir,pattern="\\.shp")
        shp_file <- ""
        for (s in 1:length(shp_files)){
          fn_split <- strsplit(shp_files[s],"_")
          
          if ("county" %in% fn_split[[1]][4]){
            shp_file <- shp_files[s]
            break
          }
        }
        full_shp_path <- paste0(shp_dir,"/",shp_file)
       # full_shp_path <- paste0(full_shp_path,shp_file)
        
        county_shp <- c(county_shp,full_shp_path)
        samp_data <- c(samp_data,csv_path)
      }
      
      c_df <- data.frame(county_shp = county_shp,samp_data = samp_data)
      
      return(c_df)
  }
  
  
  
  
  county_path <- "../Counties-of-Interest/"
  all_paths <- FindAllCountyCSV(county_path)
  for(a in 1:length(all_paths$county_shp)){
    samp_file <- all_paths[a,]$samp_data[[1]]
    sample.design <- read.csv(samp_file)
    county_shapefile <- all_paths[a,]$county_shp[[1]]
    sample.design$shp_file <- county_shapefile
    
    all_lands <- (file.path("Landscape_Setup","Pipeline_SSF_Weekly","4_Output", "all_plands"))
    names(sample.design)[1] <- "dates"  # Rename the first column to 'dates'
    sample.design$dates <- as.Date(sample.design$dates, format = "%m/%d/%Y") # change to standard date format
    
   # tile_fn <- FindSurveillanceTiles(all_lands,sample.design)
  
    county_name <- strsplit(county_shapefile,"/")[[1]][3]
    save_name <- paste(county_name,".png")
    sample <- 1
    print("control loop")
    print(all_lands)
    ras <- JoinTogetherTiles(all_lands,county_shapefile,sample.design)
    custom_crs <- crs(ras)
    print("dimensions of .tif file")
    print(dim(crds(ras)))
    
    #custom_crs <- crs("epsg:4269")
    #ras <- terra::project(ras,"epsg:4269")
    len <- dim(ras)[1]
    inc <- res(ras)[1]/1000
    grid.opt <- "heterogeneous"
    
    if(dim(ras)[1]!=dim(ras)[2]){stop("raster needs to be square")}
    if(dim(ras)[1]!=len){
      stop("dimensions of raster do not match input len")
    }
    
    # Transform county shapefile to declared CRS which is county dependent
    # resolution in meters
    county_data <- st_read(county_shapefile)
    county_data <- st_as_sf(sample.design, coords = c("longitude","latitude"))
    st_crs(county_data) <- "epsg:4269"
    county_data <- st_transform(county_data,custom_crs)

    # --- Define fixed grid parameters ---
    grid_size_km <- 100
    grid_size_m <- grid_size_km * 1000  # Total grid = 82,000 meters
    grid_cell_res_m <- inc * 1000
    
    inc <- 0.5
    len <- 200
    
    #get number of cells in grid
    cells <- len^2
    
    #if grid homogeneous-- will later enter ability to alter LULC
    if("homogeneous" %in% grid.opt & sample != 1){
      #initialize empty grid matrix
      grid <- matrix(nrow=round(cells), ncol=7)
    } else {
      #initialize empty grid matrix
      grid <- matrix(nrow=round(cells), ncol=8)
    }
    
    #first column is just cell indices
    grid[, 1] <- 1:cells
    
    #Top left X coordinate of each cell
    grid[, 2] <- rep(seq(0, ((inc*len)-inc), inc), times=len)
    
    #Top left Y coordinate of each cell
    grid[, 3] <- rep(seq(0, ((inc*len)-inc), inc), each=len)
    
    #Top right X coordinate of each cell
    grid[, 4] <- rep(seq(inc, (inc*len), inc), times=len)
    
    #Top right Y coordinate of each cell
    grid[, 5] <- rep(seq(inc, (inc*len), inc), each=len)
    
    #Center X coordinate of each cell
    grid[, 6] <- rep(seq(((0+inc)/2), (((inc*len)-inc)+(inc*len))/2, inc), times=len)
    
    #Center Y coordinate of each cell
    grid[, 7] <- rep(seq(((0+inc)/2), (((inc*len)-inc)+(inc*len))/2, inc), each=len)
    
    #get centroids-only object
    centroids <- grid[, c(6, 7)]
    grid[, 8] <- round(values(ras), 2)
    #assign to centroids
    centroids <- cbind(centroids, grid[, 8])
    name <- names(ras)
    grid.list <- list("cells"=cells, "grid"=grid, "centroids"=centroids, "names"=name)
    
    grid_centroids_sf <- st_as_sf(
      data.frame(
        id = 1:nrow(grid),
        x = grid[, 6] ,
        y = grid[, 7],
        sampled = grid[, 8]  # sampling flag 0 or 1
      ),
      coords = c("x", "y"),
      crs = st_crs(custom_crs)  # use the same CRS as county_proj
    )
    
    # grab sample design now! 
    
    names(sample.design)[1] <- "dates"  # Rename the first column to 'dates'
    names(sample.design)[2] <- "latitude"
    names(sample.design)[3] <- "longitude"
    
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
    
    sample_coords <- st_as_sf(sample.design, coords = c("longitude", "latitude"), crs=4269) # pull out x and y coords from sample.design, convert to sf object
    sample_sf_transformed <- st_transform(sample_coords, custom_crs) # transform coordinates to meter based CRS
    sample_coords_transformed <- st_coordinates(sample_sf_transformed) # extract the coordinates only
    # Convert to kilometers
    #sample_coords_transformed <- sample_coords_transformed
    
    county_sf <- st_read(county_shapefile)
    st_crs(county_sf) <- 4269
    county_transformed <- st_transform(county_sf,custom_crs)
  
    
    sample_points_sf <- st_as_sf(
      data.frame(x = sample_coords_transformed[,1] , y = sample_coords_transformed[,2] ),
      coords = c("x", "y"),
      crs = st_crs(custom_crs)
    )
    
    num_samples <- sum(sample.design$quantity)
    print("tile values!")
    tile_centroids <- crds(ras)
    fig_title <- paste(county_name," Surveillance Map, Samples: ",num_samples)

   
    ggplot() +
     # geom_sf(data = grid_centroids_sf, aes(color = factor(sampled)), size = 1, alpha = 0.7) +
      geom_tile(aes(x=tile_centroids[,1],y=tile_centroids[,2],fill=values(ras))) +
    #  scale_fill_gradient(low="#e7e1ef",high="#dd1c77") + 
      scale_fill_viridis_c()+
      guides(fill=guide_colourbar(barwidth=0.5,barheight=20),title="Land Preference") +
      geom_sf(data = county_transformed,border.color="black", size = 0.9,alpha=0.4) +
      geom_sf(data = sample_points_sf, color = "magenta", size = 1, shape = 19,size=0.1,alpha=0.4) +
      geom_sf_label(border.colour="black")+
      coord_sf(crs=custom_crs)+
      labs(title=fig_title)+
      xlab("Latitude") +
      ylab("Longitude") + 
      
      theme_minimal()
    
    ggsave(save_name)
  }
  
