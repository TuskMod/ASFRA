    library(tidyverse)
    library(sf)
    library(terra)
    
    
    apply_color <- function(val){
      
      color_vals <- c("24565" = "#FF0000", 
                      "614" = as.character("#228B22"),
                      "3760" = as.character("#9370DB"),
                      "27099" = as.character("#0000FF"),
                      "5601" = as.character("#A52A2A"),
                      "503" = as.character("#FFFF00"),
                      "1382" = as.character("#8B3E2F"),
                      "500" = as.character("#FF1493"),
                      "12000" = as.character("#FF8C00"),
                      "6136" = as.character("#000000"))
      
      return (color_vals[as.character(val)][[1]][[1]])
    }
    # will need to obsfucate these names...
    county_df = read.csv("county_data.csv")
    
    county_shapefile <- county_df$county_shp[[4]]
    sample.design <- read.csv(county_df$sample_design[[4]])
    
    #  color where sounders get initialized
    #  color the common areas making up where the surveillance locations are 
    
    summarized_sample <- sample.design %>%
      group_by(acres) %>%
      summarize(total_sampled = sum(quantity))
    
    all_acres <- unique(sample.design$acres)
    
   
    latlong_props <- sample.design %>%
      count(acres) 
    
    sample_coords <- st_as_sf(sample.design, coords = c("longitude", "latitude"), crs=4269) # pull out x and y coords from sample.design, convert to sf object
    #sample_sf_transformed <- st_transform(sample_coords, custom_crs) # transform coordinates to meter based CRS
    sample_coords_transformed <- st_coordinates(sample_coords) # extract the coordinates only
    
    custom_crs <- crs(sample_coords)
    
    acre_coords <- st_as_sf(
      data.frame(x = sample_coords_transformed[,1] , y = sample_coords_transformed[,2] ),
      coords = c("x", "y"),
      crs = custom_crs
    )
    
    acre_coords$vals <- sample.design$acres
    
    print(summarized_sample)
    print(latlong_props)
    real_color_vals <- c(as.character("#FF0000"), 
                    as.character("#228B22"),
                    as.character("#6a3d9a"),
                    as.character("#0000FF"),
                    as.character("#fdbf6f"),
                    as.character("#cab2d6"),
                    as.character("#8B3E2F"),
                    as.character("#FF1493"),
                    as.character("#FF8C00"),
                    as.character("#000000"))
    county_shp <- st_read(county_shapefile)
    
    acre_coords$color_val <- lapply(acre_coords$vals,apply_color)
    
    acre_coords$vals <- as.character(acre_coords$vals)
    
    ggplot() +
      geom_sf(data = county_shp) + 
      geom_sf(data = acre_coords, aes(color=acre_coords$vals), size = 3)+ 
      scale_color_manual(name="Acres",values = real_color_vals)
    
   # ggsave("county_save.png")
    #+
    #  geom_sf_label(data= acre_coords, aes(label=vals))
