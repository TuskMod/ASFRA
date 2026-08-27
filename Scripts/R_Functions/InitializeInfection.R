InitializeInfection <- function(pop, centroids, grid, parameters){
	
######################################
######## Initialize Infection ######## 
######################################
#num_inf_0=1 #how many pigs to infect starting off
	num_inf_0 <- parameters$num_inf_0

	if(parameters$spawn_type == "random_pref"){
	  navigable_land <- which(grid[,8] > 0)
	  possible_cells <- (grid[navigable_land,1])
	  infect_cell <- sample(possible_cells,1)
	}
	
	if(parameters$spawn_type == "randomhypercube_pref"){
	  # create a lat/long matrix to do the hypercube sampling
	  latlong_matrix <- randomLHS(500000,2)
	  
	}
	
	if(parameters$spawn_type == "land_pref"){
	  #initialize needed objects
	  cells <- nrow(centroids)
	  
	  #Get the initializing number of sounders
	  N0 <- pop_init_args[1]
	  ss <- pop_init_args[2]
	  
	  sn_i <- N0/ss
	  pref.wt <- sum(grid[,8])/cells
	  #assign to cells with weighted preference according to column 8 values
	  ## this means maximum of 1 sounder per cell in starting configuration
	  infect_cell <- rbinom(cells, 1, grid[,8] * ((sn_i / cells) / pref.wt))
	}
	
	if(parameters$spawn_type != "land_pref" && parameters$spawn_type != "random_pref"){
	
  	location_keywords <- c("n", "nw", "ne", 
  	                       "e", "w", 
  	                       "se", "sw")
  	
  	first_phrase <- split(parameters$spawn_type,"_")[1]
  	correct_phrase <- ""
  	for(i in 1:length(location_keywords)){
  	  if(str_detect(first_phrase,location_keywords[[i]])){
  	    correct_phrase <- location_keywords[[i]]
  	    break
  	  }
  	}
  	print(correct_phrase)
	if(correct_phrase != ""){
	  print("normal choosing.")
	  # determine bounding box of grid using centroids from grid matrix
	  x_vals <- grid[, 6]  # centroid X
	  y_vals <- grid[, 7]  # centroid Y
	  
	  x_min <- min(x_vals)
	  x_max <- max(x_vals)
	  y_min <- min(y_vals)
	  y_max <- max(y_vals)
	  
	  x_mid <- (x_min + x_max) / 2
	  x_left_center <- x_min + (x_mid - x_min) / 2
	  x_right_center <- x_mid + (x_max - x_mid) / 2
	  
	  y_third <- (y_max - y_min) / 3
	  y_top <- y_max - y_third / 2
	  y_mid <- y_min + y_third + y_third / 2
	  y_bot <- y_min + y_third / 2
	  location_centers <- list(
	    "n" = c(x_mid, (y_min + y_max) / 2),
	    "ne" = c(x_left_center, y_top),
	    "nw" = c(x_right_center, y_top),
	    "e" = c(x_left_center, y_mid),
	    "w" = c(x_right_center, y_mid),
	    "se" = c(x_left_center, y_bot),
	    "sw" = c(x_right_center, y_bot)
	  )
	  
	  # here we assign the coordinates to infect_loc to be used later
	  target <- location_centers[[correct_phrase]]
	  
	  # calculate distance between all centroids and infection location
	  # distance = sqrt((x2 - x1)^2 + (y2-y1)^2)
	  distances <- sqrt((x_vals - target[1])^2 + (y_vals - target[2])^2)
	  
	  # Whichever centroid/cell has the smallest distance is the point of infection
	  infect_cell <- which.min(distances)
	  
	  
	  # future code needs the cell # where infection starts
	  infect_coords <- grid[infect_cell, 6:7]
	 
	  if(grid[infect_cell, 8] == 0){
	    stop("Cannot initialze infection in the given parametrized location")
	  }
	  
	}
	}
	#find the midpoint of the grid where the infected sounder will end up
#	midpoint <- c(median(centroids[, 1]), median(centroids[, 2]))
#	id <- which(centroids[, 1] >= midpoint[1] & centroids[, 2] >= midpoint[2])[1] #location on grid closest to midpoint

    # generate a sounder of one infected individual at num_inf_0 points (usually 1)
	  # manually set to be 1, since it is a sounder of size one
	  # VR: changed to heterogeneous here
	
    infected <- InitializeSounders(centroids, grid, c(infect_cell, num_inf_0), pop_init_type="init_single", pop_init_grid_opts=parameters$pop_init_grid_opts)
  
    # (manually change state values so S=0 and I = 1)
   # pop[pop[,2] == infected[,2],] <- infected 
    infected[, 8] <- 0
    infected[, 10] <- 1
    new_pop <- as.data.frame(pop) %>%
                filter(cell != infect_cell)
    new_pop <- as.matrix(new_pop)
    ## could have the option to change initial infected sounder size
    #combine infected pig with pop matrix
    pop <- rbind(new_pop, infected)
    #print(pop[pop[,3] == infect_cell,])
    #pop[pop[,3] == infect_cell,] <- infected
	return(pop)
}
