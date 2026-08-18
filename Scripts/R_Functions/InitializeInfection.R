InitializeInfection <- function(pop, centroids, grid, parameters){
	
######################################
######## Initialize Infection ######## 
######################################
#num_inf_0=1 #how many pigs to infect starting off
	num_inf_0 <- parameters$num_inf_0

	if(parameters$spawn_type == "random_pref"){
	  navigable_land <- which(centroids[,3] > 0)
	  random_row <- sample(navigable_land,1)
	}
	if(parameters$spawn_type == "lc_pref"){
	  #initialize needed objects
	  cells <- nrow(centroids)
	  
	  #Get the initializing number of sounders
	  N0 <- pop_init_args[1]
	  ss <- pop_init_args[2]
	  
	  sn_i <- N0/ss
	  pref.wt <- sum(grid[,8])/cells
	  #assign to cells with weighted preference according to column 8 values
	  ## this means maximum of 1 sounder per cell in starting configuration
	  assigns <- rbinom(cells, 1, grid[,8] * ((sn_i / cells) / pref.wt))
	}
	#find the midpoint of the grid where the infected sounder will end up
	midpoint <- c(median(centroids[, 1]), median(centroids[, 2]))
	id <- which(centroids[, 1] >= midpoint[1] & centroids[, 2] >= midpoint[2])[1] #location on grid closest to midpoint

    # generate a sounder of one infected individual at num_inf_0 points (usually 1)
	  # manually set to be 1, since it is a sounder of size one
	  # VR: changed to heterogeneous here
    infected <- InitializeSounders(centroids, grid, c(id, num_inf_0), pop_init_type="init_single", pop_init_grid_opts=parameters$pop_init_grid_opts)
  
    # (manually change state values so S=0 and I = 1)
    infected[, 8] <- 0
    infected[, 10] <- 1
    ## could have the option to change initial infected sounder size
    #combine infected pig with pop matrix
    pop <- rbind(pop, infected)

	return(pop)
}
