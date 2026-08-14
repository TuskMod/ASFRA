## VR: Madison's Surveillance function! 
# takes in a yearly sample design, and returns 
# results of initial surveillance - and pop data

Surveillance <- function(pop, i, sample.design, parameters) {

  # Check for current week in sample.design
  # rows_to_sample <- sample.design[sample.design$fiscal_week == i, ]
  fiscal_week_current <- ((i - 1) %% 52) + 1
  # so here, rows_to_sample comes from real data
  rows_to_sample <- sample.design[sample.design$fiscal_week == fiscal_week_current, ]
  
  
  # Initialize outputs
  live_infectious_sampled <- 0
  live_recovered_sampled <-0
  dead_infected_sampled <- 0 #POSdead
  live_infected_sampled_locs <- character(0)
  dead_infected_sampled_locs <- character(0)
  pigs_sampled <- 0
  cells_sampled <- character(0)
  infectious_on_property <- 0 # currently just I
  recovered_on_property <- 0
  living_on_property <- 0
  infectiousC_on_property <- 0
  dead_on_property <- 0
  true_infected_sampled_cells <- 0 # currently just I
  true_recovered_sampled_cells <- 0
  sampled_cells_this_week <- character(0)
  total_living_pigs_sampled_cells <- 0
  true_infected_dead_sampled_cells <- 0
  total_dead_pigs_sampled_cells <- 0
  
  # If no sampling is planned for this week, return zeroes for all below variables
  if (nrow(rows_to_sample) == 0) {
    surv_summary <- data.frame(
      timestep = i,
      live_infectious_sampled = live_infectious_sampled,
      live_recovered_sampled = live_recovered_sampled,
      dead_infected_sampled = dead_infected_sampled,
      live_infected_sampled_locs = NA,
      dead_infected_sampled_locs = NA,
      pigs_sampled = pigs_sampled,
      cells_sampled = NA,
      infectious_on_property = infectious_on_property,
      recovered_on_property = recovered_on_property,
      living_on_property = living_on_property,
      infectiousC_on_property = infectiousC_on_property,
      dead_on_property = dead_on_property,
      true_infected_sampled_cells = true_infected_sampled_cells,
      true_recovered_sampled_cells = true_recovered_sampled_cells,
      total_living_pigs_sampled_cells = total_living_pigs_sampled_cells,
      true_infected_dead_sampled_cells = true_infected_dead_sampled_cells,
      total_dead_pigs_sampled_cells = total_dead_pigs_sampled_cells
    )
    return(surveillance_data = surv_summary)
  }
  
  # Loop through locations for this timestep
  for (row_index in 1:nrow(rows_to_sample)) {
    sampling_cells <- rows_to_sample$sampling_loc[[row_index]]  # cells that could be sampled
    current_quantity <- rows_to_sample$quantity[[row_index]]  # how many need to be sampled
    cells_with_pigs <- character(0)

    # Check which cells have pigs
    for (sampling_cell in sampling_cells) {
      matching_rows <- which(pop[, 3] == sampling_cell)
      if (length(matching_rows) > 0) {
        cells_with_pigs <- unique(c(cells_with_pigs, sampling_cell))
      }
    }
    rows_in_sampling_cells <- which(pop[, 3] %in% cells_with_pigs)
    # Track true number of I's in all cells that could be sampled
    infectious_on_property <- infectious_on_property + sum(pop[rows_in_sampling_cells, 10] > 0)
    # Track true number of R's in all cells that could be sampled
    recovered_on_property <- recovered_on_property + sum(pop[rows_in_sampling_cells, 11] > 0)
    # Track total living population in all cells that could be sampled
    living_on_property <- living_on_property + sum(rowSums(pop[rows_in_sampling_cells, 8:11, drop = FALSE]) > 0)
    # Track true number of C's in all cells that could be sampled
    infectiousC_on_property <- infectiousC_on_property + sum(pop[rows_in_sampling_cells, 12] > 0)
    # Track total dead population in all cells that could be sampled
    dead_on_property <- dead_on_property + sum(rowSums(pop[rows_in_sampling_cells, 12:13, drop = FALSE]) > 0)
    
    
    pigs_available_to_sample <- sum(pop[, 3] %in% cells_with_pigs)
    if (pigs_available_to_sample == 0) next
    
    ## --- Carcass Surveillance ---
    # doesn't count towards quantity needing to be sampled based on county
    carcass_rows <- which(pop[, 3] %in% cells_with_pigs)
    carcasses_C <- carcass_rows[which(pop[carcass_rows, 12] > 0)]  # infectious carcasses in area
    carcasses_Z <- carcass_rows[which(pop[carcass_rows, 13] > 0)]  # non-infectious carcasses in area
    total_carcasses <- length(carcasses_C) + length(carcasses_Z)  # total available to find
    
    n_carcasses_to_sample <- floor(parameters$CarcassDet * total_carcasses)  # round down
    all_carcasses <- c(carcasses_C, carcasses_Z)
    # Sampling carcasses if any are available
    if (n_carcasses_to_sample > 0 && length(all_carcasses) > 0) {
      sampled_carcasses <- sample(all_carcasses, min(n_carcasses_to_sample, length(all_carcasses)))
      for (row in sampled_carcasses) {
        cell <- pop[row, 3]
        
        if (pop[row, 12] > 0) {
          if (rbinom(1, 1, parameters$Sensitivity) == 1) {
            dead_infected_sampled <- dead_infected_sampled + 1
            dead_infected_sampled_locs <- unique(c(dead_infected_sampled_locs, cell))
          }
        } else if (pop[row, 13] > 0) {
          if (rbinom(1, 1, 1 - parameters$Specificity) == 1) {
            dead_infected_sampled <- dead_infected_sampled + 1
            dead_infected_sampled_locs <- unique(c(dead_infected_sampled_locs, cell))
          }
        }
        
        # Instead of removing the row here, we just track it
        cells_sampled <- unique(c(cells_sampled, cell))
      }
    }
    
    ## --- Live Pig Surveillance ---
    # counts toward quantity needing sampled
    while (pigs_sampled < current_quantity && length(cells_with_pigs) > 0) {
    #  print(cells_with_pigs)
      selected_cell <- sample(cells_with_pigs, 1, replace=FALSE)
    
      matching_rows <- which(pop[, 3] == selected_cell)
      #print("matching row found!")
      pigs_found <- FALSE
      #infec_cell, recov_cell, living_pigs, true_dead, 
      for (row in matching_rows) {
        # Infectious (I)
        if (pop[row, 10] > 0) {
          #print("I")
          pigs_sampled <- pigs_sampled + 1
          if (rbinom(1, 1, parameters$Sensitivity) == 1) {
            live_infectious_sampled <- live_infectious_sampled + 1
            live_infected_sampled_locs <- unique(c(live_infected_sampled_locs, selected_cell))
            
          }
          pigs_found <- TRUE
          sampled_cells_this_week <- unique(c(sampled_cells_this_week, selected_cell))
          pop[row,10] <- pop[row,10] - 1
          
          # Recovered (R)
        } else if (pop[row, 11] > 0) {
          #print("R")
          pigs_sampled <- pigs_sampled + 1
          if (rbinom(1, 1, parameters$Specificity) == 1) {
            live_recovered_sampled <- live_recovered_sampled + 1
            live_infected_sampled_locs <- unique(c(live_infected_sampled_locs, selected_cell))
          }
          pigs_found <- TRUE
          sampled_cells_this_week <- unique(c(sampled_cells_this_week, selected_cell))
          pop[row,11] <- pop[row,11] - 1
          # Susceptible (S)
        } else if ((pop[row, 8] > 0)) {
          #print("S")
          pigs_sampled <- pigs_sampled + 1
          if (rbinom(1, 1, 1 - parameters$Specificity) == 1) {
            false_positives <- false_positives + 1  # Optional: Track separately
            live_infected_sampled_locs <- unique(c(live_infected_sampled_locs, selected_cell))
          }
          #print("passed")
          pigs_found <- TRUE
          sampled_cells_this_week <- unique(c(sampled_cells_this_week, selected_cell))
          #print("sampled cells this week")
          #print(sampled_cells_this_week)
          pop[row,8] <- pop[row,8] - 1
          #print("S Done")
        } else if ((pop[row, 9] > 0)) {
          #print("E")
          pigs_sampled <- pigs_sampled + 1
          if (rbinom(1, 1, 1 - parameters$Specificity) == 1) {
            false_positives <- false_positives + 1  # Optional: Track separately
            live_infected_sampled_locs <- unique(c(live_infected_sampled_locs, selected_cell))
          }
          pigs_found <- TRUE
          sampled_cells_this_week <- unique(c(sampled_cells_this_week, selected_cell))
          pop[row,9] <- pop[row,9] - 1
        } 
        
        if (pigs_sampled >= current_quantity) break
      }
      
      if (!pigs_found) {
        cells_with_pigs <- setdiff(cells_with_pigs, selected_cell)
      }
    }
  }
  
  # Count how many I or R pigs are in the cells that were actually sampled
  if (length(sampled_cells_this_week) > 0) {
    # Use the sampled cells list to calculate the number of infected pigs (I + R) in the sampled cells
    rows_in_sampled_cells <- which(pop[, 3] %in% sampled_cells_this_week)
    true_infected_sampled_cells <- sum(rowSums(pop[rows_in_sampled_cells, 10, drop = FALSE]) > 0)
    true_recovered_sampled_cells <- sum(rowSums(pop[rows_in_sampled_cells, 11, drop = FALSE]) > 0)
    total_living_pigs_sampled_cells <- sum(rowSums(pop[rows_in_sampled_cells, 8:11, drop = FALSE]) > 0)
    true_infected_dead_sampled_cells <- sum(rowSums(pop[rows_in_sampled_cells, 12, drop = FALSE]) > 0)
    total_dead_pigs_sampled_cells <- sum(rowSums(pop[rows_in_sampled_cells, c(12, 13), drop = FALSE]) > 0)
   # print("passed through")
    
  } else {
    true_infected_sampled_cells <- 0
    total_living_pigs_sampled_cells <- 0
    true_infected_dead_sampled_cells <- 0
    true_recovered_sampled_cells <- 0
    total_dead_pigs_sampled_cells <- 0
  }
 # print("tidy output!")
  ## Build tidy output for this timestep
  if(length(live_infected_sampled_locs) == 0){
    live_infected_sampled_locs <- c(0)
  }
  else
  {
    live_infected_sampled_locs <- paste(live_infected_sampled_locs, collapse = ",")
  }
  if(length(dead_infected_sampled_locs) == 0){
    dead_infected_sampled_locs <- c(0)
  }
  else
  {
    dead_infected_sampled_locs <- paste(dead_infected_sampled_locs, collapse = ",")
  }
  if(length(cells_sampled) == 0){
    cells_sampled <- c(0)
  }
  else
  {
    cells_sampled <- paste(cells_sampled, collapse = ",")
  }
  surv_summary <- data.frame(
    timestep = i,
    live_infectious_sampled = live_infectious_sampled, #POSlive
    live_recovered_sampled = live_recovered_sampled,
    dead_infected_sampled = dead_infected_sampled, #POSdead
    # POSlive_locs
    live_infected_sampled_locs = live_infected_sampled_locs,
    # POSdead_locs
    dead_infected_sampled_locs = dead_infected_sampled_locs,
    pigs_sampled = pigs_sampled,
    cells_sampled = ifelse(length(cells_sampled) == 0, "none", paste(unique(cells_sampled), collapse = ",")),
    infectious_on_property = infectious_on_property,
    recovered_on_property = recovered_on_property,
    living_on_property = living_on_property,
    infectiousC_on_property = infectiousC_on_property,
    dead_on_property = dead_on_property,
    true_infected_sampled_cells = true_infected_sampled_cells,
    true_recovered_sampled_cells = true_recovered_sampled_cells,
    total_living_pigs_sampled_cells = total_living_pigs_sampled_cells,
    true_infected_dead_sampled_cells = true_infected_dead_sampled_cells,
    total_dead_pigs_sampled_cells = total_dead_pigs_sampled_cells
  )
  
  return(surv_summary)
}

