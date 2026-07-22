# Primary function to manage single simulation runs (setup, execution, output)

RunSimulationReplicates <- function(land_grid_list, parameters, variables, mv.parms, lvtable2){

    # movement parameters from NND landscape selection
    # convert to datatable
    setDT(mv.parms)

    # loops over combinations of variables, lands, and reps
    lgl.index <- unlist(lapply(land_grid_list, function(x) x$names))

    # funky mis-naming (fix at some point)
    names(variables)[names(variables) == "density"] <- "dens"

    # read in parameters
    list2env(parameters, .GlobalEnv)

    # Read in c++ functions (set here so that dynamic branching process on cluster would do this for each root and reference the correct cache locations)
    # these essentially just grab the compiled functions
    Rcpp::sourceCpp(file.path("Scripts", "cpp_Functions", "Movement_Fast_Generalized.cpp"), cacheDir = './cppcache_mv', rebuild=FALSE)
    Rcpp::sourceCpp(file.path("Scripts", "cpp_Functions", "Fast_FOI_Matrix.cpp"), cacheDir = './cppcache_ffoi', rebuild=FALSE)

    # lvtable2 is read in by dynamic branching of targets, so it only has one row
    v.val <- unlist(lvtable2[,vars])
    l.val <- unlist(lvtable2[,land])
    r.val <- unlist(lvtable2[,rep])
    if (v.val == 0 & l.val == 0 & r.val == 0) return(lvtable2)

    # read in vars
    vars <- variables[v.val,]
    vars <- as.list(vars)
    list2env(vars, .GlobalEnv)

    # calc vals based on variables
    # initial pigs total
    N0 <- dens*area
    # carrying capacity
    K <- N0*1.5

    # add vars to parameters
    parameters <- c(parameters, vars)
    parameters$K <- K

    # movement parameters from landscape tile
    lgl.entry <- which(lgl.index == l.val)
#     parameters$alpha <- 1/as.numeric(mv.parms[lgl.entry, sigdisp])
#     parameters$theta <- as.numeric(mv.parms[lgl.entry, disp])/parameters$alpha
    parameters$alpha <- mv.parms[lgl.entry, gamma.shape]
    parameters$theta <- mv.parms[lgl.entry, gamma.scale]

    # grab landscape data
    centroids <- land_grid_list[[lgl.entry]]$centroids
    grid <- land_grid_list[[lgl.entry]]$grid
    lname <- land_grid_list[[lgl.entry]]$names

    # create sounders in starting locations according to N0 and ss parameters
    pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)

    # add an infected individual near the center of the simulation space
    pop <- InitializeInfection(pop, centroids, grid, parameters)

    # pre-create outputs to catch output data
    outputs <- Initialize_Outputs(parameters)

    # Run simulation
    out.list <- SimulateOneRun(outputs, pop, centroids, grid, parameters, K, v.val, l.val, r.val)

    # Handle outputs, including writing storage files (returns NULL)
    summ.vals <- rep_outputs(out.list, v.val, l.val, r.val, parameters, out.opts)

    return(lvtable2)
}










