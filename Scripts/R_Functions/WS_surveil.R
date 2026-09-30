calib_moves <- function(pop, i, sample.design, parameters){
  # initialize simulation with most parameters fixed except for
  # land tile value
  # density rate 
  # contact rate
  # value for b1,b2 
  # Check for current week in sample.design
  # Median infection should be between 100 - 800 new infections a month
  # spatial infection should be between 0.6 - 54 km a month 
  
  # Calibrating for a county of choice 
  # -> vary density and contact rate
  # -> 500 replicates
  # -> different values for b1
  parameters_txt <- file.path("Parameters.txt") #, cue=tar_cue(mode='always')),
  parameters0 <- FormatSetParameters(parameters_txt)
  variables <- SetVarParms(parameters0)
  sample.prep <- LoadSurveillanceDesign(parameters0)
  all_lands <- (file.path("Landscape_Setup", "custom_tile2","custom_tile","4_Output","land_tiles"))
  lands_data <- FindSurveillanceTiles(parameters0,all_lands,sample.prep)
  land_grid_list <- InitializeGrids(all_lands,lands_data,  parameters0)

  
  land_grid_list <- InitializeGrids(all_lands,lands_data,  parameters0)
  parameters <- GetSurfaceParms(parameters0, land_grid_list[[2]], land_grid_list[[3]])
  Rcpp::sourceCpp(file.path("Scripts", "cpp_Functions", "Movement_Fast_Generalized.cpp"), cacheDir = './cppcache_mv', rebuild=TRUE)
  Rcpp::sourceCpp(file.path("Scripts", "cpp_Functions", "Fast_FOI_Matrix.cpp"), cacheDir = './cppcache_ffoi', rebuild=TRUE)
  Rcpp::sourceCpp("./Scripts/cpp_Functions/Movement_Fast_Generalized.cpp", cacheDir = './cppcache_mv', rebuild=TRUE)
  Rcpp::sourceCpp("./Scripts/cpp_Functions/Fast_FOI_Matrix.cpp", cacheDir = './cppcache_ffoi', rebuild=TRUE)
  
  land_grid_list <- InitializeGrids(all_lands,lands_data,  parameters0)
  parameters <- GetSurfaceParms(parameters0, land_grid_list[[2]], land_grid_list[[3]])
  sample.design <- MatchGridstoCell(sample.prep,parameters,land_grid_list,lands_data)
  mv.params <- ReadTileFolders(lands_data[[1]])
  # for calibration we really care about 
  # the fact that  we have 500 replicates for 
  # different 
  lvtable <- combo.plans(parameters, variables, parameters$nrep, mv.params)
  

  
  RunSimulationReplicates(land_grid_list = land_grid_list[[1]],
                          parameters = parameters,
                          variables = variables,
                          mv.parms = mv.params,
                          sample.design = sample.design,
                          lvtable2 = lvtable2)
  

    
  
}