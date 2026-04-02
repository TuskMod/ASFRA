# _targets.R

#Friday May 16, get parameters file set up to be able to loop
    #create function to make matrices with a row for each parameter setting
        #inputs: parameters, list of inputs, structured so that each unit is its own setting
            #check that list of inputs, names match parms with 'input' label.. otherwise do warning/stop
            #create matrix with all combinations of parameters
            #output variable parameter matrix

# Targets setup --------------------
setwd(this.path::this.dir())

#load libraries for targets script
library(targets)
library(tarchetypes)
library(crew)
library(crew.cluster)
library(clustermq)

# This hardcodes the absolute path in _targets.yaml, so to make this more
# portable, we rewrite it every time this pipeline is run (and we don't track
# _targets.yaml with git)
tar_config_set(
    store = file.path(this.path::this.dir(),("_targets")),
    script = file.path(this.path::this.dir(),("_targets.R"))
)

#Source functions in pipeline
lapply(list.files(file.path("Scripts","R_Functions"), full.names = TRUE, recursive = TRUE), source)

#set options
options(clustermq.scheduler="slurm", clustermq.template = 'asfra_tmpl.tmpl')

#Load packages
tar_option_set(packages = c("Rcpp",
                            "pracma",
                            "rdist",
                            "tidyverse",
                            "RcppArmadillo",
                            "RcppParallel",
                            "stringr",
                            "dplyr",
                            "sf",
                            "raster",
                            "terra",
                            "EnvStats",
                            "clustermq",
                            "deSolve",
                            "colorspace",
                            "data.table",
                            "pdist"),
               seed = 12345, ## can set the seed for reproducibility, or NA for non-reproducible totally stochastic -- see targets manual section 9.2
               controller = crew.cluster::crew_controller_slurm(workers=12,
                                                        seconds_idle=120,
                                                        options_cluster=crew_options_slurm(verbose=TRUE,
                                                                                           memory_gigabytes_per_cpu=8,
                                                                                           cpus_per_task=12,
                                                                                           n_tasks=20,
                                                                                           partition='scicomp-compute'),
#                                                                                            partition='scicomp-high-memory'),
                                                        profile='ASFRA'
                                                        ),
#                 error = 'stop') # for troubleshooting
               garbage_collection=TRUE,
               error = 'continue'
)

# tar_source()
# Pipeline ---------------------------------------------------------

list(
  
    ## Input raw data files -----
    ### Input parameters file: -----------
    tar_target(parameters_txt, file.path("Parameters.txt"), format="file"),

    ### Input landscapes directory: -----------
    tar_target(lands_path, file.path("Landscape_Setup","NND_Lands","4_Output", "sel_plands"), format="file"),

    ### Read and format parameters file: -----------
    tar_target(parameters0, FormatSetParameters(parameters_txt)),

    tar_target(variables, SetVarParms(parameters0)),

    ## Input cpp scripts as files to enable tracking -----
    tar_target(Fast_FOI_Matrix_script, file.path("Scripts","cpp_Functions","Fast_FOI_Matrix.cpp"), format="file"),
    tar_target(FindCellfromCentroid_script, file.path("Scripts","cpp_Functions","FindCellfromCentroid.cpp"), format="file"),
    tar_target(Movement_Fast_Generalized_script, file.path("Scripts","cpp_Functions","Movement_Fast_Generalized.cpp"), format="file"),
    tar_target(Movement_Fast_RSFavail_script, file.path("Scripts","cpp_Functions","Movement_Fast_RSFavail.cpp"), format="file"),
    tar_target(SpatialZones_fast_script, file.path("Scripts","cpp_Functions","SpatialZones_fast.cpp"), format="file"),
  
  ## Initialize surface -----
  ### Initialize grid(s): ---------------
    tar_target(land_grid_list, InitializeGrids(lands_path, parameters0)),

    ### Get surface parameters: ---------------
    tar_target(parameters, GetSurfaceParms(parameters0, land_grid_list[[2]], land_grid_list[[3]])),

    ## Get landscape-specific movement parameters
    tar_target(mv.params, if(parameters$grid.opts=='ras'){ readRDS('./Landscape_Setup/NND_Lands/4_Output/ldsel.rds') } else { return(NA)}),

    ## Run Model ---------------
    tar_force(out.list,
        RunSimulationReplicates(land_grid_list = land_grid_list[[1]],
                                parameters = parameters,
                                variables = variables,
                                cpp_functions = list(Fast_FOI_Matrix_script, Movement_Fast_Generalized_script),
                                reps = parameters$nrep,
                                mv.parms = mv.params
        )
        ,force=TRUE
    ),
      ## Copy paste everything in the {} including the {} to run simulations using targets outputs without running targets so you can read the error messages and outputs! :)
      ## {lapply(list.files('./Scripts/R_Functions/', full.names=TRUE), source) ;RunSimulationReplicates(tar_read(land_grid_list), tar_read(parameters), tar_read(variables), list(tar_read(Fast_FOI_Matrix_script), tar_read(Movement_Fast_Generalized_script)), tar_read(parameters)$nrep, tar_read(mv.params))}#, tar_read(burn.list)) }


    tar_target(plot_outputs, VisualOutputs(out.list, variables, land_grid_list, parameters))
      ## Copy paste everything in the {} including the {} to run simulations using targets outputs without running targets so you can read the error messages and outputs! :)
      ## {lapply(list.files('./Scripts/R_Functions/', full.names=TRUE), source) ; VisualOutputs(tar_read(out.list), tar_read(variables), tar_read(land_grid_list), tar_read(parameters)) }


) # end targets list

## run tar_visnetwork() to visualize targets dependencies

