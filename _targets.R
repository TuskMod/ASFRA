# _targets.R

# Helpful HPC resources:
## example for how to do targets and hpc
# https://multimeric.github.io/targets-workshop/hpc.html
## setup example and explanation of some options
# https://github.com/ropensci/targets/discussions/1468
## crew.cluster implementation with a workaround (slurm script calls tar_make())
# https://github.com/cct-datascience/targets-uahpc/blob/main/_targets.R
## overall slurm tutorial
# https://supercomputingwales.github.io/SCW-tutorial/04-running-jobs/
## crew.cluster manual
# https://cran.r-project.org/web/packages/crew.cluster/refman/crew.cluster.html#crew.cluster-package
## the updated method, including crew.cluster
# https://books.ropensci.org/targets/crew.html
## old methods, like clustermq and future
# https://books.ropensci.org/targets/hpc.html

# Targets setup --------------------
# setwd(this.path::this.dir())

#load libraries for targets script
library(targets)
# library(tarchetypes)
# library(crew)
# library(crew.cluster)
# library(clustermq)
# library(Rcpp)
# library(data.table)
library(geotargets)

# This hardcodes the absolute path in _targets.yaml, so to make this more
# portable, we rewrite it every time this pipeline is run (and we don't track
# _targets.yaml with git)
# i.e. _targets.yaml is an optional file that sets defaults for tar_make() via
# tar_config_get(), so with cluster hpc whatever computing these change the
# default values to those that would be helpful for that...
# tar_config_set(
#     store = file.path(this.path::this.dir(),("_targets")),
#     script = file.path(this.path::this.dir(),("_targets.R"))#,
# )

# Source functions in pipeline (could also be handled by tar_source())
lapply(list.files(file.path("Scripts","R_Functions"), full.names = TRUE, recursive = TRUE), source)


## if statement here doesn't play well with cluster option settings
# options(
#     clustermq.scheduler="slurm",
#     clustermq.scheduler="multicore"#,
#         clustermq.template = './asfra_tmpl.tmpl', # handles the attributes of the workers, things in defaults below are defaults (if resoures changes template files)
#         clustermq.defaults = list(
#                                     conda="asfra",
#                                     log_file = './log/o%j.out',
#                                     error='./log/e%j.err',
# #                                     memory=6*8129,
# #                                     memory=49152,
# #                                     n_jobs=1,
# #                                     cores=4,
# #                                     cores=12,
#                                     walltime=7200)
#         clustermq.worker.timeout = 7200)
#         )

#Load packages
tar_option_set(packages = c("Rcpp",
#                             "pracma",
                            "rdist",
#                             "tidyverse",
                            "RcppArmadillo",
                            "RcppParallel",
                            "stringr",
                            "dplyr",
                            "sf",
                            "raster",
                            "terra",
#                             "EnvStats",
#                             "deSolve",
#                             "colorspace",
                            "data.table",
                            "pdist"
# #                             "MuMIn",
#                             "lme4",
                            ,"glmmTMB"
#                             "glmnet",
#                             "glmmLasso"
),
               seed = 12345 ## can set the seed for reproducibility, or NA for non-reproducible totally stochastic -- see targets manual section 9.2
                ,error = 'stop' # for troubleshooting
#                ,error = 'null' # for production runs -- stops the errored cases and flags them for re-running later, but lets other things continue
#                ,deployment='worker'
               ,garbage_collection=TRUE
               ,workspace_on_error=TRUE
)

# Pipeline ---------------------------------------------------------

list(
#
#     ## Input raw data files -----
#     ### Input parameters file: -----------
    tar_target(parameters_txt, file.path("Parameters.txt"), format="file"),#, cue=tar_cue(mode='always')),
#
#     ### Input landscapes directory: -----------
# #     tar_target(lands_path, file.path("Landscape_Setup","NND_Lands","4_Output", "plandstest"), format="file"),
    tar_target(lands_path, file.path("Landscape_Setup","NND_Lands","4_Output", "sel_plands"), format="file"),
#
#     ### Read and format parameters file: -----------
    tar_target(parameters0, FormatSetParameters(parameters_txt)),
#
#     ### Pull variables out of parameters list (parameters with multiple values = variables, for flexibility)
    tar_target(variables, SetVarParms(parameters0)),
#
#     ## Input cpp scripts as files to enable tracking -----
#     # this target is a dead end, just returns "TRUE" if complete; precompiles c++ scripts for later use so it only does it once
    tar_target(cpp.compile, {
        Rcpp::sourceCpp("./Scripts/cpp_Functions/Movement_Fast_Generalized.cpp", cacheDir = './cppcache_mv', rebuild=TRUE)
        Rcpp::sourceCpp("./Scripts/cpp_Functions/Fast_FOI_Matrix.cpp", cacheDir = './cppcache_ffoi', rebuild=TRUE)
        TRUE}
    # set cue to "always" so this target always runs, in case cpp scripts were changed or cache was cleared/corrupted
#     ,cue=tar_cue(mode='always')
    ),
#
#   ## Initialize surface -----
#     ### Initialize grid(s): ---------------
#     # gives a list with (1) a list of all the landscape values (very large), (2) 'inc' resolution in km, (3) 'km_len' length of landscape side in km
    tar_target(land_grid_list, InitializeGrids(lands_path, parameters0)),
#
#     ### Get surface parameters from raster and add to parameters list: ---------------
    tar_target(parameters, GetSurfaceParms(parameters0, land_grid_list[[2]], land_grid_list[[3]])),
#
#     ### Get landscape-specific movement parameters
# #     tar_target(mv.params, if(parameters$grid.opts=='ras'){ readRDS(file.path('Landscape_Setup', 'NND_Lands', '4_Output', 'ldsel_test.rds')) } else { return(NA)}),
    tar_target(mv.params, if(parameters$grid.opts=='ras'){ readRDS(file.path('Landscape_Setup', 'NND_Lands', '4_Output', 'ldsel.rds')) } else { return(NA)})#,

    ### Build table of all desired landscape and parameter combinations
    ,tar_target(lvtable, combo.plans(parameters, variables, parameters$nrep, mv.params)),
# #
# #     ### Check for existing outputs and make a table of only what remains to be done
    tar_target(lvtable2, {lvnew <- check.existing(lvtable, parameters$out.repl)
                          if (nrow(lvnew) == 0){
                              return(data.table(vars=0, land=0, rep=0))
                          } else {
#                               return(lvnew[vars==1, ])
                              return(lvnew)
                          }}
                          , cue = tar_cue(mode='always'))
#
#     ## Run Model ---------------
#         # as currently configured, looks for existing output files and skips those that already exist
#         # (Becase seed is set, outputs should be the same for a given landscape, variable set, and replicate)
#         ## except that's not the case -- targets does weird things with seeds and they reset inside the model
#         # so using tar_cue(mode = 'always') really just forces it to check for existing previous runs if parameter out.repl = 0
#         # if out.repl = 1 (i.e. "replace outputs") it will run everything again (i.e. if something changes in the model, set out.repl to 1 to get all new results)
#         # out.repl = 1 does NOT delete existing files, so shorter or incomplete re-runs will not leave only the legitimate outputs
    ,tar_target(out.list, RunSimulationReplicates(land_grid_list = land_grid_list[[1]],
                                                parameters = parameters,
                                                variables = variables,
                                                mv.parms = mv.params,
                                                lvtable2 = lvtable2)
        , pattern = map(lvtable2) # tells targets to branch nodes by lines in lvtable2, which is the list of land, variable, and replicate combinations
        , iteration = 'vector'
        , cue = tar_cue(mode = 'always')
    )
#
#     , tar_target(centroid.matrix, extract.centroids(land_grid_list))
#
#     , tar_target(tm.out.files, {print(out.list); return(list.files('./Output/tm.mat/'))})#, cue = tar_cue(mode = 'always'))
#
#     , tar_target(result.outputs, data.looper(tm.out.files, centroid.matrix, lvtable)
#                , pattern = map(tm.out.files),
#                , iteration = 'list',
#                  cue = tar_cue(mode='always')) # faster than having it check for existing outputs (at least on scicomp)
#
#     , tar_target(rslt.file, fwrite(rbindlist(result.outputs), 'result.outputs.csv'))

    ## put together result table
#     , tar_target(rslt, make_rslt(rslt.file, variables, mv.params)) # in ridge_lasso.R
    , tar_target(rslt, make_rslt(variables, mv.params)) # in ridge_lasso.R

    ## ridge/lasso analysis (lasso, because there's a glmm function for it already)
#     , tar_target(lkf.table, make_lkf(rslt, kf.sz=1)) # for leave one out cross validation
#     , tar_target(lasslist, make.lasslist(rslt, c('inf.spd','est','sounder.weeks','max.inc','prop.infd'), lkf.table, lambdas=15))
#     , tar_target(cv.out, glmm.lasso(lasslist, rslt, lkf.table, cv=TRUE),
#                  pattern=map(lasslist),
#                  iteration='vector')
#     , tar_target(cv.lambdas, organize.lasso(cv.out),
#                  pattern=map(cv.out),
#                  iteration='list')
#     , tar_target(cv.plt, cv.curves(cv.lambdas))

    ## predictions based on cv outputs
    , tar_target(l.table, make_lkf(rslt, kf.sz=0)) # no cross validation
#     , tar_target(lasslist2, make.lasslist(rslt, c('inf.spd','est','sounder.weeks','max.inc','prop.infd'), l.table, lambdas=1, cvlbd=cv.lambdas))
#     , tar_target(mods.out, glmm.lasso(lasslist2, rslt, l.table, cv=FALSE),
#                  pattern=map(lasslist2),
#                  iteration='vector')

    ## analysis of outputs with just a messy glmm process
    , tar_target(lasslist2, c('inf.spd','est','sounder.weeks','max.inc','prop.infd'))
#     , tar_target(lasslist2, c('prop.infd'))
    , tar_target(mods.out, glmm.vanilla(lasslist2, rslt),
                 pattern=map(lasslist2),
                 iteration='list')


    ## generate out of sample predictions for the rest of the country
    , tar_target(tiledat, './Landscape_Setup/NND_Lands/4_Output/all_tile_attribs.csv', format='file')
    , tar_target(grid_path, './Landscape_Setup/Pipeline_SSF_Weekly/tile_grid.shp', format='file')
    , tar_target(tile_tab, data.table(values(vect(grid_path, crs='EPSG:5070'))[,c('ID','wb_inside','wb_edge')]))
    , tar_target(oos.tile.dat, oos.tiles(tiledat, tile_tab, variables))
    , tar_target(oos.tile.dat2, edge.correction(oos.tile.dat))
#     , tar_target(model.defs, return(mods.out[seq(2, length(mods.out), by=2)]))
#     , tar_target(model.names, return(mods.out[seq(1, length(mods.out), by=2)]))
#     , tar_target(preds.out, oos.preds(model.defs, lasslist2, oos.tile.dat),
#                  pattern=map(model.defs),#, lasslist2),
#                  iteration='list')

    , tar_target(preds.out, oos.preds.glmm(mods.out, lasslist2, oos.tile.dat2),
                 pattern=map(mods.out, lasslist2),#, lasslist2),
                 iteration='list')
#     , tar_target(preds.compiled, preds.compile(preds.out, oos.tile.dat))
    , tar_target(preds.compiled, Reduce(function(...) merge(..., all=TRUE), preds.out))


    ## Plotting function(s) for outputs.
    , tar_target(varset, seq(8))
    , tar_terra_vect(vector_lyr, make_vector_lyr(varset, preds.compiled, grid_path),
                     pattern=map(varset))
    , tar_target(contus_dir, file.path('Input','contus_shp','usmapplot_best.shp'), format='file')
#     , tar_target(map, maps.plot(vector_lyr, contus_dir, variables),
#                  pattern=map(vector_lyr))#, cue=tar_cue(mode='always'))
    , tar_target(map2, new.map(lasslist2, vector_lyr, contus_dir, variables, inside.only = FALSE), pattern=map(vector_lyr))
#     , tar_target(hmap, heatmap(model.defs, oos.tile.dat, model.names),#, preds.compiled),
#                  pattern=map(model.defs, model.names),
#                  iteration='list')#, cue=tar_cue(mode='always'))


) # end targets list

