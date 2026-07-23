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
setwd(this.path::this.dir())

#load libraries for targets script
library(targets)
library(tarchetypes)
# library(crew)
# library(crew.cluster)
library(clustermq)
library(Rcpp)
library(data.table)

# This hardcodes the absolute path in _targets.yaml, so to make this more
# portable, we rewrite it every time this pipeline is run (and we don't track
# _targets.yaml with git)
# i.e. _targets.yaml is an optional file that sets defaults for tar_make() via
# tar_config_get(), so with cluster hpc whatever computing these change the
# default values to those that would be helpful for that...
tar_config_set(
    store = file.path(this.path::this.dir(),("_targets")),
    script = file.path(this.path::this.dir(),("_targets.R"))#,
)

# Source functions in pipeline (could also be handled by tar_source())
lapply(list.files(file.path("Scripts","R_Functions"), full.names = TRUE, recursive = TRUE), source)


# set options if using clustermq (slightly more traditional hpc approach with config file and whatnot)
# if on scicomp
# if (Sys.info()[4] == 'aapksmanlogin1'){
## if statement here doesn't play well with cluster option settings
options(
#     clustermq.scheduler="slurm",
#         clustermq.template = './asfra_tmpl.tmpl', # handles the attributes of the workers, things in defaults below are defaults (if resoures changes template files)
        clustermq.defaults = list(#conda="asfra_run",
                                    log_file = './log/o%j.out',
                                    error='./log/e%j.err',
                                    memory=6*8129,
#                                     memory=49152,
                                    n_jobs=1,
                                    cores=4,
#                                     cores=12,
                                    walltime=7200),
        clustermq.worker.timeout = 7200)

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
                            "deSolve",
                            "colorspace",
                            "data.table",
                            "pdist",
                            "MuMIn",
                            "lme4",
                            "glmmTMB",
                            "glmnet"),
               seed = 12345 ## can set the seed for reproducibility, or NA for non-reproducible totally stochastic -- see targets manual section 9.2
                ,error = 'stop' # for troubleshooting
#                ,error = 'null' # for production runs -- stops the errored cases and flags them for re-running later, but lets other things continue
               ,deployment='worker'
               ,garbage_collection=TRUE
               ,workspace_on_error=TRUE
)

# Pipeline ---------------------------------------------------------

list(
#
#     ## Input raw data files -----
#     ### Input parameters file: -----------
    tar_target(parameters_txt, file.path("Parameters.txt"), format="file"),#, cue=tar_cue(mode='always')),
  

    #     ### Read and format parameters file: -----------
    tar_target(parameters0, FormatSetParameters(parameters_txt)),
    #
    #     ### Pull variables out of parameters list (parameters with multiple values = variables, for flexibility)
    tar_target(variables, SetVarParms(parameters0)),
  
    tar_target(sample.prep,PrepSurveillance(variables$inc)),
    tar_target(lands_path, FindSurveillanceTiles((file.path("Landscape_Setup","NND_Lands","4_Output", "sel_plands")),sample.prep)),
#
#
#     ### Input landscapes directory: -----------
#+    tar_target(lands_path, file.path("Landscape_Setup","NND_Lands","4_Output", "sel_plands"), format="file"),
#

#     ## Input cpp scripts as files to enable tracking -----
#     # this target is a dead end, just returns "TRUE" if complete; precompiles c++ scripts for later use so it only does it once
     tar_target(cpp.compile, {
         Rcpp::sourceCpp(file.path("Scripts", "cpp_Functions", "Movement_Fast_Generalized.cpp"), cacheDir = './cppcache_mv', rebuild=TRUE)
         Rcpp::sourceCpp(file.path("Scripts", "cpp_Functions", "Fast_FOI_Matrix.cpp"), cacheDir = './cppcache_ffoi', rebuild=TRUE)
         Rcpp::sourceCpp("./Scripts/cpp_Functions/Movement_Fast_Generalized.cpp", cacheDir = './cppcache_mv', rebuild=TRUE)
         Rcpp::sourceCpp("./Scripts/cpp_Functions/Fast_FOI_Matrix.cpp", cacheDir = './cppcache_ffoi', rebuild=TRUE)
         TRUE}
#     # set cue to "always" so this target always runs, in case cpp scripts were changed or cache was cleared/corrupted
     ,cue=tar_cue(mode='always')
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

#   Now, add proper tiles to sample.design
    tar_target(sample.design,FindSurveillanceTiles(lands_path,sample.prep)),

#     ### Get landscape-specific movement parameters
    tar_target(mv.params, if(parameters$grid.opts=='ras'){ readRDS(file.path('Landscape_Setup', 'NND_Lands', '4_Output', 'ldsel.rds')) } else { return(NA)}),
#
#    ### Build table of all desired landscape and parameter combinations
     ### Each row is an input to the model containing replicate ID, and different parameter values
     tar_target(lvtable, combo.plans(parameters, variables, parameters$nrep, mv.params)),
#
#     ### Check for existing outputs and make a table of only what remains to be done
     tar_target(lvtable2, {lvnew <- check.existing(lvtable, parameters$out.repl)
                           if (nrow(lvnew) == 0){
                               return(data.table(vars=0, land=0, rep=0))
                           } else {
                               return(lvnew)
                           }}
                           , cue = tar_cue(mode='always')),

#     ## Run Model ---------------
#         # as currently configured, looks for existing output files and skips those that already exist
#         # (Becase seed is set, outputs should be the same for a given landscape, variable set, and replicate)
#         ## except that's not the case -- targets does weird things with seeds and they reset inside the model
#         # so using tar_cue(mode = 'always') really just forces it to check for existing previous runs if parameter out.repl = 0
#         # if out.repl = 1 (i.e. "replace outputs") it will run everything again (i.e. if something changes in the model, set out.repl to 1 to get all new results)
#         # out.repl = 1 does NOT delete existing files, so shorter or incomplete re-runs will not leave only the legitimate outputs
     tar_target(out.list, RunSimulationReplicates(land_grid_list = land_grid_list[[1]],
                                                 parameters = parameters,
                                                 variables = variables,
                                                 mv.parms = mv.params,
                                                 lvtable2 = lvtable2)
         , pattern = map(lvtable2) # tells targets to branch nodes by lines in lvtable2, which is the list of land, variable, and replicate combinations
         , iteration = 'vector'
         , cue = tar_cue(mode = 'always')
     )
#
     , tar_target(centroid.matrix, extract.centroids(land_grid_list))
#
     , tar_target(tm.out.files, {print(out.list); return(list.files('./Output/tm.mat/'))})#, cue = tar_cue(mode = 'always'))
#
     , tar_target(result.outputs, data.looper(tm.out.files, centroid.matrix)
                , pattern = map(tm.out.files),
                , iteration = 'list')
#
#

    ## put together result table
    , tar_target(rslt, make_rslt(result.outputs, variables, mv.params))
    ## get results where establishment took place
    , tar_target(rslt1, rslt[est == 1,])
    ## get survival analysis table set up
    , tar_target(tmtab3, tm.tab.3(rslt1, variables, mv.params))

#     ## establishemnt response
#     # create list of all possible models
#     , tar_target(est.dredge.list, est.glmer.dredge(rslt))
# #     evaluate each model to get quality of fit
#     , tar_target(est.dredge, model.eval(est.dredge.list, rslt),
#                  pattern=map(est.dredge.list),
#                  iteration='list')
# #     make a table of model outputs
#     , tar_target(est.modlist, ref.to.list(est.dredge, is.esc=0),
#                  pattern=map(est.dredge),
#                  iteration='list')#, cue=tar_cue(mode='never'))
#     # find models within 2 AIC points of the best fit
#     , tar_target(est.table, list.to.table(est.modlist))
#     , tar_target(est.table2, table.filter(est.table, 2, 2))
#     # get only the models within the threshold of performance
#     , tar_target(est.tablesub, est.dredge.list[est.table2[,id]])
#     , tar_target(est.modelsub, model.eval(est.tablesub, rslt),
#                  pattern=map(est.tablesub),
#                  iteration='list')
#     # average the models with delta AIC less than 2
#     , tar_target(est.avg, model.avgs(est.table2, est.modelsub))
#
#     ## infection wave speed response
#     # create list of all possible models
#     , tar_target(inf.spd.dredge.list, inf.spd.glmer.dredge(rslt1))
#     # evaluate each model to get quality of fit
#     , tar_target(inf.spd.dredge, model.eval1(inf.spd.dredge.list, rslt1),
#                  pattern=map(inf.spd.dredge.list),
#                  iteration='list')
#     # make a table of model outputs
#     , tar_target(inf.spd.modlist, ref.to.list(inf.spd.dredge, is.esc=0),
#                  pattern=map(inf.spd.dredge),
#                  iteration='list')
#     # find models within 2 AIC points of the best fit
#     , tar_target(inf.spd.table, list.to.table(inf.spd.modlist))
#     , tar_target(inf.spd.table2, table.filter(inf.spd.table, 2, 2))
#     # get only the models within the threshold of performance
#     , tar_target(inf.spd.tablesub, inf.spd.dredge.list[inf.spd.table2[,id]])
#     , tar_target(inf.spd.modelsub, model.eval1(inf.spd.tablesub, rslt1),
#                  pattern=map(inf.spd.tablesub),
#                  iteration='list')
#     # average the models with delta AIC less than 2
#     , tar_target(inf.spd.avg, model.avgs(inf.spd.table2, inf.spd.modelsub))
#
#     ## sounder weeks response
#     # create list of all possible models
#     , tar_target(sounder.weeks.dredge.list, sounder.weeks.glmer.dredge(rslt1))
#     # evaluate each model to get quality of fit
#     , tar_target(sounder.weeks.dredge, model.eval1(sounder.weeks.dredge.list, rslt1),
#                  pattern=map(sounder.weeks.dredge.list),
#                  iteration='list')
#     # make a table of model outputs
#     , tar_target(sounder.weeks.modlist, ref.to.list(sounder.weeks.dredge, is.esc=0),
#                  pattern=map(sounder.weeks.dredge),
#                  iteration='list')
#     # find models within 2 AIC points of the best fit
#     , tar_target(sounder.weeks.table, list.to.table(sounder.weeks.modlist))
#     , tar_target(sounder.weeks.table2, table.filter(sounder.weeks.table, 2, 2))
#     # get only the models within the threshold of performance
#     , tar_target(sounder.weeks.tablesub, sounder.weeks.dredge.list[sounder.weeks.table2[,id]])
#     , tar_target(sounder.weeks.modelsub, model.eval1(sounder.weeks.tablesub, rslt1),
#                  pattern=map(sounder.weeks.tablesub),
#                  iteration='list')
#     # average the models with delta AIC less than 2
#     , tar_target(sounder.weeks.avg, model.avgs(sounder.weeks.table2, sounder.weeks.modelsub))
#
#     ## proportion infected response
#     # create list of all possible models
#     , tar_target(prop.infd.dredge.list, prop.infd.beta.dredge(rslt1))
#     # evaluate each model to get quality of fit
#     , tar_target(prop.infd.dredge, model.eval1(prop.infd.dredge.list, rslt1),
#                  pattern=map(prop.infd.dredge.list),
#                  iteration='list')
#     # make a table of model outputs
#     , tar_target(prop.infd.modlist, ref.to.list(prop.infd.dredge, is.esc=0),
#                  pattern=map(prop.infd.dredge),
#                  iteration='list')
#     # find models within 2 AIC points of the best fit
#     , tar_target(prop.infd.table, list.to.table(prop.infd.modlist))
#     , tar_target(prop.infd.table2, table.filter(prop.infd.table, 2, 2))
#     # get only the models within the threshold of performance
#     , tar_target(prop.infd.tablesub, prop.infd.dredge.list[prop.infd.table2[,id]])
#     , tar_target(prop.infd.modelsub, model.eval1(prop.infd.tablesub, rslt1),
#                  pattern=map(prop.infd.tablesub),
#                  iteration='list')
#     # average the models with delta AIC less than 2
#     , tar_target(prop.infd.avg, model.avgs(prop.infd.table2, prop.infd.modelsub))
#
#     ## maximum incidence response
#     # create list of all possible models
#     , tar_target(max.inc.dredge.list, max.inc.glmer.dredge(rslt1))
#     # evaluate each model to get quality of fit
#     , tar_target(max.inc.dredge, model.eval1(max.inc.dredge.list, rslt1),
#                  pattern=map(max.inc.dredge.list),
#                  iteration='list')
#     # make a table of model outputs
#     , tar_target(max.inc.modlist, ref.to.list(max.inc.dredge, is.esc=0),
#                  pattern=map(max.inc.dredge),
#                  iteration='list')
#     # find models within 2 AIC points of the best fit
#     , tar_target(max.inc.table, list.to.table(max.inc.modlist))
#     , tar_target(max.inc.table2, table.filter(max.inc.table, 2, 2))
#     # get only the models within the threshold of performance
#     , tar_target(max.inc.tablesub, max.inc.dredge.list[max.inc.table2[,id]])
#     , tar_target(max.inc.modelsub, model.eval1(max.inc.tablesub, rslt1),
#                  pattern=map(max.inc.tablesub),
#                  iteration='list')
#     # average the models with delta AIC less than 2
#     , tar_target(max.inc.avg, model.avgs(max.inc.table2, max.inc.modelsub))
#
#     ## time to escape response
#     # create list of all possible models
#     , tar_target(tm.esc.dredge.list, tm.esc.binom.dredge(tmtab3, variables, mv.params))
# #     evaluate each model to get quality of fit
#     , tar_target(tm.esc.dredge, model.eval.esc(tm.esc.dredge.list, tmtab3),
#                  pattern=map(tm.esc.dredge.list),
#                  iteration='list')
#     # make a table of model outputs
#     , tar_target(tm.esc.modlist, ref.to.list(tm.esc.dredge, is.esc=0),
#                  pattern=map(tm.esc.dredge),
#                  iteration='list')
#     # find models within 2 AIC points of the best fit
#     , tar_target(tm.esc.table, list.to.table(tm.esc.modlist))
#     , tar_target(tm.esc.table2, table.filter(tm.esc.table, 2, 2))
#     # get only the models within the threshold of performance
#     , tar_target(tm.esc.tablesub, tm.esc.dredge.list[tm.esc.table2[,id]])
#     , tar_target(tm.esc.modelsub, model.eval.esc(tm.esc.tablesub, tmtab3),
#                  pattern=map(tm.esc.tablesub),
#                  iteration='list')
#     # average the models with delta AIC less than 2
#     , tar_target(tm.esc.avg, model.avgs(tm.esc.table2, tm.esc.modelsub))
#
#     ## aggregate into lolo predictions
#     # combine best-performing dredge models and average models
#     , tar_target(all.mods, combine.models( est.avg ,
#                                             inf.spd.avg ,
#                                             sounder.weeks.avg ,
#                                             prop.infd.avg ,
#                                             max.inc.avg ,
#                                             tm.esc.avg
#     ))
#     # set up combinations of models and cv land exclusions
#     , tar_target(lolo.combos, lolo.table(mv.params, all.mods, rslt))
#
#     # get yhat values for each model within 2 AIC of the best one, for each leave one landscape out cross validation
#     , tar_target(lolo.test.fit, lolo.predict(lolo.combos, rslt, rslt1, tmtab3),
#                  pattern=map(lolo.combos),
#                  iteration='vector')

    # calculate rmse for each candidate prediction model for each response variable
#     , tar_target(lolo.rmse, calc.rmse(lolo.test.fit, lolo.combos))
    # generate predictive models for the best rmse values by response variable
#     , tar_target(predmods.list, get.pred.mods(lolo.rmse, all.mods, rslt, rslt1, tmtab3))

    # generate out of sample predictions for the rest of the country
#    , tar_target(tiledat, './Landscape_Setup/NND_Lands/all_tile_attribs.csv', format='file')
#    , tar_target(tiledat_edge, './Landscape_Setup/NND_Lands/all_tile_attribs_edge.csv', format='file')
#    , tar_target(oos.tile.dat, oos.tiles(tiledat, rslt, variables))
#    , tar_target(oos.tile.dat.edge, oos.tiles(tiledat_edge, rslt, variables))
#    , tar_target(model.defs, return(ridge_lasso[[2]]))
#    , tar_target(model.coefs, return(t(ridge_lasso[[1]])))
#    , tar_target(lambda.mins, as.numeric(model.coefs[1:6,28]))
#    , tar_target(preds.out, oos.preds(model.defs, lambda.mins, oos.tile.dat),
#                 pattern=map(model.defs, lambda.mins),#, model.coefs[1:6,28]),
#                 iteration='list')
#    , tar_target(preds.out.edge, oos.preds(model.defs, lambda.mins, oos.tile.dat.edge),
#                 pattern=map(model.defs, lambda.mins),
#                 iteration='list')
#    , tar_target(preds.compiled, preds.compile(preds.out, preds.out.edge, oos.tile.dat, oos.tile.dat.edge))
#
#     ## Plotting function(s) for outputs.
#     , tar_target(map, maps.plot(preds.compiled, rlst, rslt1, tmtab3, variables))
#    , tar_target(hmap, heatmap(model.defs, oos.tile.dat, oos.tile.dat.edge, lambda.mins),#, preds.compiled),
#                 pattern=map(model.defs, lambda.mins),
#                 iteration='list')
    ## move some of these further up (landscape selection map could go to NND_Lands pipeline)
#     ,tar_target(plot_outputs, VisualOutputs(out.list, variables, land_grid_list, parameters))
      ## Copy paste everything in the {} including the {} to run simulations using targets outputs without running targets so you can read the error messages and outputs! :)
      ## {lapply(list.files('./Scripts/R_Functions/', full.names=TRUE), source) ; VisualOutputs(tar_read(lvtable), tar_read(variables), tar_read(land_grid_list), tar_read(parameters)) }


) # end targets list


