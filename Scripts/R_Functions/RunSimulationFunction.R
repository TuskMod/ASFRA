
### Initialize grid(s): ---------------
    #Method for class 'list'
        #Initialize_Grids(object, parameters=parameters, movement=c(parameters$shape,parameters$rate))
    #Arguments
                #movement
                    #default is vector with gamma distribution shape and rate fed from parameters file

      #grid.opts- 
        #"homogenous" or "heterogeneous", default for class numeric is "homogenous" and default for type SpatRaster and SpatRasterCollection is "ras"
          #ras- 
            #use input raster to create grid
          #homogeneous-
            #creates grid with 7 columns, no land class variables
          #heterogeneous-
            #creates a neutral random landscape model with X lc variables
    #Value
      #a nested list of grid parameters
RunSimulationReplicates <- function(land_grid_list, parameters, variables, cpp_functions, reps, mv.parms){


## check if there are already outputs for a given combination; these should be skipped

    dirlist <- list.dirs()
    if ('./Output/tm.mat' %in% dirlist == FALSE){dir.create('./Output/tm.mat')}
    if ('./Output/summ.vals' %in% dirlist == FALSE){dir.create('./Output/summ.vals')}
    if ('./Output/incidence' %in% dirlist == FALSE){dir.create('./Output/incidence')}
#     if ('./Output/detections' %in% dirlist == FALSE){dir.create('./Output/detections')}
#     if ('./Output/allzone' %in% dirlist == FALSE){dir.create('./Output/allzone')}
    if ('./Output/solocs.all' %in% dirlist == FALSE){dir.create('./Output/solocs.all')}

    names(variables)[names(variables) == "density"] <- "dens"

    # Filter variables (parameters with >1 value) out of parameters
    ## selecting variables is done in SetVarParms.R
    parameters <- parameters[names(parameters) %in% names(variables) == FALSE]
    parameters <- parameters[-grep('.\\_\\_.', names(parameters))]

    #Pull needed parms from parameters for all reps
    list2env(parameters, .GlobalEnv)

    setDT(mv.parms)
    # looping table for mapply
#     lvtable <- CJ(vars = seq(nrow(variables)), land = 1:4, rep = seq(reps))
#     lvtable <- CJ(vars = seq(nrow(variables)), land = seq(length(land_grid_list)), rep = seq(reps))
    lvtable <- CJ(vars = seq(nrow(variables)), land = unique(mv.parms[,index]), rep = seq(reps))

#     lvtable[,':='(tm.mat = 0, summ.vals=0, incidence=0, solocs.all=0)]

    tm.mat.in <- list.files('./Output/tm.mat')
    summ.vals.in <- list.files('./Output/summ.vals')
    incidence.in <- list.files('./Output/incidence')
    solocs.all.in <- list.files('./Output/solocs.all')
    if (out.repl != TRUE & length(c(tm.mat.in, summ.vals.in, incidence.in, solocs.all.in)) != 0){

        splt.check <- function(nlst, nm){
            instr <- as.data.table(tstrsplit(nlst, '_', keep=2:4))
            instr <- instr[,lapply(.SD, function(x) unlist(regmatches(x, gregexpr('[0-9]', x))))]
            setnames(instr, c('r','l','v'))
            instr[,tmp := 1]
            setnames(instr, 'tmp', nm)
            return(instr)
        }
        tm.tab <- unique(as.data.table(rbindlist(lapply(tm.mat.in, splt.check, nm = 'tm.mat'))))
        summ.tab <- unique(as.data.table(rbindlist(lapply(summ.vals.in, splt.check, nm = 'summ.vals'))))
        incid.tab <- unique(as.data.table(rbindlist(lapply(incidence.in, splt.check, nm = 'incidence'))))
        solocs.tab <- unique(as.data.table(rbindlist(lapply(solocs.all.in, splt.check, nm = 'solocs.all'))))
        bndtab <- merge(tm.tab, summ.tab, by=c('v','l','r'), all=TRUE)
        bndtab <- merge(bndtab, incid.tab, by=c('v','l','r'), all=TRUE)
        bndtab <- merge(bndtab, solocs.tab, by=c('v','l','r'), all=TRUE)
        bndtab[is.na(bndtab)] <- 0
        bndtab[,names(.SD) := lapply(.SD, as.numeric)]

        lvtable2 <- merge(lvtable, bndtab, by.x=c('vars', 'land', 'rep'), by.y=c('v', 'l', 'r'), all=TRUE)
        lvtable2 <- lvtable2[is.na(tm.mat) | is.na(summ.vals) | is.na(incidence) | is.na(solocs.all),]
        lvtable2 <- lvtable2[,.(vars, land, rep)]
        print(nrow(lvtable2))
    } else {
        lvtable2 <- lvtable
    }


    # movement parameters from NND landscape selection
    setDT(mv.parms)

    # loops over combinations of variables, lands, and reps
#     rep.list <- mapply(function(v.val, l.val, r.val){
    lgl.index <- unlist(lapply(land_grid_list, function(x) x[4]))

    mapply(function(v.val, l.val, r.val){

        # add in vars
        vars <- variables[v.val,]
        vars <- as.list(vars)
        list2env(vars, .GlobalEnv)

        #calc vals based on variables
        N0 <- dens*area
        K <- N0*1.5
        parameters <- c(parameters, vars)
        parameters$K <- K

        # movement parameters
        lgl.entry <- which(lgl.index == l.val)
        parameters$shape <- as.numeric(mv.parms[lgl.entry, gamma.shape])
        parameters$rate <- as.numeric(mv.parms[lgl.entry, gamma.rate])

        #loop through landscapes
        centroids <- land_grid_list[[lgl.entry]]$centroids
        grid <- land_grid_list[[lgl.entry]]$grid
        lname <- land_grid_list[[lgl.entry]]$names
        print(paste('l.val == lname:', l.val==lname))

        # create sounders in starting locations according to N0 and ss parameters
        pop <- InitializeSounders(centroids, grid, c(N0, ss), pop_init_grid_opts)
        # add an infected individual near the center of the simulation space
        pop <- InitializeInfection(pop, centroids, grid, parameters)
        # pre-create outputs to catch output data
        outputs <- Initialize_Outputs(parameters)
        # Run simulation
        out.list <- SimulateOneRun(outputs, pop, centroids, grid, parameters, cpp_functions, K, v.val, l.val, r.val)
        # Handle outputs, including writing storage files
        rep_outputs(out.list, v.val, l.val, r.val, parameters, out.opts)
#         rep.out <- rep_outputs(out.list, v.val, l.val, r.val, parameters, out.opts)
#         return(rep.out)
        gc()
        return(NULL)
    },
    v.val=lvtable2[,vars], l.val=lvtable2[,land], r.val = lvtable2[,rep])

#     tm.mat <- rbindlist(rep.list[1,])
#     summ.vals <- rbindlist(lapply(rep.list[2,], as.data.table))
#     incidence <- rbindlist(rep.list[3,which(lapply(rep.list[3,], ncol) > 1)])
#     detections <- rbindlist(lapply(rep.list[4,][!is.na(rep.list[4,])], as.data.table))
#     allzone <- rbindlist(lapply(rep.list[5,][!is.na(rep.list[5,])], as.data.table))
#     solocs.all <- rbindlist(rep.list[6,])
#
#     return(list('tm.mat' = tm.mat, 'summ.vals' = summ.vals, 'incidence' = incidence, 'detections' = detections, 'allzone' = allzone, 'solocs.all' = solocs.all))#, 'wv.speed' = wv.speed))
    gc()
    return(lvtable)
}










