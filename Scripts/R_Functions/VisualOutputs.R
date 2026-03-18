
VisualOutputs <- function(lvtable, variables, land_grid_list, parameters){
    ## NEED TO:
    ## -- add land image to background of gifs?
    ## -- get ldsel table from NND calc for attributes of landscapes


    # Generates output plots for examples
    library(data.table)
    library(parallel)
    library(akima)
    library(rgl)

    setDT(variables)

    core.dir <- './Output/'
#     core.dir <- 'C:/Users/eric.sodja/Documents/testrun_outputs/Output/'

    # make things into data.tables and name columns for sanity
    mapp.to.dt <- function(x, dt){
        out <- data.table(matrix(unlist(dt[,x]), ncol=length(dt[,x]), byrow=FALSE))
        setnames(out, names(dt[,x]))
    }
    print('Reading in timestep demographics...')
    tm.mat <- mcmapply(function(v, l, r){
        filenm <- paste0('tm.mat_r', r,'_l',l, '_v',v,'.csv')
        if (filenm %in% list.files(paste0(core.dir, 'tm.mat/'))){
            filenm <- paste0(paste0(core.dir, 'tm.mat/', filenm))
            tm.mat <- fread(filenm)
            setnames(tm.mat, unlist(lapply(strsplit(names(tm.mat), 'tm.mat.'), function(x) unlist(x)[2])))
            tm.mat <- tm.mat[!is.na(timestep),]
            return(tm.mat)
        }
    }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])

    tm.mat <- rbindlist(lapply(1:ncol(tm.mat), mapp.to.dt, dt=tm.mat))
    print('Reading in summary statistics...')
    summ.vals <- mcmapply(function(v, l, r){
#         summ.vals <- as.data.table(out.list["summ.vals"])
        filenm <- paste0('summ.vals_r', r,'_l',l, '_v',v,'.csv')
        if (filenm %in% list.files(paste0(core.dir, 'summ.vals/'))){
#             print(filenm)
            filenm <- paste0(paste0(core.dir, 'summ.vals/', filenm))
            summ.vals <- fread(filenm)
            setnames(summ.vals, unlist(lapply(strsplit(names(summ.vals), 'summ.vals.'), function(x) unlist(x)[2])))
            summ.vals <- summ.vals[!is.na(endtime)]
            return(summ.vals)
        }
    }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
    summ.vals <- rbindlist(lapply(1:ncol(summ.vals), mapp.to.dt, dt=summ.vals))
    print('Reading in incidence measures...')
    incidence <- mcmapply(function(v, l, r){
        filenm <- paste0('incidence_r', r,'_l',l, '_v',v,'.csv')
        if (filenm %in% list.files(paste0(core.dir, 'incidence/'))){
            filenm <- paste0(paste0(core.dir, 'incidence/', filenm))
            incidence <- fread(filenm)
            setnames(incidence, unlist(lapply(strsplit(names(incidence), 'incidence.'), function(x) unlist(x)[2])))
            incidence[, max.time := max(timestep), by=.(var, rep, land)]
            return(incidence)
        }
    }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
    incidence <- rbindlist(lapply(1:ncol(incidence), mapp.to.dt, dt=incidence))

    ## for later
#     detections <- rbindlist(mcmapply(function(v, l, r){
# #         detections <- as.data.table(out.list["detections"])
#         detections <- fread(paste0(paste0(core.dir, 'detections/detections_r', r,'_l',l, '_v',v,'.csv'))
#         setnames(detections, unlist(lapply(strsplit(names(detections), 'detections.'), function(x) unlist(x)[2])))
#         detections[, max.time := 0]
#         if(nrow(detections) > 0) {
#             detections[,max.time := max(time), by=.(var, rep, land)]
#             unq.det <- unique(detections[, .(var, land, rep, timestep, loc, max.time, code, detected)])
#         }
#         return(detections)
#     }
#         allzones <- as.data.table(out.list["allzone"])
#         setnames(allzones, unlist(lapply(strsplit(names(allzones), 'allzone.'), function(x) unlist(x)[2])))

    print('Reading in sounderlocations...')
    solocs.all <- mcmapply(function(v, l, r){
        filenm <- paste0('solocs.all_r', r,'_l',l, '_v',v,'.csv')
        if (filenm %in% list.files(paste0(core.dir, 'solocs.all/'))){
            filenm <- paste0(paste0(core.dir, 'solocs.all/', filenm))
            solocs.all <- fread(filenm)
            setnames(solocs.all, unlist(lapply(strsplit(names(solocs.all), 'solocs.all.'), function(x) unlist(x)[2])))
            solocs.all[,nlive := S+E+I+R]
            return(solocs.all)

        }
    }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
    solocs.all <- rbindlist(mclapply(1:ncol(solocs.all), mapp.to.dt, dt=solocs.all))


    print('Calculating wavespeed metrics...')
    edge.interactions <- unique(solocs.all[x==0.25 | y==0.25 | x==99.75 | y==99.75, .(v,l,r,time)])[,edge := 0]

    wv.speed <- wave_speed(solocs.all)

    print('Reading in selected land attributes...')
    ldsel.nnd <- readRDS('./Landscape_Setup/NND_Lands/4_Output/ldsel.rds')
    setDT(ldsel.nnd)
    ldsel.nnd[,land:=NULL]
    setnames(ldsel.nnd, 'index', 'land')


    print('Generating state variable temporal plots')
    seircz.plot(tm.mat, variables)


    ## wavespeed stuff
    print('Plotting wave speed metrics')
    tm.mat.edge <- wave.plot(wv.speed, tm.mat, land_grid_list, variables, edge.interactions)

    ## Proportion of simulations that ASF establishes
    print('Generating plots of establishment success...')
    unq.parms <- estab.plot(tm.mat, ldsel.nnd, variables)


##     max observed incidence for var/land combos
    print('Generating heatmaps of incidence...')
    maxinc.heat(tm.mat.edge, unq.parms, ldsel.nnd)

    ## max.dist.wave
    print('Generating heatmaps of extent of infection...')
    wvdist.plot(wv.speed, unq.parms, ldsel.nnd, wvtime=60)

    ## spd.wave
    print('Generating heatmaps of wavespeed')
    wvspd.plot(wv.speed, unq.parms, ldsel.nnd)

    ## inf.cells
    print('Generating heatmaps of infected cells')
    infcell.plot(wv.speed, unq.parms, ldsel.nnd)

    ## create map image of selected lands

    print('Generating national map of selected lands...')
    sel.lands(ldsel.nnd)

    ## landscape attributes
    # land_grid_list
    # 1: number of cells
    # 2: table:
    #     1: cell id
    #     2: dim edge low
    #     3: dim edge high
    #     4: dim edge low
    #     5: dim edge high
    #     6: center point x?
    #     7: center point y?
    #     8: cell preference value
#     lapply(seq(length(land_grid_list)), function(x){
#         in.rast <- rast(land_grid_list[[x]][[3]], type='xyz')
#         corr <- autocor(in.rast, global=TRUE)
#         return(corr)
#     })
}
