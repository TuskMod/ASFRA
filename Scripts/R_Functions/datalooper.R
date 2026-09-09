
## inputs -- names of tm.mat files, table of cell number and centroid points
## cent = land_grid_list[[1]][[1]][[2]][,c(1,6,7)] for cell, x, y

data.looper <- function(infile, cent, lvtable){
    library(data.table)
    # Grab the files with the infile name/numbers
    tm.mat <- fread(paste0('./Output/tm.mat/', infile), select=c('E', 'I'))

    vlr_end <- unlist(tstrsplit(infile, 'tm.mat_', keep=2))
    solocs.all <- fread(paste0('./Output/solocs.all/solocs.all_', vlr_end), select=c('time','cell','E','I','C'))

    if (solocs.all[time == 1 & (E > 0 | I > 0 | C > 0), cell] != 20101){
        # if not default introduction point, get the introduction point
        inf.pt <- solocs.all[time == 1 & (E > 0 | I > 0 | C > 0),]
        xy.pt <- cent[cell==inf.pt[,cell]]
        cent[,dist := sqrt((x-xy.pt[,x])^2 + (y-xy.pt[,y])^2)]
    }

    # get the vlr values
    vlr <- tstrsplit(vlr_end, '[v|l|r]+|_|.gz', keep=c(2,4,6))
    vval <- as.numeric(vlr[[3]])
    lval <- as.numeric(vlr[[2]])
    rval <- as.numeric(vlr[[1]])
browser()
    # get wavespeed metrics (also trim data to a distance equal to the closest edge to introduction point)
    eic.edge <- wave_speed(solocs.all, cent)

    # do things with the files
    # (1) establishment proportion
    max.inf <- tm.mat[,max(I)]
    tot.exp <- tm.mat[,sum(E)]
    est <- 0 # failed to establish
    if (max.inf > 1 & tot.exp >= 10) est <- 2 # established
    # (3) max spread distance
    max.dist <- eic.edge[,max(dist)]
    if (max(eic.edge[,time]) < 78 & max.dist < 50 & est == 2) est <- 1 # fizzled out part way

    # 5% escape 10km radius
    tm.escape <- eic.edge[order(time, dist)][prop10plus >= 0.05, .SD[1]][,time]

    # (4) total area of epidemic
    inf.area <- length(eic.edge[,unique(cell)])
    # (5) speed of epidemic spread
    inf.spd <- eic.edge[,unique(avg.dist.diff), by=time][,mean(V1)]

    # sounder-weeks of infection (Exclude C) over the first 20 weeks
    sounder.weeks <- eic.edge[E > 0 | I > 0,.N]

    # cells inside of max infection distance
    infpot.cells <- cent[dist <= max.dist,]
    prop.infd <- inf.area/nrow(infpot.cells)

    # incidence
#     eic.edge[, deltaE := E-data.table::shift(E,1)][is.na(deltaE) & timestep==1, deltaE := 0]
#     eic.edge[,incidence := deltaE/S]

    max.inc <- max(tm.mat[1:max(eic.edge[,time]), E - data.table::shift(E,1)], na.rm=TRUE)


    outs <- data.table(v = vval, l = lval, r = rval,
                est = est,
                edge.tm = max(eic.edge[,time]),
                max.dist = max.dist,
                inf.area = inf.area,
                inf.spd = inf.spd,
                sounder.weeks = sounder.weeks,
#                 infpot.cells = nrow(infpot.cells),
                prop.infd = prop.infd,
                max.inc = max.inc,
                tm.esc = tm.escape)

    return(outs)
}


## data.looper('tm.mat_r5_l65_v5.gz', tar_read(land_grid_list)[[1]][[1]][[2]][,c(1,6,7)])
