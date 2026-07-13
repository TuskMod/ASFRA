# calculates the wave speed after simulations are done

wave_speed <- function(solocs.all, cent){
    library(data.table)
    # columns are time, cell, pref, S, E, I, R, C, Z

    # get introduction point
#     setnames(solocs.all, c('ctX', 'ctY'), c('x', 'y'))
    solocs.all <- solocs.all[cent, on='cell']
    setorder(solocs.all, time, cell)
#     init.pts <- solocs.all[E > 0 | I > 0 | C > 0, .SD[1]]
    solocs.all <- solocs.all[E > 0 | I > 0 | C > 0, ]
#     setnames(init.pts, c('x', 'y'), c('xi', 'yi'))

    # get incidence by cell by time
    # time during infection should be contiguous b/c simulations end if there are no E, I, or C in the population
#     solocs.all[, tot.inf := E + I + C]

    # grab the distance of each cell from the origin (could be cached for speed)
    # connect initial points to infected points
#     solocs.all[,icell := init.pts[,cell]][,xi := init.pts[,xi]][,yi := init.pts[,yi]]
#     solocs.all <- unique(solocs.all[, .(time, tot.inf, cell, x, y)][init.pts[, .(cell, xi, yi)], on=.(v, l, r)])
#     solocs.all[, dist := sqrt((x - xi)^2 + (y - yi)^2)]

    # get mean, median, sd of infected cell distances
    edge.zone.time <- solocs.all[,max(dist), by=time][round(V1, 2) <= 50.25, .N]+1
    solocs.all <- solocs.all[time <= edge.zone.time]
    solocs.all[, avg := mean(dist), by=.(time)]
    solocs.all[, var := var(dist),  by=.(time)]
    solocs.all[, max := max(dist),  by=.(time)]
    solocs.all[, min := min(dist),  by=.(time)]
    solocs.all[, prop10plus := length(which(dist>10))/.N, by=.(time)]

    # get speed of average distance
    avg.dist.diff <- unique(solocs.all[,.(time, avg)])
    avg.dist.diff[, avg.dist.diff := avg - data.table::shift(avg, 1),]
    # na values in shifted average distances at time 1 are from not having time 0 data to shift from, so 0 makes sense
    avg.dist.diff[is.na(avg.dist.diff) & time == 1, avg.dist.diff := 0]

    # connect average wave speed to output data table
    solocs.all <- solocs.all[avg.dist.diff, on=.NATURAL]

#     solocs.all[,edge := 0]
#     solocs.all[x==0.25 | y==0.25 | x==99.75 | y==99.75, edge := 1] # alternative: first to reach distance for closest edge to point of introduction (50km)
#     solocs.all[, vedge := cumsum(edge), by=.(time)]
#     solocs.all <- solocs.all[vedge==0,]
#
#     solocs.all[,vedge := NULL]
#     solocs.all[,edge := NULL]
#     solocs.all[,xi := NULL]
#     solocs.all[,yi := NULL]


    return(solocs.all)
}
