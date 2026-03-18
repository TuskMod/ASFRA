## Proportion of simulation runs that establish for a given parameter/land combination

estab.plot <- function(tm.mat, ldsel.nnd, variables){

    # takes anything with (1) 2 or more cases at a given time, OR (2) at least 10 cases total
    ## should be establishment criteria from madison's paper
    tm.mat[,max.inf := max(I), by=.(var, land, rep)]
    tm.mat[,tot.exp := sum(E), by=.(var, land, rep)]
    tm.mat[,est := 0]
    tm.mat[max.inf > 1 & tot.exp >= 10, est := 1]
    #     established <-  unique(tm.mat[,est.no := sum(est, na.rm=TRUE), by=.(var, land)][,est.pct := est.no/.N, by=.(var, land)][,.(var, land, est.no, est.pct)])
    established <- unique(unique(tm.mat[,.(var, land, rep, est)])[,est := sum(est)/.N, by=.(var, land)][,-'rep'])
    established <- established[variables[,.(contact, variant, density, id)], on=.(var=id)]
    est.lands <- established[ldsel.nnd, on=.(land=land), nomatch=NULL]


    #     established[, cat.name := paste(land, state, variant, density, sep='_')]
    #     est.wide.pct <- dcast(established, var ~ land, value.var='est.pct')
    #     est.wide.no <- dcast(established, var ~ land, value.var='est.no')
#     png('./Output/figures/PCt_estab.png', width=1000, height=800)
#     barplot(est ~ as.factor(land) + as.factor(var), data=established, beside=TRUE, legend=TRUE, names=established[,paste(contact, variant, density, sep='_')], ylim=c(0,1), main='% Establishment by land', xlab="Movement_Strain_Density", ylab='% established')
#     dev.off()
#     ## will want to have this as heatmap for landscape attributes for each of the 8 variable combinations

    unq.parms <- unique(est.lands[,.(var, contact, variant, density)])


    est.lands[,disp := as.numeric(disp)]
    # find ranges of values to keep all plots on the same scale
    x0 <- est.lands[,disp]
    y0 <- est.lands[,nnd_med]
    z0 <- est.lands[,est]
    int.pts0 <- interp(x0, y0, z0, xo=seq(min(x0), max(x0), length=200), yo=seq(min(y0), max(y0), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

    xmin <- round(min(int.pts0$x)-50, -2)
    xmax <- round(max(int.pts0$x)+50, -2)
    ymin <- round(min(int.pts0$y)-50, -2)
    ymax <- round(max(int.pts0$y)+50, -2)

    lapply(seq(nrow(unq.parms)), function(row){
        png(paste0('./Output/figures/land_parcombos_estab',row,'.png'), width=1000, height=800)
        sub.parms <- unq.parms[row,]
        sub.dat <- est.lands[var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]]

        # define inputs for interpolation
        x <- as.numeric(sub.dat[,disp])
        y <- sub.dat[,nnd_med]
        z <- sub.dat[,est]

        # interpolate
        int.pts <- interp(x, y, z, xo=seq(min(x), max(x), length=200), yo=seq(min(y), max(y), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

        # plot contours within range of input data
        if(diff(range(int.pts$z, na.rm=TRUE) > 0)) {
        filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, xlim=c(xmin, xmax), ylim=c(ymin, ymax), zlim=c(0, 1), plot.axes={
            points(y ~ x, pch=10)
            contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(0, 1))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=2)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
        })} else {
        filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, xlim=c(xmin, xmax), ylim=c(ymin, ymax), zlim=c(0, 1), plot.axes={
            points(y ~ x, pch=10)
#             contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(0, 1))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=2)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
        })}

        mtext('Proportion of simulated epidemics established', side=3, outer=FALSE, line=2, font=2)
        mtext('Sounder dispersal km', side=1, outer=FALSE, line=3, font=1)
        mtext('Median NND km', side=2, outer=FALSE, line=3, font=1, las=3)

        dev.off()
    })

    return(unq.parms)
}
