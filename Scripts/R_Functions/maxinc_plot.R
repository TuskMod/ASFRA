## establishment heatmap

maxinc.heat <- function(tm.mat.edge, unq.parms, ldsel.nnd){
    max.incidence <- tm.mat.edge[, max(incidenceEI, na.rm=TRUE), by=.(var, land)][unq.parms, on='var'][ldsel.nnd, on='land']
    setnames(max.incidence, 'V1','max.inc')
    max.incidence[,disp := as.numeric(disp)]

    # find ranges of values to keep all plots on the same scale
    x0 <- max.incidence[,disp]
    y0 <- max.incidence[,nnd_med]
    z0 <- max.incidence[,max.inc]

    int.pts0 <- interp(x0, y0, z0, xo=seq(min(x0), max(x0), length=200), yo=seq(min(y0), max(y0), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

    xmin <- round(min(int.pts0$x)-50, -2)
    xmax <- round(max(int.pts0$x)+50, -2)
    ymin <- round(min(int.pts0$y)-50, -2)
    ymax <- round(max(int.pts0$y)+50, -2)
    inc.min <- round(min(int.pts0$z, na.rm=TRUE)-0.1, 1)
    inc.max <- round(max(int.pts0$z, na.rm=TRUE)+0.3, 1)

    lapply(seq(nrow(unq.parms)), function(row){
        png(paste0('./Output/figures/land_parcombos_maxincd',row,'.png'), width=1000, height=800)
        sub.parms <- unq.parms[row,]
        sub.dat <- max.incidence[var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]]


        # define inputs for interpolation
        x <- as.numeric(sub.dat[,disp])
        y <- sub.dat[,nnd_med]
        z <- sub.dat[,max.inc]

        # interpolate
        int.pts <- interp(x, y, z, xo=seq(min(x), max(x), length=200), yo=seq(min(y), max(y), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

        # plot contours within range of input data
        filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, xlim=c(xmin, xmax), ylim=c(ymin, ymax), zlim=c(inc.min, inc.max), plot.axes={
            points(y ~ x, pch=10)
            contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(inc.min, inc.max))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=2)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
        })

        mtext('Incidence maximum', side=3, outer=FALSE, line=2, font=2)
        mtext('Sounder dispersal km', side=1, outer=FALSE, line=3, font=1)
        mtext('Median NND km', side=2, outer=FALSE, line=3, font=1, las=3)

        dev.off()
    })
}
