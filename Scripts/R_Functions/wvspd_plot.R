## wavespeed.plot over time

wvspd.plot <- function(wv.speed, unq.parms, ldsel.nnd){
    ## spd.wave
    wvspd.avgdistdiff <- unique(wv.speed[is.na(edge), mean(avg.dist.diff), by=.(v,l)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(wvspd.avgdistdiff, 'V1','avg.spd')
    wvspd.avgdistdiff[,disp := as.numeric(disp)]

    # find ranges of values to keep all plots on the same scale
    x0 <- wvspd.avgdistdiff[,disp]
    y0 <- wvspd.avgdistdiff[,nnd_med]
    z0 <- wvspd.avgdistdiff[,avg.spd]
    int.pts0 <- interp(x0, y0, z0, xo=seq(min(x0), max(x0), length=200), yo=seq(min(y0), max(y0), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

    xmin <- round(min(int.pts0$x)-50, -2)
    xmax <- round(max(int.pts0$x)+50, -2)
    ymin <- round(min(int.pts0$y)-50, -2)
    ymax <- round(max(int.pts0$y)+50, -2)
    spd.min <- round(min(int.pts0$z, na.rm=TRUE)-0.1, 1)
    spd.max <- round(max(int.pts0$z, na.rm=TRUE)+0.3, 1)

    lapply(seq(nrow(unq.parms)), function(row){
        png(paste0('./Output/figures/land_parcombos_wvspd',row,'.png'), width=800, height=600)
        # subset data to specific variable combinations
        sub.parms <- unq.parms[row,]
        sub.dat <- wvspd.avgdistdiff[v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density],]

        # define inputs for interpolation
        x <- as.numeric(sub.dat[,disp])
        y <- sub.dat[,nnd_med]
        z <- sub.dat[,avg.spd]

        # interpolate
        int.pts <- interp(x, y, z, xo=seq(min(x), max(x), length=200), yo=seq(min(y), max(y), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

        # plot contours within range of input data
        filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, xlim=c(xmin, xmax), ylim=c(ymin, ymax), zlim=c(spd.min, spd.max), plot.axes={
            points(y ~ x, pch=10)
            contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(spd.min, spd.max))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=2)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
        })

        mtext('Average wave speed from introduction point', side=3, outer=FALSE, line=2, font=2)
        mtext('Sounder dispersal km', side=1, outer=FALSE, line=3, font=1)
        mtext('Median NND km', side=2, outer=FALSE, line=3, font=1, las=3)

    dev.off()
    })
}
