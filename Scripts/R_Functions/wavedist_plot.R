## Maximum distance of the wave before it hits the edge of the simulation space/at time 60 weeks (could be variable)

wvdist.plot <- function(wv.speed, unq.parms, ldsel.nnd, wvtime=60){
    ## dist.wave
    wvdist.t60 <- unique(wv.speed[is.na(edge), maxtime := max(time), by=.(v,l)][time <=wvtime,][time == wvtime | time == maxtime, .(v,l,r,time, max)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]


#     browser()
    # find ranges of values to keep all plots on the same scale
    wvdist.t60[,disp:=as.numeric(disp)]
    x0 <- wvdist.t60[,disp]
    y0 <- wvdist.t60[,nnd_med]
    z0 <- wvdist.t60[,max]

    int.pts0 <- interp(x0, y0, z0, xo=seq(min(x0), max(x0), length=200), yo=seq(min(y0), max(y0), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

    xmin <- round(min(int.pts0$x)-50, -2)
    xmax <- round(max(int.pts0$x)+50, -2)
    ymin <- round(min(int.pts0$y)-50, -2)
    ymax <- round(max(int.pts0$y)+50, -2)
    dist.min <- 0
#     dist.min <- round(min(int.pts0$z, na.rm=TRUE)-0.1, 1)
    dist.max <- round(max(int.pts0$z, na.rm=TRUE)+0.3, 1)

    lapply(seq(nrow(unq.parms)), function(row){
        png(paste0('./Output/figures/land_parcombos_maxdist',row,'.png'), width=1000, height=800)
        sub.parms <- unq.parms[row,]
        sub.dat <- wvdist.t60[v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]]


        # define inputs for interpolation
        x <- as.numeric(sub.dat[,disp])
        y <- sub.dat[,nnd_med]
        z <- sub.dat[,max]

        # interpolate
        int.pts <- interp(x, y, z, xo=seq(min(x), max(x), length=200), yo=seq(min(y), max(y), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

        # plot contours within range of input data
        filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, xlim=c(xmin, xmax), ylim=c(ymin, ymax), zlim=c(dist.min, dist.max), plot.axes={
            points(y ~ x, pch=10)
            contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(dist.min, dist.max))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=2)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
        })

        mtext('Wave distance maximum t=60', side=3, outer=FALSE, line=2, font=2)
        mtext('Sounder dispersal km', side=1, outer=FALSE, line=3, font=1)
        mtext('Median NND km', side=2, outer=FALSE, line=3, font=1, las=3)

        dev.off()
    })



#     png(paste0('./Output/figures/land_parcombos_maxdist_t',wvtime,'.png'), width=1000, height=800)
#     par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
#     lapply(seq(nrow(unq.parms)), function(y){
#         sub.parms <- unq.parms[y,]
#         plot(nnd_med ~ disp, data=wvdist.t60[!is.na(max.norm) & v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
#              pch=15, col=rgb(max.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
#     })
#     dev.off()
}
