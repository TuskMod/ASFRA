## Maximum distance of the wave before it hits the edge of the simulation space/at time 60 weeks (could be variable)

wvdist.plot <- function(wv.speed, unq.parms, ldsel.nnd, wvtime=60){
    ## dist.wave
    wvdist.t60 <- unique(wv.speed[is.na(edge), maxtime := max(time), by=.(v,l)][time <=wvtime,][time == wvtime | time == maxtime, .(v,l,r,time, max)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]

    # find ranges of values to keep all plots on the same scale
    wvdist.t60[,disp:=as.numeric(disp)]
    x0 <- wvdist.t60[,disp]
    y0 <- wvdist.t60[,nnd_med]
    z0 <- wvdist.t60[,max]

    hl <- chull(x0, y0)
    hl.tab <- data.table(x=x0[hl], y=y0[hl])
    hl.tab[,x1 := shift(x,1, type='cyclic')]
    hl.tab[,y1 := shift(y,1, type='cyclic')]
    int.pts0 <- interp(x0, y0, z0, xo=seq(min(x0), max(x0), length=200), yo=seq(min(y0), max(y0), length=200), linear=FALSE, extrap=TRUE, duplicate='median')

    xmin <- round(min(int.pts0$x)-50, -2)
    xmax <- round(max(int.pts0$x)+50, -2)
    ymin <- round(min(int.pts0$y)-50, -2)
    ymax <- round(max(int.pts0$y)+50, -2)
    dist.min <- 0
#     dist.max <- round(max(int.pts0$z, na.rm=TRUE)+0.3, 1)
    dist.max <- round(1.5*max(z0))

    png('./Output/figures/maxdist.png', width=1000, height=800)
    par(mfrow=c(3,3), oma=c(3,3,2,0))
    lapply(seq(nrow(unq.parms)), function(row){
        sub.parms <- unq.parms[row,]
        sub.dat <- wvdist.t60[v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]]

        # define inputs for interpolation
        x <- as.numeric(sub.dat[,disp])
        y <- sub.dat[,nnd_med]
        z <- sub.dat[,max]

        # interpolate
        int.pts <- interp(x, y, z, xo=seq(min(x), max(x), length=200), yo=seq(min(y), max(y), length=200), linear=FALSE, extrap=TRUE, duplicate='median')

        par(mar=c(2,2,1,1))
        if (row %in% c(1,4,7)){
            par(mar=c(2,5,1,1))

        }
        if (row %in% c(7:8)){
            par(mar=c(5,2,1,1))
        }
        if (row == 7) par(mar=c(5,5,1,1))

        # plot contours within range of input data
        range.val <- range(int.pts$x, na.rm=TRUE)
        if(abs(range.val[1]-range.val[2]) > 0) {
            plot(y ~ x, data=int.pts, col=NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), axes=FALSE, xlab='', ylab='')
            .filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, levels=pretty(c(dist.min,dist.max), 20), col=hcl.colors(20, 'Roma', rev=TRUE))
            points(y ~ x, pch=10, cex=3, lwd=2)
            contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(dist.min, dist.max))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=1)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
            polygon(hl.tab, col='darkgrey', lwd=2, lty=3, density=0)

        } else { # if there is only one level for contour() lines, throws an error
            plot(y ~ x, data=int.pts, col=NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), axes=FALSE, xlab='', ylab='')
            .filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, levels=pretty(c(dist.min,dist.max), 20), col=hcl.colors(20, 'Roma', rev=TRUE))
            points(y ~ x, pch=10, cex=3, lwd=2)
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=1)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
            polygon(hl.tab, col='darkgrey', lwd=2, lty=3, density=0)
        }

        if (row %in% c(1,4,7)){
            mtext('Median NND km', side=2, outer=FALSE, line=3, font=1, las=3)

        }
        if (row %in% c(6:8)){
            mtext('Sounder dispersal km', side=1, outer=FALSE, line=3, font=1)
        }

        text(x=850, y=8000, labels=row, cex=2.5, font=2)
    })
    mtext('Maximum distance of epidemic spread', side=3, outer=TRUE, line=0.8, font=2)

    colnms <- hcl.colors(20, 'Roma', rev=TRUE)
    xl <- 1
    yb <- 1
    xr <- 1.1
    yt <- 2
    par(mar=c(5.1,1.5,4.1,0.5))
    plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
    rect(
        xl,
        head(seq(yb,yt,(yt-yb)/20),-1),
        xr,
        tail(seq(yb,yt,(yt-yb)/20),-1),
        col=colnms
    )

    mtext(round(seq(dist.min,dist.max,length=20),2),side=2,at=tail(seq(yb,yt,(yt-yb)/20),-1)-0.02,las=2,cex=0.7)
    legend('center', c('original land point','boundary of original land points'), pch=c(10,NA), pt.cex=c(3, 1), lty=c(NA, 3), col=c(1,'darkgray'), lwd=c(0, 2))
    dev.off()

#         # plot contours within range of input data
#         filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, xlim=c(xmin, xmax), ylim=c(ymin, ymax), zlim=c(dist.min, dist.max), plot.axes={
#             points(y ~ x, pch=10)
#             contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(dist.min, dist.max))
#             axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=2)
#             axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
#         })
#
#         mtext('Wave distance maximum t=60', side=3, outer=FALSE, line=2, font=2)
#         mtext('Sounder dispersal km', side=1, outer=FALSE, line=3, font=1)
#         mtext('Median NND km', side=2, outer=FALSE, line=3, font=1, las=3)
#
#         dev.off()
#     })



#     png(paste0('./Output/figures/land_parcombos_maxdist_t',wvtime,'.png'), width=1000, height=800)
#     par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
#     lapply(seq(nrow(unq.parms)), function(y){
#         sub.parms <- unq.parms[y,]
#         plot(nnd_med ~ disp, data=wvdist.t60[!is.na(max.norm) & v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
#              pch=15, col=rgb(max.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
#     })
#     dev.off()
}
