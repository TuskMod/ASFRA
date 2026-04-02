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

    unq.parms <- unique(est.lands[,.(var, contact, variant, density)])

    est.lands[,disp := as.numeric(disp)]
    # find ranges of values to keep all plots on the same scale
    x0 <- est.lands[,disp]
    y0 <- est.lands[,nnd_med]
    z0 <- est.lands[,est]
    hl <- chull(x0, y0)
    hl.tab <- data.table(x=x0[hl], y=y0[hl])
    hl.tab[,x1 := shift(x,1, type='cyclic')]
    hl.tab[,y1 := shift(y,1, type='cyclic')]
    int.pts0 <- interp(x0, y0, z0, xo=seq(min(x0), max(x0), length=200), yo=seq(min(y0), max(y0), length=200), linear=FALSE, extrap=FALSE, duplicate='median')

    xmin <- round(min(int.pts0$x)-50, -2)
    xmax <- round(max(int.pts0$x)+50, -2)
    ymin <- round(min(int.pts0$y)-50, -2)
    ymax <- round(max(int.pts0$y)+50, -2)

    png('./Output/figures/estab.png', width=1000, height=800)
    par(mfrow=c(3,3), oma=c(3,3,2,0))
    lapply(seq(nrow(unq.parms)), function(row){
        sub.parms <- unq.parms[row,]
        sub.dat <- est.lands[var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]]

        # define inputs for interpolation
        x <- as.numeric(sub.dat[,disp])
        y <- sub.dat[,nnd_med]
        z <- sub.dat[,est]

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
        if(diff(range(int.pts$z, na.rm=TRUE) > 0)) {
            plot(y ~ x, data=int.pts, col=NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), axes=FALSE, xlab='', ylab='')
            .filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, levels=pretty(c(0,1), 20), col=hcl.colors(20, 'Roma', rev=TRUE))
            points(y ~ x, pch=10, cex=3, lwd=2)
            contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, nlevels=20, add=TRUE, zlim=c(0, 1))
            axis(1, seq(xmin, xmax, by=100), labels = seq(xmin, xmax, by=100)/1000, las=1)
            axis(2, seq(ymin, ymax, by=500), labels = seq(ymin, ymax, by=500)/1000, las=2)
            polygon(hl.tab, col='darkgrey', lwd=2, lty=3, density=0)

        } else { # if there is only one level for contour() lines, throws an error
            plot(y ~ x, data=int.pts, col=NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), axes=FALSE, xlab='', ylab='')
            .filled.contour(x=int.pts$x, y=int.pts$y, z=int.pts$z, levels=pretty(c(0,1), 20), col=hcl.colors(20, 'Roma', rev=TRUE))
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
    mtext('Proportion of simulated epidemics established', side=3, outer=TRUE, line=0.8, font=2)

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

    mtext(round(seq(0,1,length=20),2),side=2,at=tail(seq(yb,yt,(yt-yb)/20),-1)-0.02,las=2,cex=0.7)
    legend('center', c('original land point','boundary of original land points'), pch=c(10,NA), pt.cex=c(3, 1), lty=c(NA, 3), col=c(1,'darkgray'), lwd=c(0, 2))
    dev.off()

    return(unq.parms)
}
