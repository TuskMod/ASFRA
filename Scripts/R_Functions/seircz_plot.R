## Temporal plots of full populations SEIRCZ

seircz.plot <- function(tm.mat, variables){

    tm.pop.unq <- unique(tm.mat[,.(var, land)])
    rows.plt <- round(sqrt(length(unique(tm.pop.unq[,var]))))
    cols.plt <- ceiling(sqrt(length(unique(tm.pop.unq[,var]))))
    for (l in unique(tm.pop.unq[,land])){
        png(paste0('./Output/figures/tm.plots_land_', l, '.png'), width=450*cols.plt, height=500*rows.plt)
        ## make plot dimensions dynamics based on parameter inputs
        pdef <- par(mfrow=c(rows.plt, cols.plt),
                    oma=c(8, 1.2, 0.2, 0.2),
                    lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
        # par(mfrow=c(1,1), oma=c(8,0.2,0.2,0.2),lwd=2.3, cex.lab=1.8, cex.axis=1.6, cex.main=2, cex.sub=1.4)
        # par(mfcol=c(length(unique(tm.mat[,var])), length(unique(tm.mat[,land]))))

        # loop over the parameter combinations for each plot panel
        mcmapply(seirczbb.temporal, v=tm.pop.unq[land==l, var], l=tm.pop.unq[land==l,land], plt.i=seq(nrow(tm.pop.unq[land==l,])), MoreArgs = list(rows.plt=rows.plt, cols.plt=cols.plt, dat=tm.mat, vlist=variables))
        # create a legend under the panels
        par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend("bottom", legend=c('S','E','I','R','C','Z','Births'),#'rep1','rep2'),
        lty=c(rep(1,8), 2),
        col=c('olivedrab','orange','red','blue','purple','black','pink'),#'grey','grey'),
        horiz=TRUE, bty='n', cex=2.3, lwd=3)
        dev.off()
        par(pdef)
    }

}

seirczbb.temporal <- function(v, l, plt.i, rows.plt, cols.plt, dat=tm.mat, vlist = variables){
    # set plot ranges based on maximum of everything that will be on the multiplot figure
    xrng = range(dat[,timestep])
    yrng = log1p(range(dat[,.(BB, S, E, I, R, C, Z)]))

    vardat <- paste(vlist[v,], collapse=' ') # quick and dirty parameter inclusion
    # xlabel default
    xlabi <- ''
    # ylabel default
    ylabi <- ''
    # default plot panel margin sizes
    bmar <- lmar <- 2
    if (plt.i / rows.plt >= rows.plt){ # first row of panels
                                      xlabi <- 'time (weeks)'
                                      bmar <- 4
                                      }
    if (plt.i %% cols.plt == 1){ # first column of panels
                                ylabi <- 'log(1+individuals)'
                                lmar <- 4
                                }
    pdef2 <- par(mar = c(bmar, lmar, 1, 1))
    # create a blank plot
    plot(0,0,xlim=xrng, ylim=yrng, col=NULL, ann=FALSE)
    mtext(xlabi, 1, line=2.5, cex=1.8)
    mtext(ylabi, 2, line=2.5, cex=1.8)
    mtext(paste('vars',v,'| land',l), 3, line=-1.5, cex=1.7, font=2)
    mtext(vardat, 3, line=-2.7)
    # subset with things that have only the land tile and variable combination
    subdat <- dat[var==v & land==l,]
    # count up reps
    reps <- unique(subdat[,rep])
    # draw a line of each type for each rep
    lapply(reps, function(r) lines(log1p(S) ~ timestep, data=subdat[rep==r | rep == 0,], col='olivedrab', lty=1))
    lapply(reps, function(r) lines(log1p(E) ~ timestep, data=subdat[rep==r | rep == 0,], col='orange', lty=1))
    lapply(reps, function(r) lines(log1p(I) ~ timestep, data=subdat[rep==r | rep == 0,], col='red', lty=1))
    lapply(reps, function(r) lines(log1p(R) ~ timestep, data=subdat[rep==r | rep == 0,], col='blue', lty=1))
    lapply(reps, function(r) lines(log1p(C) ~ timestep, data=subdat[rep==r | rep == 0,], col='purple', lty=1))
    lapply(reps, function(r) lines(log1p(Z) ~ timestep, data=subdat[rep==r | rep == 0,], col='black', lty=1))
    lapply(reps, function(r) lines(log1p(BB) ~ timestep, data=subdat[rep==r | rep == 0,], col='pink', lty=1))
    lapply(reps, function(r) abline(v=subdat[rep==r,max(timestep)], col='grey'))
    lapply(reps, function(r) abline(v=subdat[rep==0,max(timestep)], col='black', lty=3))
    par(pdef2)
}
