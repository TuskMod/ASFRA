## Maximum distance of the wave before it hits the edge of the simulation space/at time 60 weeks (could be variable)

wvdist.plot <- function(wv.speed, unq.parms, ldsel.nnd, wvtime=60){
    ## max.dist.wave
    wvdist.t60 <- unique(wv.speed[is.na(edge), maxtime := max(time), by=.(v,l)][time <=wvtime,][time == wvtime | time == maxtime, .(v,l,r,time, max)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    wvdist.t60[,max.norm := max/max(max, na.rm=TRUE)]
    png(paste0('./Output/figures/land_parcombos_maxdist_t',wvtime,'.png'), width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=wvdist.t60[!is.na(max.norm) & v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(max.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()
}
