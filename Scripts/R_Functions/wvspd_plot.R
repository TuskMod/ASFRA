## wavespeed.plot over time

wvspd.plot <- function(wv.speed, unq.parms, ldsel.nnd){
    ## spd.wave
    wvspd.avgdistdiff <- unique(wv.speed[is.na(edge), mean(avg.dist.diff), by=.(v,l)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(wvspd.avgdistdiff, 'V1','avg.spd')
    wvspd.avgdistdiff[,avg.spd.norm := avg.spd/max(avg.spd, na.rm=TRUE)]
    png('./Output/figures/land_parcombos_wvspd.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=wvspd.avgdistdiff[!is.na(avg.spd.norm) & v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(avg.spd.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()
}
