## establishment heatmap

maxinc.heat <- function(tm.mat.edge, unq.parms, ldsel.nnd){
    max.incidence <- tm.mat.edge[, max(incidenceEI, na.rm=TRUE), by=.(var, land)][unq.parms, on='var'][ldsel.nnd, on='land']
    setnames(max.incidence, 'V1','max.inc')
    max.incidence[, max.inc.norm := max.inc/max(max.inc, na.rm=TRUE)]
    png('./Output/figures/land_parcombos_maxincd.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=max.incidence[!is.na(max.inc.norm) & var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
            pch=15, col=rgb(max.inc.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()
}
