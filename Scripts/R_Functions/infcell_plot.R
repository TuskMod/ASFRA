## infected cells over time

infcell.plot <- function(wv.speed, unq.parms, ldsel.nnd){
    ## inf.cells
    max.inf.cells <- wv.speed[is.na(edge), max(tot.inf.cells), by=.(v,l)][unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(max.inf.cells, 'V1', 'max.inf.cells')
    max.inf.cells[,max.inf.cells.norm := max.inf.cells/max(max.inf.cells, na.rm=TRUE)]
    png('./Output/figures/land_parcombos_infcells.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=max.inf.cells[!is.na(max.inf.cells.norm) & v == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(max.inf.cells.norm,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

}
