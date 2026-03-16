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


    #     established[, cat.name := paste(land, state, variant, density, sep='_')]
    #     est.wide.pct <- dcast(established, var ~ land, value.var='est.pct')
    #     est.wide.no <- dcast(established, var ~ land, value.var='est.no')
    png('./Output/figures/PCt_estab.png', width=1000, height=800)
    barplot(est ~ as.factor(land) + as.factor(var), data=established, beside=TRUE, legend=TRUE, names=established[,paste(contact, variant, density, sep='_')], ylim=c(0,1), main='% Establishment by land', xlab="Movement_Strain_Density", ylab='% established')
    dev.off()
    ## will want to have this as heatmap for landscape attributes for each of the 8 variable combinations

    unq.parms <- unique(est.lands[,.(var, contact, variant, density)])
    png('./Output/figures/land_parcombos_estab.png', width=1000, height=800)
    par(mfrow=c(length(unique(unq.parms[,density])), nrow(unq.parms)/length(unique(unq.parms[,density]))))
    lapply(seq(nrow(unq.parms)), function(y){
        sub.parms <- unq.parms[y,]
        plot(nnd_med ~ disp, data=est.lands[var == sub.parms[,var] & contact == sub.parms[,contact] & variant == sub.parms[,variant] & density == sub.parms[,density]],
             pch=15, col=rgb(est,0,0), cex=4, main=paste(unlist(sub.parms), collapse='_'))
    })
    dev.off()

    return(unq.parms)
}
