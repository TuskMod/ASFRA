## plotting wavespeed metrics

wave.plot <- function(wv.speed, tm.mat, land_grid_list, variables, edge.interactions){
#     edge.interactions <- unique(wv.speed[x==0.25 | y==0.25 | x==99.75 | y==99.75, .(v,l,r,time)])[,edge := 0]
    tm.mat.edge <- edge.interactions[tm.mat, on=.(v=var, l=land, r=rep, time=timestep)]
    setnames(tm.mat.edge, c('v','l','r','time'), c('var','land','rep','timestep'))
    if(any(tm.mat.edge[,is.na(edge)])) {
        tm.mat.edge[is.na(edge), seq.test := 1:.N, by=.(var,land,rep)]
        tm.mat.edge <- tm.mat.edge[timestep == seq.test,][,seq.test := NULL]
    }

    #     tm.mat.edge[,EIC := E+I+C]
    #     tm.mat.edge[,dEIC := EIC-data.table::shift(EIC,1), by=.(var,land,rep)]
    #     tm.mat.edge[is.na(dEIC) & timestep==1,dEIC := 0]
    #     tm.mat.edge[,incidence := dEIC/S]
    #     max.inc <- tm.mat.edge[,max(incidence), by=.(var,land,rep)]
    #     matplot((dcast(tm.mat.edge, timestep ~ var + land + rep, value.var = 'incidence', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]), type='l', lty=1)

    tm.mat.edge[,EI := E+I]
    tm.mat.edge[,dEI := EI-data.table::shift(EI,1), by=.(var,land,rep)]
    tm.mat.edge[is.na(dEI) & timestep==1,dEI := 0]
    tm.mat.edge[,incidenceEI := dEI/S]
    #     matplot((dcast(tm.mat.edge, timestep ~ var + land + rep, value.var = 'incidenceEI', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]), type='l', lty=1, main="Incidence", ylab="(E+I_t - E+I_t-1)/S_t", xlab='Week')

    max.inc <- tm.mat.edge[,max(incidenceEI), by=.(var,land,rep)]
    max.inc <- max.inc[variables[,id := 1:.N], on=.(var=id)]
    #     boxplot(V1 ~ contact + variant + as.factor(land) + as.factor(density), data=max.inc, ylab='maximum incidence')

    ## wavespeed stuff
    ## 4 plots for temporal dynamics -- incidence/time, max distance/time, peak speed, area/time
    png('./Output/figures/wavespeed.png', width=1200, height=1400)
    par.og <- par(mfrow=c(2,2), oma=c(8,1.2,0.2,0.2), mar=c(4,4,1,1))

    tm.mat.inc <- dcast(tm.mat.edge, timestep ~ var + land, value.var = 'incidenceEI', fun.aggregate=mean, na.rm=TRUE, fill=NA)[,-c(1)]
    matplot(tm.mat.inc, type='b', lty=1, main="Incidence", ylab="(E+I_t - E+I_t-1)/S_t", xlab='Week', col=unlist(tstrsplit(names(tm.mat.inc), '_', keep=1)), pch=as.numeric(unlist(tstrsplit(names(tm.mat.inc), '_', keep=2))))

    #     avg.wave <- dcast(wv.speed, time~ v+l, value.var='avg', fun.aggregate=mean, fill=NA)
    #     matplot(avg.wave[,2:ncol(avg.wave)], type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
                #             ylab='avg dist from source', xlab='week', main='Avg dist from introduction')

    max.dist.wave <- dcast(wv.speed, time~ v+l, value.var='max', fun.aggregate=max, fill=NA)
    yrng <- range(max.dist.wave[,2:ncol(max.dist.wave)], na.rm=TRUE)
    matplot(max.dist.wave[,2:ncol(max.dist.wave)], type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
            ylab='max dist from source', xlab='week', main='Max dist from introduction', ylim=yrng)
    #     matplot(max.dist.wave.edge[,2:ncol(max.dist.wave.edge)], type='l', lty=2, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)), add=TRUE)

    spd.wave <- dcast(unique(wv.speed[is.na(edge),.(v,l,r,time,avg.dist.diff)]), time~ v+l, value.var='avg.dist.diff', fun.aggregate=mean, fill=NA)
    matplot(spd.wave[,2:ncol(spd.wave)], type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
            ylab='km/week', xlab='week', main='Avg. wave speed from introduction point')
    abline(h=0, col='grey')

    #     pk.dist <- dcast(wv.speed[order(v,l,r,time,-tot.inf)][,.SD[1], by=.(v,l,r,time)], time~ v+l, value.var=c('tot.inf','dist'), fun.aggregate=mean, fill=NA)
    #     matplot(log1p(pk.dist[,.SD,.SDcols=patterns('^dist')]), log1p(pk.dist[,.SD,.SDcols=patterns('^tot')]), type='b', col=rep(1:nrow(variables), each=length(land_grid_list)), lty=1, pch=seq(length(land_grid_list)),
                #             ylab='peak intensity (E+I+C)', xlab='peak distance', main='Peak intensity vs. distance by timestep')

    inf.cells1 <- unique(wv.speed[is.na(edge),tot.inf.cells := length(tot.inf), by=.(v,l,r,time)][,.(v,l,r,time,tot.inf.cells)])
    inf.cells <- dcast(inf.cells1, time ~ v+l, value.var='tot.inf.cells', fun.aggregate=mean, fill=NA)
    matplot(log1p(inf.cells[,2:ncol(inf.cells)]), type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
            ylab='log(1+tot. infected cells)', xlab='week', main='Infected cells over time')

    #     inf.indiv1 <- unique(wv.speed[,tot.inf.indiv := sum(tot.inf), by=.(v,l,r,time)][,.(v,l,r,time,tot.inf.indiv)])
    #     inf.indiv <- dcast(inf.indiv1, time ~ v+l, value.var='tot.inf.indiv', fun.aggregate=mean, fill=NA)
    #     matplot(log1p(inf.indiv[,2:ncol(inf.indiv)]), type='b', lty=1, col=rep(1:nrow(variables), each=length(land_grid_list)), pch=seq(length(land_grid_list)),
                #             ylab='log(1+tot. infected individuals E+I+C)', xlab='week', main='Infected individuals over time')

    par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend("bottom", legend=c('vars:', 1:nrow(variables), 'land:', seq(length(land_grid_list))),
        lty=c(NA, rep(1,nrow(variables)), NA, rep(NA,length(land_grid_list))),
        col=c(NA, seq(nrow(variables)), NA, rep(1, length(land_grid_list))),
        pch=c(rep(NA, 2+nrow(variables)), seq(length(land_grid_list))),
        horiz=TRUE, bty='n', cex=1.3, lwd=2)

    dev.off()

    par(par.og)
    return(tm.mat.edge)
}
