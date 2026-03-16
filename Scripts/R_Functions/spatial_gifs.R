## animation of spatial populations over time

spatial.gif <- function(parameters, detections, incidence, land_grid_list, var.in=NA, rep.in=NA, lnd.in=NA){


    grid.centers.out <- rbindlist(lapply(seq(length(land_grid_list)), function(i){
        grid.key <- as.data.table(land_grid_list[[i]][[2]])
        setnames(grid.key, c('cell','tlX','tlY','trX','trY','ctX','ctY','rasval')[1:ncol(grid.key)]) ## see Make_Grid.R
        grid.centers <- cbind(grid.key[,.(cell, ctX, ctY, rasval)], land=i)
        return(grid.centers)
    }))

    radius = parameters['Rad']
    #     unq.combos <- unique(unq.eic[max.time > 10,.(var, land, rep)])
    if (radius != 0 & parameters[['detectday']] < parameters[['thyme']]){
        ## plotting incidence with zone of control over time for a bunch of plots
        # detection locations
        unq.detections <- unique(detections[,.(var,land,rep,time,loc,max.time,detected)])

        unq.incidence <- unique(incidence[,.(var,land,rep,timestep,loc,max.time)])
        unq.incidence[,is.inf := 1]
        unq.incidence[,loc.min := min(timestep), by=.(var, land, rep, loc)]
        unq.incidence[,loc.max := max(timestep), by=.(var, land, rep, loc)]
        unq.join <- CJ(timestep=seq(parameters[['thyme']]),
                    var=unique(unq.incidence[,var]),
                    land=unique(unq.incidence[,land]),
                    rep=unique(unq.incidence[,rep]),
                    loc=unique(unq.incidence[,loc]))
        unq.join <- unq.join[!unq.incidence, on=.(var, land, rep, timestep, loc)]
        unq.join[, is.inf := 0]
        unq.join <- unq.join[unique(unq.incidence[,.(var, land, rep, loc, loc.min, loc.max, max.time)]), on=.(var, land, rep, loc)]
        unq.join <- unq.join[timestep >= loc.min,]
        #         unq.join <- unq.join[timestep <= loc.max & timestep >= loc.min,]
        unq.incidence <- rbind(unq.incidence, unq.join)

        unq.incidence[,loc.min := NULL]
        unq.incidence[,loc.max := NULL]
        setorder(unq.incidence, var, land, rep, timestep, loc)
        ## gets all runs than last at least 50 weeks, infects at least 5 cells during the course of the run, and then takes the first and last rep of each var/land combo to use for gif example
        unq.combos <- unique(unq.incidence[max.time > 50,sum(is.inf), by=.(var, land, rep, timestep)][,max(V1), by=.(var, land, rep)][V1 > 5, .(var,land,rep)])[, .SD[.N], by=.(var, land)]

        #input var and rep values determine which variable combinations (1:8) and which reps (1:100) to create gifs for. lands are trickier...
        if (!is.na(var.in)){
            unq.combos <- unq.combos[var %in% var.in,]
        }
        if (!is.na(rep.in)){
            unq.combos <- unq.combos[rep %in% rep.in,]
        }
        ## land numbers may not correspond with the name of the input tif file...
        if (!is.na(lnd.in)){
            unq.combos <- unq.combos[land %in% lnd.in,]
        }

        mcmapply(function(i, j, k){
            sub.incidence <- unique(unq.incidence[var == i & land == j & rep == k & loc != 0,.(timestep, loc, is.inf)])[grid.centers.out[land==j,], on=.(loc = cell), nomatch=NULL] # gets rid of code column
            if(nrow(unq.detections) > 0){ sub.detections <- unique(unq.detections[var == i & land == j & rep == k & detected != 0,.(time, loc)])[grid.centers.out[land==j,], on=.(loc = cell), nomatch=NULL] } # gets rid of code column
            if(nrow(allzones) >0){sub.zones <- allzones[var == i & land == j & rep == k & loc != 0,][grid.centers.out[land==j,], on=.(loc=cell), nomatch=NULL]}
            sub.solocs <- solocs.all[v == i & l == j & r == k,]
            paneldim <- ceiling(sqrt(max(sub.incidence[,timestep])))
            png(paste0('./Output/figures/sptmplot_vars_', i,'_land_', j, '_rep_', k,'.png'), width=2000, height=2000)
            par(mfrow = c(paneldim, paneldim), oma=c(0,0,0,0), mar=c(0,0,0,0))
            lapply(seq(max(sub.incidence[,timestep])), function(x){
                plot(y ~ x, data=sub.solocs[time == x,], pch='.', col='gray', xlim=c(0,100), ylim=c(0,100))
                if(nrow(allzones) > 0){points(ctY ~ ctX, data=sub.zones[timestep <= x,], pch='.', col='yellow')}#, cex=1.2, xlim=c(0,100), ylim=c(0,100))
        #                 points(y ~ x, data=sub.solocs[timestep == x,], pch=1, col='gray')
        #                 points(ctY ~ ctX, data=sub.incidence[timestep <= x & is.inf == 0,], pch='.', col='red')
                points(ctY ~ ctX, data=sub.incidence[timestep == x & is.inf == 1,], pch=3, col='red')
                if (x > parameters['detectday'] & exists('sub.detections')) points(ctY ~ ctX, data=sub.detections[time <= x,], pch=2, col='blue')
            })
            dev.off()

            library(animation)
            plot.step <- function(x){
                panels <- layout(matrix(c(1, 1, 1, 2, 3, 4), nrow=3, ncol=2), widths=c(3, 1), heights=c(1, 1, 1), respect=FALSE)
                plot(y ~ x, data=sub.solocs[time == x,], pch='.', col='gray', xlim = c(0, 100), ylim = c(0, 100), main=paste('week', x), cex=log1p(nlive))
                if(nrow(allzones) > 0) points(ctY ~ ctX, data=sub.zones[time <= x,], pch=3, col='yellow')
                points(ctY ~ ctX, data=sub.incidence[timestep == x & is.inf == 1,], pch=3, col='red')
                if (x > parameters['detectday'] & exists('sub.detections')) points(ctY ~ ctX, data=sub.detections[time <= x,], pch=2, col='blue')
                plot(allnlive ~ time, data=sub.solocs[, allnlive := sum(nlive), by=time][order(time)], main='Live pigs', type='l')
                points(allnlive ~ time, data=sub.solocs[time==x,], pch=1, cex=2)
                plot(allnlive ~ timestep, data=unique(sub.incidence[, allnlive := sum(is.inf), by=timestep][,.(timestep, allnlive)])[order(timestep)], main='Infected cells', type='l', col='red')
                points(allnlive ~ timestep, data=sub.incidence[timestep==x,], pch=3, col=2, cex=2)
                if (exists('sub.detections')){
                    temp.detect <- unique(sub.detections[, allnlive := length(unique(loc)), by=time][,.(time, allnlive)])[order(time)][,cs:=cumsum(allnlive)]
                    plot(cs ~ time, data=temp.detect, main='Detections', type='l', col='blue')
                    points(cs ~ time, data=temp.detect[time==x, ], pch=2, col='blue', cex=2)
                }
            }
            out.gif <- paste0('testgif_vars_', i, '_land_', j, '_rep_', k, '.gif')
            saveGIF({
                lapply(seq(min(sub.incidence[,timestep]), max(sub.incidence[,timestep])), plot.step)
            }, movie.name = out.gif, ani.width=850, ani.height=600, interval=0.1, imgdir='./Output/figures') ## imgdir argument doesn't work, so there is a line outside the function to move all gives to the right place
        }, i=unq.combos[,var], j=unq.combos[,land], k=unq.combos[,rep])

#         moves gifs to the proper folder
        lapply(list.files(pattern='*.gif'), function(x) file.rename(x, paste0('./Output/figures/',x)))
    }
}
