
VisualOutputs <- function(lvtable, variables, land_grid_list, parameters){
    ## NEED TO:
    ## -- add land image to background of gifs?
    ## -- get ldsel table from NND calc for attributes of landscapes


    # Generates output plots for examples
    library(data.table)
    library(parallel)
    library(akima)
    library(rgl)

    setDT(variables)

    core.dir <- './Output/'
#     core.dir <- 'C:/Users/eric.sodja/Documents/testrun_outputs/Output/'

    ## essentially, if something has changed, re-run the compilation processes in this if statement; otherwise, just read in the existing files (could be targets-ed, but this is simpler...)
    if('lvtable.csv' %in% list.files('./Output/')){
        lvtable.in <- fread(paste0(core.dir, 'lvtable.csv'))
    } else {
        fwrite(lvtable, paste0(core.dir, 'lvtable.csv'))
        lvtable.in <- data.table()
    }
    if (nrow(lvtable.in) != nrow(lvtable)){
        fwrite(lvtable, paste0(core.dir, 'lvtable.csv'))
        # make things into data.tables and name columns for sanity
        mapp.to.dt <- function(x, dt){
            out <- data.table(matrix(unlist(dt[,x]), ncol=length(dt[,x]), byrow=FALSE))
            setnames(out, names(dt[,x]))
        }
        print('Reading in timestep demographics...')
        tm.mat <- mcmapply(function(v, l, r){
            filenm <- paste0('tm.mat_r', r,'_l',l, '_v',v,'.csv')
            if (filenm %in% list.files(paste0(core.dir, 'tm.mat/'))){
                filenm <- paste0(paste0(core.dir, 'tm.mat/', filenm))
                tm.mat <- fread(filenm)
                setnames(tm.mat, unlist(lapply(strsplit(names(tm.mat), 'tm.mat.'), function(x) unlist(x)[2])))
                tm.mat <- tm.mat[!is.na(timestep),]
                return(tm.mat)
            }
        }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
        tm.mat <- rbindlist(lapply(1:ncol(tm.mat), mapp.to.dt, dt=tm.mat))
        fwrite(tm.mat, paste0(core.dir, 'tm.mat.agg.csv'))

        print('Reading in summary statistics...')
        summ.vals <- mcmapply(function(v, l, r){
    #         summ.vals <- as.data.table(out.list["summ.vals"])
            filenm <- paste0('summ.vals_r', r,'_l',l, '_v',v,'.csv')
            if (filenm %in% list.files(paste0(core.dir, 'summ.vals/'))){
    #             print(filenm)
                filenm <- paste0(paste0(core.dir, 'summ.vals/', filenm))
                summ.vals <- fread(filenm)
                setnames(summ.vals, unlist(lapply(strsplit(names(summ.vals), 'summ.vals.'), function(x) unlist(x)[2])))
                summ.vals <- summ.vals[!is.na(endtime)]
                return(summ.vals)
            }
        }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
        summ.vals <- rbindlist(lapply(1:ncol(summ.vals), mapp.to.dt, dt=summ.vals))
        fwrite(summ.vals, paste0(core.dir, 'summ.vals.agg.csv'))
        print('Reading in incidence measures...')
        incidence <- mcmapply(function(v, l, r){
            filenm <- paste0('incidence_r', r,'_l',l, '_v',v,'.csv')
            if (filenm %in% list.files(paste0(core.dir, 'incidence/'))){
                filenm <- paste0(paste0(core.dir, 'incidence/', filenm))
                incidence <- fread(filenm)
                setnames(incidence, unlist(lapply(strsplit(names(incidence), 'incidence.'), function(x) unlist(x)[2])))
                incidence[, max.time := max(timestep), by=.(var, rep, land)]
                return(incidence)
            }
        }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
        incidence <- rbindlist(lapply(1:ncol(incidence), mapp.to.dt, dt=incidence))
        fwrite(incidence, paste0(core.dir, 'incidence.agg.csv'))

        ## for later
    #     detections <- rbindlist(mcmapply(function(v, l, r){
    # #         detections <- as.data.table(out.list["detections"])
    #         detections <- fread(paste0(paste0(core.dir, 'detections/detections_r', r,'_l',l, '_v',v,'.csv'))
    #         setnames(detections, unlist(lapply(strsplit(names(detections), 'detections.'), function(x) unlist(x)[2])))
    #         detections[, max.time := 0]
    #         if(nrow(detections) > 0) {
    #             detections[,max.time := max(time), by=.(var, rep, land)]
    #             unq.det <- unique(detections[, .(var, land, rep, timestep, loc, max.time, code, detected)])
    #         }
    #         return(detections)
    #     }
    #         allzones <- as.data.table(out.list["allzone"])
    #         setnames(allzones, unlist(lapply(strsplit(names(allzones), 'allzone.'), function(x) unlist(x)[2])))

        print('Reading in sounderlocations...')
        wv.speed <- mcmapply(function(v, l, r){
            filenm <- paste0('solocs.all_r', r,'_l',l, '_v',v,'.csv')
            if (filenm %in% list.files(paste0(core.dir, 'solocs.all/'))){
                filenm <- paste0(paste0(core.dir, 'solocs.all/', filenm))
                solocs.all <- fread(filenm)
                setnames(solocs.all, unlist(lapply(strsplit(names(solocs.all), 'solocs.all.'), function(x) unlist(x)[2])))
                solocs.all[,nlive := S+E+I+R]
                wvspd <- wave_speed(solocs.all)
                return(wvspd)

            }
        }, v= lvtable[,vars], l= lvtable[,land], r=lvtable[,rep])
    #     solocs.all <- rbindlist(mclapply(1:ncol(solocs.all), mapp.to.dt, dt=solocs.all))
        wv.speed <- rbindlist(mclapply(1:ncol(wv.speed), mapp.to.dt, dt=wv.speed))
        fwrite(wv.speed, paste0(core.dir, 'wv.speed.agg.csv'))
    } else {
        tm.mat <- fread(paste0(core.dir, 'tm.mat.agg.csv'))
        summ.vals <- fread(paste0(core.dir, 'summ.vals.agg.csv'))
        incidence <- fread(paste0(core.dir, 'incidence.agg.csv'))
        wv.speed <- fread(paste0(core.dir, 'wv.speed.agg.csv'))
    }

    print('Reading in selected land attributes...')
    ldsel.nnd <- readRDS('./Landscape_Setup/NND_Lands/4_Output/ldsel.rds')
    setDT(ldsel.nnd)
    ldsel.nnd[,land:=NULL]
    setnames(ldsel.nnd, 'index', 'land')

    print('Generating state variable temporal plots')
    seircz.plot(tm.mat, variables)

    ## wavespeed stuff
    print('Plotting wave speed metrics')

    tm.mat.edge <- wave.plot(wv.speed, tm.mat, land_grid_list, variables)

    ## Proportion of simulations that ASF establishes
    print('Generating plots of establishment success...')
    unq.parms <- estab.plot(tm.mat, ldsel.nnd, variables)

    # takes anything with (1) 2 or more cases at a given time, OR (2) at least 10 cases total
    ## should be establishment criteria from madison's paper
    tm.mat[,max.inf := max(I), by=.(var, land, rep)]
    tm.mat[,tot.exp := sum(E), by=.(var, land, rep)]
    tm.mat[,est := 0]
    tm.mat[max.inf > 1 & tot.exp >= 10, est := 1]
    #     established <-  unique(tm.mat[,est.no := sum(est, na.rm=TRUE), by=.(var, land)][,est.pct := est.no/.N, by=.(var, land)][,.(var, land, est.no, est.pct)])
    established <- unique(tm.mat[,.(var, land, rep, max.inf, tot.exp, est)])
    established <- established[variables[,.(contact, variant, density, id)], on=.(var=id)]
    established <- established[ldsel.nnd, on=.(land=land), nomatch=NULL]
    established[,density := as.character(density)][density == '1.5', density := 'Lo'][density == '5', density := 'Hi']
    established[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('moranI', 'gearyC', 'disp')]
    established[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('est','contact','variant','land','density')]

    prcc.table <- established[,.(var, land, rep, contact, variant, density, nnd_med, nnd_range, disp, moranI, est, max.inf, tot.exp)]
    library(epiR)
    prcc.est <- epi.prcc(established[,.(contact, variant, density, nnd_med, disp, moranI, est)])
    prcc.maxinf <- epi.prcc(prcc.table[,c(4:10,12)])
    prcc.totexp <- epi.prcc(prcc.table[,c(4:10,13)])

#     estab.mod <- glmer(est ~ contact + variant + density + moranI + nnd_med + disp + (1|land), data=established, family=binomial, na.action=na.fail, control=glmerControl(autoscale=TRUE))
    var.mods.est <- lapply(unique(established[,var]), function(x){
        idx <- established[,which(var == ..x)]
        est.subset <- established[idx,]
        estab.mod <- glm(est ~ moranI + disp, data=est.subset, family=binomial)
#         estab.mod <- glmer(est ~ morani + nnd_med + disp + (1|land), data=est.subset, family=binomial, na.action=na.fail, control=glmercontrol(autoscale=true))
        return(estab.mod)
    } )

    morani.range <- range(established[,moranI])
    dispersal.range <- range(established[,disp])
    var.mods.pred <- CJ(morani = seq(morani.range[1], morani.range[2], length=30), dispersal=seq(dispersal.range[1], dispersal.range[2], length=30))
    val.estimate <- function(mod.no){
        mod.def <- var.mods.est[[mod.no]]
        coef <- coefficients(mod.def)
        int <- coef[1]
        mi <- coef[2]
        disp <- coef[3]
        var.mods.pred[, v1 := int + mi*morani + disp*dispersal]
        setnames(var.mods.pred, 'v1', paste0('mod',mod.no))
        return(var.mods.pred[,paste0('mod',..mod.no)])
    }

    val.est.tab.out <- as.data.table(lapply(1:8, val.estimate))
    est.preds <- cbind(var.mods.pred, val.est.tab.out)

    est.pts <- unique(established[,.(moranI, disp)])

    png('./Output/figures/heat_estab.png', width=1000, height=1000)
    par(mfrow=c(3,3))
    lapply(1:8, function(x){
        plot(dispersal ~ morani, data=est.preds, col=rgb(rescale(unlist(est.preds[,..x])), 0, 0), pch=15, cex=2, main = x, asp=1)
        points(disp ~ moranI, data=est.pts, col='yellow', pch=1)
    })
    dev.off()

    ## max observed incidence for var/land combos
    print('Generating heatmaps of incidence...')
    maxinc.heat(tm.mat.edge, unq.parms, ldsel.nnd)
    max.incidence <- tm.mat.edge[, max(incidenceEI, na.rm=TRUE), by=.(var, land)][unq.parms, on='var'][ldsel.nnd, on='land']
    setnames(max.incidence, 'V1','max.inc')
    max.incidence[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('moranI', 'gearyC', 'disp')]
    max.incidence[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact','variant','land')]
    prcc.maxinc <- epi.prcc(max.incidence[,.(contact, variant, density, nnd_med, disp, moranI, max.inc)])

    var.mods.inc <- lapply(unique(max.incidence[,var]), function(x){
        idx <- max.incidence[,which(var == ..x)]
        inc.subset <- max.incidence[idx,]
        inc.mod <- lm(max.inc ~ moranI + disp, data=inc.subset)
#         estab.mod <- glmer(est ~ morani + nnd_med + disp + (1|land), data=est.subset, family=binomial, na.action=na.fail, control=glmercontrol(autoscale=true))
        return(inc.mod)
    } )

    morani.range <- range(established[,moranI])
    dispersal.range <- range(established[,disp])
    var.mods.pred <- CJ(morani = seq(morani.range[1], morani.range[2], length=30), dispersal=seq(dispersal.range[1], dispersal.range[2], length=30))
    val.estimate <- function(mod.no){
        mod.def <- var.mods.inc[[mod.no]]
        coef <- coefficients(mod.def)
        int <- coef[1]
        mi <- coef[2]
        disp <- coef[3]
        var.mods.pred[, v1 := int + mi*morani + disp*dispersal]
        setnames(var.mods.pred, 'v1', paste0('mod',mod.no))
        return(var.mods.pred[,paste0('mod',..mod.no)])
    }

    val.inc.tab.out <- as.data.table(lapply(1:8, val.estimate))
    inc.preds <- cbind(var.mods.pred, val.est.tab.out)

    inc.pts <- unique(established[,.(moranI, disp)])

    png('./Output/figures/heat_maxincidence.png', width=1000, height=1000)
    par(mfrow=c(3,3))
    lapply(1:8, function(x){
        plot(dispersal ~ morani, data=inc.preds, col=rgb(rescale(unlist(inc.preds[,..x])), 0, 0), pch=15, cex=2, main = x, asp=1)
        points(disp ~ moranI, data=inc.pts, col='yellow', pch=1)
    })
    dev.off()
    browser()

    ## max.dist.wave
    print('Generating heatmaps of extent of infection...')
    wvdist.plot(wv.speed, unq.parms, ldsel.nnd, wvtime=60)
    wvtime <- 60
    wvdist.t60 <- unique(wv.speed[is.na(edge), maxtime := max(time), by=.(v,l)][time <=wvtime,][time == wvtime | time == maxtime, .(v,l,r,time, max)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    wvdist.t60[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('moranI','gearyC','disp')]
    wvdist.t60[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact','variant','l')]
    prcc.wvdist <- epi.prcc(wvdist.t60[,.(contact, variant, density, nnd_med, disp, moranI, max)])


    ## spd.wave
    print('Generating heatmaps of wavespeed')
    wvspd.plot(wv.speed, unq.parms, ldsel.nnd)
    wvspd.plot.moran(wv.speed, unq.parms, ldsel.nnd)

    wvspd.avgdistdiff <- unique(wv.speed[is.na(edge), mean(avg.dist.diff), by=.(v,l)])[unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(wvspd.avgdistdiff, 'V1','avg.spd')
    wvspd.avgdistdiff[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('disp','sigdisp','gamma.shape','gamma.rate','moranI','gearyC')]
    wvspd.avgdistdiff[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact','variant','l')]
    prcc.wvspd <- epi.prcc(wvspd.avgdistdiff[,.(contact, variant, density, nnd_med, disp, moranI, avg.spd)])

#     summary(lm(avg.spd ~ as.factor(contact) + as.factor(variant) + density + nnd_med + disp + moranI, data=wvspd.avgdistdiff))
#     summary(lm(avg.spd ~ (as.factor(contact) + as.factor(variant) + density + nnd_med + disp + moranI)^2, data=wvspd.avgdistdiff))
#     summary(lm(avg.spd ~ nnd_med + disp + moranI, data=wvspd.avgdistdiff))
#     summary(lm(avg.spd ~ (density + nnd_med + disp + moranI)^2, data=wvspd.avgdistdiff))
#     summary(lm(avg.spd ~ as.factor(contact) + as.factor(variant) + density, data=wvspd.avgdistdiff))
#     summary(lm(avg.spd ~ (as.factor(contact) + as.factor(variant) + density)^2, data=wvspd.avgdistdiff))

    ## inf.cells
    print('Generating heatmaps of infected cells')
    infcell.plot(wv.speed, unq.parms, ldsel.nnd)
    max.inf.cells <- wv.speed[is.na(edge), max(tot.inf.cells), by=.(v,l,r)][unq.parms, on=.(v=var)][ldsel.nnd, on=.(l=land)]
    setnames(max.inf.cells, 'V1', 'max.inf.cells')
    max.inf.cells[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('disp','sigdisp','gamma.shape','gamma.rate','moranI','gearyC')]
    max.inf.cells[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact','variant','l')]
    prcc.infcell <- epi.prcc(max.inf.cells[,.(contact, variant, density, nnd_med, disp, moranI, max.inf.cells)])

    barplot.prcc <- function(prcc.table, title, color){
        bp <- barplot(prcc.table$est, names.arg=prcc.table$var, ylim=c(-1,1), main=title, ylab='Effect size')
        abline(h=0, col='gray')
        bp <- barplot(prcc.table$est, names.arg=prcc.table$var, ylim=c(-1,1), main=title, add=TRUE, col=color)
        segments(bp, prcc.table[,'lower'], bp, prcc.table[,'upper'], lend=0)
        ew <- (bp[2,1]-bp[1,1])/4
        segments(bp-ew, prcc.table[,'lower'], bp+ew, prcc.table[,'lower'], lend=0)
        segments(bp-ew, prcc.table[,'upper'], bp+ew, prcc.table[,'upper'], lend=0)
        text(bp, y=-0.9, paste0('p=',round(prcc.table$p.value, 2)))
    }
    png('./Output/figures/prccplot.png', width=1000, height=1000)
    par(mfrow=c(2,3))
    barplot.prcc(prcc.est, title = 'Proportion of established epidemics', color=2)
    barplot.prcc(prcc.maxinc, title = 'Max incidence', color=3)
    barplot.prcc(prcc.wvdist, title = 'Max wave distance by 60 weeks', color=4)
    barplot.prcc(prcc.infcell, title = 'Max infected cells', color=5)
    barplot.prcc(prcc.wvspd, title = 'Max wave speed', color=6)
    dev.off()
    ## create map image of selected lands

    print('Generating national map of selected lands...')
    sel.lands(ldsel.nnd)
    browser()

    ## landscape attributes
    # land_grid_list
    # 1: number of cells
    # 2: table:
    #     1: cell id
    #     2: dim edge low
    #     3: dim edge high
    #     4: dim edge low
    #     5: dim edge high
    #     6: center point x?
    #     7: center point y?
    #     8: cell preference value
#     lapply(seq(length(land_grid_list)), function(x){
#         in.rast <- rast(land_grid_list[[x]][[3]], type='xyz')
#         corr <- autocor(in.rast, global=TRUE)
#         return(corr)
#     })
#     library('magick')
#     contour_stack <- function(pattern, outname, dir='./Output/figures/'){
#         pattern <- paste0('^', pattern)
#         imnames <- list.files(dir, pattern=pattern, include.dirs=TRUE, full.names=TRUE)
#         rowA <- image_append(c(image_read(imnames[1]), image_read(imnames[2])))
#         rowB <- image_append(c(image_read(imnames[3]), image_read(imnames[4])))
#         rowC <- image_append(c(image_read(imnames[5]), image_read(imnames[6])))
#         rowD <- image_append(c(image_read(imnames[7]), image_read(imnames[8])))
#         imgout <- image_append(c(rowA, rowB, rowC, rowD), stack=TRUE)
#         image_write(imgout, paste0('./Output/figures/',outname,'.png'))
#     }
#
#     mapply(contour_stack,
#         pattern = c('land_parcombos_maxincd','land_parcombos_estab','land_parcombos_maxincd','land_parcombos_infcells','land_parcombos_wvspd','land_parcombos_maxdist'),
#         outname = c('maxincd','estab','maxincd','infcells','wvspd','maxdist'))
}


