
oos.tiles <- function(tiledat, tile_tab, variables){
    # the rest of the tiles' data for predictions
    tiledat <- fread(tiledat)
    setnames(tiledat, c('index', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd', 'x', 'y', 'state', 'gamma.shape', 'gamma.scale',
                        'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rd1', 'rd2','rd3','prcp', 'tmin', 'tmax', 'drt', 'contag', 'aggindex',
                        'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean', 'mvmt.cv'))
    setDT(variables)
    variables[,v := 1:.N]
    tvar <- CJ(l = unique(tiledat[,index]), var = variables[,v])
    tiledat <- tvar[tiledat, on=.(l=index)][variables, on=.(var=v)]
    tiledat[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('gamma.shape', 'gamma.scale', 'moranI', 'tc', 'mast', 'rgd', 'rd1','rd2','rd3','prcp','drt','contag','aggindex','entropy','simpindx')]
    tiledat[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact', 'variant', 'density')]
    tiledat <- tiledat[tile_tab, on=.(l=ID)]
    setnames(tiledat, c('wb_inside','wb_edge'), c('inside','edge'))
    return(tiledat)
}


oos.preds <- function(predmods, lasslist2, tiledat){
    resp.var <- as.character(predmods$fix[2])
    lambda.mins <- lasslist2[rsp.var == resp.var,lambda]
#     if(resp.var != lasslist2[,rsp.var]) stop('response variables are not matching with lambdas')
#     if (resp.var == 'exc') {
#         tiledat <- tiledat[CJ(l=unique(tiledat[,l]), tm=1:78), on=.(l), allow.cartesian=TRUE]
#     }
    tiledat.trim <- tiledat[,.(l, contact, variant, density, nnd_med, nnd_range, nnd_mean, nnd_sd, moranI, gearyC, tc, mast, rgd, rds, prcp, tmin, tmax, drt, contag, aggindex, entropy, simpindx, nnd_cv, mvmt.mean, mvmt.cv)]
    tiledat.trim[,l := as.factor(l)]
    tiledat.trim[,as.character(resp.var) := 0]
    pred <- predict(predmods, newdata=tiledat.trim)
    coefs.vals <- coef(predmods)
    coefs.used <- names(coefs.vals)[which(coefs.vals != 0)]
    pred.table <- cbind(tiledat, resp.variable=resp.var, coefs.used=paste(coefs.used, collapse=' '), lambda.nom=lambda.mins, lambda.mod=predmods$lambda.max, pred)
    return(pred.table)
}

preds.compile <- function(preds, tiledat){
    # aggregate outputs from oos.preds
    preds <- rbindlist(preds, fill=TRUE, use.names=TRUE)
    # dcast to wide-ish format
    preds.out <- dcast(preds, l + var ~ resp.variable, value.var=c('pred'), fun.aggregate = mean)
    # join with tiledat
    setnames(preds.out, names(preds.out)[-c(1,2)], paste0(names(preds.out[,-c(1,2)]), '_pred'))
    preds.out <- preds.out[tiledat, on=.(l, var)]
    return(preds.out)
}




maps.plot <- function(preds.table, contus_dir, variables){
    # national maps of response variables
    # not looped by targets
    library(sf)
    library(raster)
    library(terra)

    # get the contiguous us state border shapefile
    contus <- vect(contus_dir)
    contus.transform <- project(contus, 'epsg:5070')

    preds.table <- terra::crop(preds.table, aggregate(contus.transform, dissolve=TRUE))
    # create bounding box of contiguous US
    bounds <- ext(contus.transform)

    ## unified gradient across non-landscape parameters
    # loop through response variables
#     lapply(c('est','inf.spd','sounder.weeks','prop.infd','max.inc'), function(pred.var){
#         # get the values predicted by the models to set the ranges necessary for the map colors
#         pred.var.cols <- names(preds.table)[c(1:3, grep(paste0('^',pred.var,'_[0-9]'), names(preds.table)))]
#         pred.var.vals <- data.table(values(preds.table[,pred.var.cols]))
#         if (pred.var == 'sounder.weeks'){
#             pred.var.vals[, names(.SD) := lapply(.SD, log10), .SDcols=4:ncol(pred.var.vals)]
#             quants <- quantile(unlist(pred.var.vals[inside==1,-c(1:3)]), probs=c(0, 1))
#             bias <- 1#.2
#         } else if (pred.var == 'max.inc'){
#             pred.var.vals[, names(.SD) := lapply(.SD, log10), .SDcols=4:ncol(pred.var.vals)]
#             quants <- quantile(unlist(pred.var.vals[inside==1,-c(1:3)]), probs=c(0, 1))
#             bias <- 1
#         } else {
#             quants <- quantile(unlist(pred.var.vals[inside==1,-c(1:3)]), probs=c(0, 1))
#             bias <- 1
#         }
#         max.val <- quants[2]
#         min.val <- quants[1]
#         # scaled response variables from 0 to 1
#         pixels <- pred.var.vals[, lapply(.SD, function(x) max(0,(x-min.val)/(max.val-min.val))), .SDcols=patterns(paste0('^',pred.var)), by=ID]
# #         setnames(pixels, c('ID', names(pred.var.vals[,-c(1:3)])))
#         # define color gradient
#         color.grads <- colorRampPalette(c('darkgreen','yellowgreen','orange', 'darkred'), bias=bias)
#         # match colors to response variable
#         pixrank <- as.factor(as.numeric(cut(unlist(pixels[,.SD,.SDcols=patterns(paste0('^',pred.var))]), 255)))
#         pixrank <- data.table(matrix(pixrank, nrow=nrow(pixels)))
#         pixrank[,names(.SD) := lapply(.SD, as.numeric)]
#         setnames(pixrank, paste0(names(pixels[,-'ID']), '_r'))
#         pixrank <- cbind(pixels[,1], pixrank)
#         setDT(variables)
#         variables[,v:=1:.N]
#         # plot
#         png(paste0('./Output/figures/', pred.var, '_lm_tile_map.png'), width=2000, height=1000)
#         par(oma=c(0.2,0.2,3,0.2))
#         layout(matrix(c(seq(1,8),rep(9,4)), nrow=3, byrow=TRUE), widths=c(1,1,1,1), heights=c(1,1,0.3))
#         lapply(c(3,1,4,2,7,5,8,6), function(x){
#             variable.line <- unlist(variables[v==x, 1:3])
#             preds.tablei <- preds.table[,paste0(pred.var, '_', x)]
#             pixelsi <- pixrank[,.SD, .SDcols=paste0(pred.var, '_', x, '_r')]
#             terra::plot(preds.tablei, col=color.grads(255)[unlist(pixelsi)], pch=1, cex=1.3, ext=bounds, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n', border=NA)
#             mtext(paste('contact', variable.line[1], '| variant', variable.line[2], '| density', variable.line[3]), 3, -1.5, cex=1.5)
#             plot(contus.transform, col=NA, lwd=1.5, alpha=0, fill=NA, border='black', add=TRUE)
#         })
#         # add a legend
#         mtext(resp.translate(pred.var), 3, 0, cex=2, outer=TRUE)
#         legendimg <- as.raster(matrix(color.grads(255), nrow=1))
#         plot(c(0,1), c(0,1), type='n', axes=F, xlab = '', ylab='', main=resp.translate(pred.var), cex.main=1.5)
#         min.val2 <- floor(log10(abs(min.val)))
#         if (pred.var %in% c('max.inc','sounder.weeks')){
#             mtext(round((seq(min.val, max.val, length.out=5)), 1), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
#         } else {
#             mtext(round(seq(min.val, max.val, length.out=5), -min.val2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
#         }
#         rasterImage(legendimg, 0, 0, 1, 1)
#         dev.off()
#     })

    # loop through response variables
#     lapply(c('est','inf.spd','sounder.weeks','prop.infd','max.inc'), function(pred.var){
    lapply(c('max.inc'), function(pred.var){
        # get the values predicted by the models to set the ranges necessary for the map colors
        # plot
#         lapply(c(4), function(x){
        lapply(c(3,1,4,2,7,5,8,6), function(x){
            png(paste0('./Output/figures/', pred.var, '_lm_tile_map_',x,'.png'), width=1400, height=1000)
            par(oma=c(0.2,0.2,3,0.2))
    #         layout(matrix(c(seq(1,8),rep(9,4)), nrow=3, byrow=TRUE), widths=c(1,1,1,1), heights=c(1,1,0.3))
    #         layout(matrix(seq(1,16), nrow=4, byrow=FALSE), widths=c(1,1,1,1), heights=c(0.8,0.2,0.8,0.2))
            layout(matrix(seq(1,2), nrow=2, byrow=FALSE), widths=c(1), heights=c(0.8,0.15))
            pred.var.cols <- names(preds.table)[c(1:3, grep(paste0('^',pred.var,'_', x), names(preds.table)))]
            pred.var.vals <- data.table(values(preds.table[,pred.var.cols]))
            if (pred.var == 'sounder.weeks'){
                pred.var.vals[, names(.SD) := lapply(.SD, log), .SDcols=4:ncol(pred.var.vals)]
                quants <- quantile(unlist(pred.var.vals[inside==1,-c(1:3)]), probs=c(0, 1))
                bias <- 0.7
            } else if (pred.var == 'max.inc'){
                pred.var.vals[, names(.SD) := lapply(.SD, log10), .SDcols=4:ncol(pred.var.vals)]
                quants <- quantile(unlist(pred.var.vals[inside==1,-c(1:3)]), probs=c(0, 1))
                bias <- 0.7
            } else {
                quants <- quantile(unlist(pred.var.vals[inside==1,-c(1:3)]), probs=c(0, 1))
                bias <- 1
            }
            max.val <- quants[2]
            min.val <- quants[1]
            # scaled response variables from 0 to 1
            pixels <- pred.var.vals[, pmin(pmax(0,(.SD -min.val)/(max.val-min.val)), 1), .SDcols=patterns(paste0('^',pred.var)), by=ID]
            setnames(pixels, c('ID', names(pred.var.vals[,-c(1:3)])))
            # define color gradient
            color.grads <- colorRampPalette(c('darkgreen','yellowgreen','orange', 'darkred'), bias=bias)
            # match colors to response variable
            pixrank <- as.factor(as.numeric(cut(unlist(pixels[,.SD,.SDcols=patterns(paste0('^',pred.var))]), 255)))
            pixrank <- data.table(matrix(pixrank, nrow=nrow(pixels)))
            pixrank[,names(.SD) := lapply(.SD, as.numeric)]
            setnames(pixrank, paste0(names(pixels[,-'ID']), '_r'))
            pixrank <- cbind(pixels[,1], pixrank)
            setDT(variables)
            variables[,v:=1:.N]

            variable.line <- unlist(variables[v==x, 1:3])
            preds.tablei <- preds.table[,paste0(pred.var, '_', x)]
            pixelsi <- pixrank[,.SD, .SDcols=paste0(pred.var, '_', x, '_r')]
            terra::plot(preds.tablei, col=color.grads(255)[unlist(pixelsi)], pch=1, cex=1.3, ext=bounds, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n', border=NA, axes=FALSE)
            mtext(paste('contact', variable.line[1], '| variant', variable.line[2], '| density', variable.line[3]), 3, -1.5, cex=1.5)
            plot(contus.transform, col=NA, lwd=1.5, alpha=0, fill=NA, border='black', add=TRUE)

            # add a legend
            mtext(resp.translate(pred.var), 3, 0, cex=2, outer=TRUE)
            legendimg <- as.raster(matrix(color.grads(255), nrow=1))
            plot(c(0,1), c(0,1), type='n', axes=F, xlab = '', ylab='', main=resp.translate(pred.var), cex.main=1.5)
            min.val2 <- floor(log10(abs(min.val)))
            if (pred.var %in% c('max.inc','sounder.weeks')){
                mtext(round((seq(min.val, max.val, length.out=5)), 2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
            } else {
                mtext(round(seq(min.val, max.val, length.out=5), 2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
            }
            rasterImage(legendimg, 0, 0, 1, 1)
            dev.off()
        })
    })
}


## heatmaps
heatmap <- function(mod, tiledats, model.names){
    # looped by targets inputs on mod and lambda.mins
    # use the 2 most important landscape attributes, vary them for the axes, and for each combination of contact, variant, and density make a colored heatmap of response value predictions
    # mark points where actual landscapes were
    # use averages of other landscape values for model inputs

    in.mod <- mod[[1]]

    lambda.mins <- model.names[[1]][,lambda]
    # get the name of the response variable
    rsp.var <- model.names[[1]][,rsp.var]
    # pull in model coefficients
    coefs <- coef(in.mod)
    # keep only the landscape-associated values
    coefs <- coefs[(names(coefs) %in% c('(Intercept)', 'as.factor(contact)Lo', 'as.factor(density)5', 'as.factor(variant)Pol')==FALSE)]
    # keep the landscape variables from the tile data and average them to make an average for the background conditions of the top two parameters
    landval.avgs <- apply(tiledats[,c(3:6,12:28)], 2, mean) ## should use names here, to be sure
    # grab the two coefficients with the greatest effect sizes for the plot axes
    # these leave negative coefficients as negative (just use abs for ordering)
    coefs.sel <- coefs[order(abs(coefs), decreasing=TRUE)][1:2]
    if (all(names(coefs.sel) %in% c('moranI','gearyC'))) {
        coefs.sel <- coefs[order(abs(coefs), decreasing=TRUE)][c(1,3)]
    }
    print(coefs.sel)
    # grab coefficient names
    cnames <- names(coefs.sel)
    cf1nm <- cnames[1]
    cf2nm <- cnames[2]
    # get ranges of coefficients for each
    coef1.range <- range(tiledats[,..cf1nm])
    coef2.range <- range(tiledats[,..cf2nm])
    # create a table with all combinations of the active coefficients and the non-landscape variables
    ht.table <- CJ(
        cf1 = seq(coef1.range[1], coef1.range[2], length.out=25)
        , cf2 = seq(coef2.range[1], coef2.range[2], length.out=25)
        , density = c(1.5,5)
        , contact = c('Hi','Lo')
        , variant = c('DR','Pol'))
    # filter out landscape variables that are included in the plot axes
    landval.avgs <- landval.avgs[names(landval.avgs) %in% names(coefs.sel) == FALSE]
    setnames(ht.table, c('cf1','cf2'), c(cnames[1], cnames[2]))
    # connect coefficient combination with average landscape variables
    ht.table <- cbind(ht.table, as.data.table(t(landval.avgs)))
    # set population/epidemiology variables as factors
    ht.table[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('density','contact','variant')]
    # predict the response for each variable value combination
    print(rsp.var)
    ht.table <- ht.table[,.SD, .SDcols=c('contact', 'variant', 'density', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd',
                                         'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rd1','rd2','rd3', 'prcp', 'tmin',
                                         'tmax', 'drt', 'contag', 'aggindex', 'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean', 'mvmt.cv')]

    ht.table[,as.character(rsp.var):=0]
    ht.table[,l := as.factor(0)]
    ht.table[, as.character(rsp.var) := predict(in.mod, newdata=.SD)]
    # make table of variaable combinations to feed into mapply function below (to generate 8 heatmaps)
    quick.vars <- unique(ht.table[,.(density, contact, variant)])

    # define color ramp
    color.grads <- colorRampPalette(c('darkgreen','yellowgreen','orange', 'darkred'), bias=1)
    # match the color gradient with response variable values
    if (rsp.var %in% c('max.inc','sounder.weeks')){
        ht.table[,rank := as.factor( as.numeric( cut(log10(unlist(.SD)), 255))), .SDcols=c(rsp.var)]
    } else {
        ht.table[,rank := as.factor( as.numeric( cut(unlist(.SD), 255))), .SDcols=c(rsp.var)]
    }
    ht.table[,colr := color.grads(255)[as.numeric(as.character(rank))]]
    # plot
    png(paste0('./Output/figures/', rsp.var, '_heat_.png'), width=1700, height=1050)
    par(mfrow=c(2,4), oma=c(1,1,4,1))
    layout(matrix(c(1:9,9,9,9), nrow=3, ncol=4, byrow=TRUE), widths=c(1,1,1,1), heights=c(1,1,0.3))
    mapply(function(dens, cont, varnt){
        rows <- ht.table[,which(density == dens & contact == cont & variant == varnt)]
        ht.sub <- ht.table[rows, ]
        plot(unlist(ht.sub[,..cf1nm]), unlist(ht.sub[,..cf2nm]), col=ht.sub[,colr], pch=15, cex=4,ann=FALSE, pty='s')
        mtext(paste('density:', dens,'| contact:', cont, '| variant:', varnt), 3, 1.2, cex=1.4)
        mtext(exp.translate(cf1nm), 1, 2)
        mtext(exp.translate(cf2nm), 2, 2)
        tile.pts <- tiledats[variant == varnt & contact == cont & density == dens,.SD, .SDcols=c(cf1nm, cf2nm)]
        points(tile.pts)
    }, dens=quick.vars[,density], cont=quick.vars[,contact], var=quick.vars[,variant])
    mtext(resp.translate(rsp.var), 3, 2, outer=TRUE, cex=2)
    # create gradient legend object
    legendimg <- as.raster(matrix(color.grads(255), nrow=1))
    plot(c(0,1), c(0,1), type='n', axes=F, xlab = '', ylab='', main=resp.translate(rsp.var), cex.main=1.5)
    min.val2 <- floor(log10(abs(min(ht.table[,..rsp.var]))))-3
    if (rsp.var %in% c('max.inc', 'sounder.weeks')){
        mtext(round(seq(range(log10(ht.table[,..rsp.var]))[1], range(log1p(ht.table[,..rsp.var]))[2], length.out=5), -min.val2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
    } else {
        mtext(round(seq(range(ht.table[,..rsp.var])[1], range(ht.table[,..rsp.var])[2], length.out=5), -min.val2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
    }
    rasterImage(legendimg, 0, 0, 1, 1)
    dev.off()
    return(NA)
}



resp.translate <- function(rv){
    # translates response variable codes into text for labeling figures
    print(rv)
    if (rv == 'est')            rv.out <- 'Epidemic Establishment Proportion'
    if (rv == 'inf.spd')        rv.out <- 'Epidemic Wave Max Speed (km/wk)'
    if (rv == 'sounder.weeks')  rv.out <- 'log10(Epidemic Intensity (Sounder-Weeks))'
    if (rv == 'prop.infd')      rv.out <- 'Proportion of Cells Infected'
    if (rv == 'max.inc')        rv.out <- 'log10(Maximum Incidence)'
    if (rv == 'exc')            rv.out <- 'Proportion of Epidemics Escaped by 52 Weeks'
    if (rv == 'tm.esc')         rv.out <- 'Proportion of Epidemics Escaped by 52 Weeks'
    return(rv.out)
}

exp.translate <- function(exp.var){
    # translates explanatory variable codes into text for labeling figures
    if (exp.var == 'nnd_med')   ev.out <- 'Median Nearest Neighbor Distance'
    if (exp.var == 'nnd_range') ev.out <- 'Range Nearest Neighbor Distance'
    if (exp.var == 'nnd_mean')  ev.out <- 'Mean Nearest Neighbor Distance'
    if (exp.var == 'nnd_cv')    ev.out <- 'Nearest Neighbor Coefficient of Variation'
    if (exp.var == 'nnd_sd')    ev.out <- 'Nearest Neighbor Standard Deviation'
    if (exp.var == 'mvmt.mean') ev.out <- 'Mean Sounder Movement'
    if (exp.var == 'mvmt.cv')   ev.out <- 'Sounder Movement Coefficient of Variation'
    if (exp.var == 'moranI')    ev.out <- 'Landscape Preference Autocorrelation (Moran\'s I)'
    if (exp.var == 'gearyC')    ev.out <- 'Landscape Preference Autocorrelation (Geary\'s C)'
    if (exp.var == 'tc')        ev.out <- 'Tree Cover'
    if (exp.var == 'mast')      ev.out <- 'No. Masting Species'
    if (exp.var == 'rgd')       ev.out <- 'Landscape Ruggedness'
    if (exp.var == 'rds')       ev.out <- 'Roads Index'
    if (exp.var == 'dayl')      ev.out <- 'Daylight'
    if (exp.var == 'prcp')      ev.out <- 'Precipitation'
    if (exp.var == 'tmax')      ev.out <- 'Mean Annual Maximum Daily Temperature'
    if (exp.var == 'tmin')      ev.out <- 'Mean Annual Minimum Daily Temperature'
    if (exp.var == 'simpindx')  ev.out <- 'Landscape Cover Simpson Diversity Index'
    if (exp.var == 'drt')       ev.out <- 'Drought'
    if (exp.var == 'contag')    ev.out <- 'Contagion'
    if (exp.var == 'aggindex')  ev.out <- 'Aggregation index'
    if (exp.var == 'entropy')   ev.out <- 'Entropy'
    return(ev.out)
}
