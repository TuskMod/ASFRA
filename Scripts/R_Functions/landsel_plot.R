## create map image of selected lands

### Could add background of land cover types and/or regions with feral pig activity?


sel.lands <- function(landlist = NULL){
    library(sf)
    library(raster)
    library(terra)

    contus <- st_read('./Input/contus_shp/') # for windows must be directory, might be something else for linux?
    sel.files <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)
    selpts.box <- bind(lapply(sel.files, function(tif){
        box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
        return(box)
    }))
    tile.potential <- list.files('./Landscape_Setup/Pipeline_SSF_Weekly/4_Output/indiv_plands', full.names=TRUE)
    pot.pts.box <- bind(lapply(tile.potential, function(tif){
        box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
        return(box)
    }))

    tile.one <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)[1]
    rast.crs <- crs(raster(tile.one))
    contus.transform <- st_transform(contus[1], rast.crs)

    bounds <- st_bbox(contus.transform)
    xrng <- bounds[c(1,3)]
    yrng <- bounds[c(2,4)]

    map.scale = 100
    ysz <- round(log(yrng[2] - yrng[1]) * map.scale)
    xsz <- round(log(xrng[2] - xrng[1]) * map.scale)

    browser()
    png('./Output/figures/tile_map.png', width=1600, height=1200)
    plot(pot.pts.box, col=NA, pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n')
    plot(do.call(merge, lapply(list.files('./Landscape_Setup/Pipeline_SSF_Weekly/4_Output/indiv_plands/', full.names=TRUE), raster)), add=TRUE)
    polys(contus.transform, col=0, lwd=2, alpha=1, fill=NA)
    plot(pot.pts.box, col=NA, border='gray', add=TRUE)
    plot(selpts.box, col='pink', border='red', fill='pink', lwd=2, alpha=0.7, add=TRUE)
    dev.off()
}



pred.lands <- function(pred.var, atype, var.set, tiledat){
    # prediction variable (est, max.dist, etc.)
    # analysis type (g/lmer, gbm)
    # vector of response variables
    # table of all tile data (averaged responses from training tiles)
    library(sf)
    library(raster)
    library(terra)

    contus <- st_read('./Input/contus_shp/') # for windows must be directory, might be something else for linux?
    sel.files <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)
    selpts.box <- bind(lapply(sel.files, function(tif){
        box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
        return(box)
    }))
    name <- tstrsplit(sel.files, '_')
    name <- name[length(name)]
    name <- unlist(tstrsplit(unlist(name), '\\.', keep=1))
    name <- as.numeric(name)
    selpts.box <- SpatialPolygonsDataFrame(selpts.box, data=data.table(l=name))

    tile.potential <- list.files('./Landscape_Setup/Pipeline_SSF_Weekly/4_Output/indiv_plands', full.names=TRUE)
    pot.pts.box <- bind(lapply(tile.potential, function(tif){
        name <- tstrsplit(tif, '_')
        name <- name[length(name)]
        name <- paste0('tile', unlist(tstrsplit(unlist(name), '\\.', keep=1)))
        box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
        names(box) <- name
        return(box)
    }))
    name <- tstrsplit(tile.potential, '_')
    name <- name[length(name)]
    name <- unlist(tstrsplit(unlist(name), '\\.', keep=1))
    name <- as.numeric(name)
    pot.pts.box <- SpatialPolygonsDataFrame(pot.pts.box, data=data.table(l=name))

    tile.one <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)[1]
    rast.crs <- crs(raster(tile.one))
    contus.transform <- st_transform(contus[1], rast.crs)

    bounds <- st_bbox(contus.transform)
    xrng <- bounds[c(1,3)]
    yrng <- bounds[c(2,4)]

    map.scale = 100
    ysz <- round(log(yrng[2] - yrng[1]) * map.scale)
    xsz <- round(log(xrng[2] - xrng[1]) * map.scale)


    # 1 tile for each varaible type
    png(paste0('./Output/figures/tile_map', pred.var, '.png'), width=1600, height=1200)

    par(mfrow=c(2,8))
    lapply(1:8, function(x){
    #     connect tiledat and prediction data to specific tiles
        pot.pts.box.sub <- merge(pot.pts.box, tiledat[v==x, .SD, .SDcols=c('l', pred.var)], by='l')
        selpts.box.sub <- merge(selpts.box, tiledat[v==x, .SD, .SDcols=c('l', pred.var)], by='l')
#         plot(selpts.box.sub, col=NA, pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n')
        plot(pot.pts.box.sub, col=NA, pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n')
        polys(contus.transform, col=0, lwd=2, alpha=1, fill=NA)
        plot(pot.pts.box.sub, col=NA, border='gray', add=TRUE)#, fill=rgb())
        plot(pot.pts.box.sub, border='gray', add=TRUE, col=rgb(unlist(as.data.table(pot.pts.box.sub@data)[,..pred.var]),0,0))
        plot(selpts.box.sub, col=rgb(unlist(as.data.table(selpts.box.sub@data)[,..pred.var]),0,0), border='red', fill=NA, lwd=2, alpha=0.7, add=TRUE)
    })
    dev.off()
}

