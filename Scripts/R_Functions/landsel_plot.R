## create map image of selected lands

sel.lands <- function(ldsel.nnd){
    library(sf)
    library(raster)
    library(terra)
    contus <- st_read('./Input/contus_shp/') # for windows must be directory, might be something else for linux?
    #     contus <- st_read('./Input/contus_shp/usamapplot_best.shp')
    #     selpts <- land_grid_list$
    sel.files <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)
    selpts <- rbindlist(lapply(sel.files, function(tif){
        extnt <- extent(raster(tif))
        x <- mean(extnt[1], extnt[2])+50000
        y <- mean(extnt[3], extnt[4])+50000
        return(data.table(x, y))
    }))
    #     selpts <- ldsel.nnd[,.(x,y)]
    tile.one <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)[1]
    #     center.coords <- as.data.table(matrix(unlist(lapply(tile.list, get.centers)), ncol=2, byrow=TRUE))
    #     setnames(center.coords, c('x','y'))
    tile.potential <- list.files('./Landscape_Setup/Pipeline_SSF_Weekly/4_Output/indiv_plands', full.names=TRUE)
    pot.pts <- rbindlist(lapply(tile.potential, function(tif){
        extnt <- extent(raster(tif))
        x <- mean(extnt[1], extnt[2])+50000
        y <- mean(extnt[3], extnt[4])+50000
        return(data.table(tif, x, y))
    }))

    rast.crs <- crs(raster(tile.one))
    contus.transform <- st_transform(contus[1], rast.crs)
    #     bounds <- raster::extent(raster(tile.one))
    bounds <- st_bbox(contus.transform)
    xrng <- bounds[c(1,3)]
    yrng <- bounds[c(2,4)]

    map.scale = 100
    ysz <- round(log(yrng[2] - yrng[1]) * map.scale)
    xsz <- round(log(xrng[2] - xrng[1]) * map.scale)
    png('./Output/figures/tile_map.png', width=1600, height=1200)
    plot(pot.pts[,.(x,y)], col=NA, pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n')
    polys(contus.transform, col=0, lwd=2, alpha=1, fill=NA)
    points(selpts, pch=16, cex=2.3, col=2)
    points(pot.pts[,.(x,y)], col='gray', pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n')
    #     points(pot.pts[,.(x,y)], col=1, pch=1, cex=1.3)
    dev.off()
}
