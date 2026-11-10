## prepare density input files for integration into the simulation workflow

dens.prep <- function(infile, target_raster){
    orig <- terra::rast(infile)
    proj <- terra::project(orig, terra::rast(target_raster))
    outname <- gsub('.tif','_5070.tif', infile)
    terra::writeRaster(proj, outname, filetype='GTiff')
    return(outname)
}


habitatmap <- function(tiles){
    tiles_merge <- terra::merge(terra::sprc(tiles))
    terra::writeRaster(tiles_merge, './Input/lwr48_habitat.tiff', filetype='GTiff')
    return('./Input/lwr48_habitat.tiff')
}
