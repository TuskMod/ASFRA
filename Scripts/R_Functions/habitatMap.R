
# compiles habitat maps for scaling preference values to determine initial sounder range
habitatmap <- function(tiles){
    tiles_merge <- terra::merge(terra::sprc(tiles))
    terra::writeRaster(tiles_merge, './Input/lwr48_habitat.tiff', filetype='GTiff')
    return('./Input/lwr48_habitat.tiff')
}
