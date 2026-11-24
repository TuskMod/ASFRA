GetSurfaceParms<-function(parameters, inc, km_len, habminmax){
    if (parameters$grid.opts == 'ras'){
        # if it's a raster, surface parameters are driven by raster attributes from a previous step
        area <- km_len^2
        parameters$inc <- inc
        parameters <- append(parameters,km_len)
        names(parameters)[length(parameters)] <- "km_len"

        parameters <- append(parameters,area)
        names(parameters)[length(parameters)] <- "area"

        parameters <- append(parameters, habminmax[1:2])
        names(parameters)[length(parameters)-1] <- "habmin"
        names(parameters)[length(parameters)] <- "habmax"
    } else {
        km_len <- parameters$inc * parameters$len
        area <- km_len^2
        parameters <- append(parameters, km_len)
        names(parameters)[length(parameters)] <- "km_len"
        parameters <- append(parameters, area)
        names(parameters)[length(parameters)] <- "area"
        parameters$inc <- inc
    }

    return(parameters)

}
