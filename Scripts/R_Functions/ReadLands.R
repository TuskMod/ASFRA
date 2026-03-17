ReadLands <- function(path){
    fs <- list.files(path,full.names=TRUE)
    plands_list <- vector(mode="list",length=length(fs))
    nm <- unlist(tstrsplit(fs, '/', keep=5))
    nm <- unlist(tstrsplit(nm, '[_.]', keep=2))
    for(f in 1:length(fs)){
        plands_list[[f]] <- terra::rast(fs[f])
        names(plands_list[[f]]) <- nm[f]
    }
    plands_sprc <- terra::sprc(plands_list)
    names(plands_sprc) <- lapply(plands_list, names)
    return(plands_sprc)
}
