#Updates needed here:
  #way inputs are put into Make Grid is clunky... want similar format of object, grid.opt
  #add sprc and ras input options
  #change checks in MakeGrid to look for SpatRaster format raster than raster
    #for now just going to convert types back/forth. shouldnt be too bad bc just run need to run this once
  #add way to trigger multiple homogenous or heterogeneous grids
  #pass len and inc to grid.list


#object=c(parameters$len,parameters$inc)
#grid.opt=parameters$grid.opts
InitializeGrids <- function(path, parameters0){
    pop_init_grid_opts <- parameters0$pop_init_grid_opts
    grid.opts <- parameters0$grid.opts
    len <- parameters0$len
    inc <- parameters0$inc
    km_len <- len*inc

    # previously in ReadLands.R
    # finds land tiles from path object and collects into sprc object
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
#     return(plands_sprc)

    # previously in _targets.R
    if (pop_init_grid_opts == 'homogeneous'){
        if(grid.opts != 'ras'){ # if grid.opts is homogeneous or heterogeneous
            # make a grid either uniform or random with even initial pig locations
            land_grid_list <- InitializeGrids_sub(c(len, inc), grid.opts)
        } else if (grid.opts == 'ras'){ # if there is an input raster
            inc <- terra::res(plands_sprc[1])[1]/1000
            km_len <- dim(plands_sprc[1])[1]*inc
            land_grid_list <- InitializeGrids_sub(plands_sprc, grid.opts)
        }
    } else if (pop_init_grid_opts == 'heterogeneous'){
        # make a grid with uneven pig initial locations...
        if (grid.opts == 'homogeneous') {
            # can't do neutral plane with random pig distribution
            stop('Cannot run homogeneous grid.opts with heterogeneous pop_init_grid_opts')
        } else if (grid.opts == 'heterogeneous'){
            # random pig distribution with random landscape
            land_grid_list <- InitializeGrids_sub(c(len, inc), grid.opts)
        } else if (grid.opts == 'ras'){
            # random pig distribution with raster landscape
            inc <- terra::res(plands_sprc[1])[1]/1000
            km_len <- dim(plands_sprc[1])[1]*inc
            land_grid_list <- InitializeGrids_sub(plands_sprc, grid.opts)
        }
    }

    return(list(land_grid_list, inc, km_len))
}

    # modified to include ReadLands.R in this targets step
    # not sure if this works for homogenous, random, or individual land tiles anymore, but skips the need for geotargets so it can operate in a conda environment (hopefully)
InitializeGrids_sub <- function(object, grid.opts){
    if(class(object)=="SpatRaster"){
        tar.grid.list <- vector(mode="list", length=1)
        grid.list <- Make_Grid(object, grid.opts)
        tar.grid.list[[1]] <- grid.list
    }
    if(class(object)=="SpatRasterCollection"){
        tar.grid.list <- vector(mode="list", length=length(object))
        for(g in 1:length(object)){
            ras <- object[g]
            grid.list <- Make_Grid(ras, grid.opts)
            tar.grid.list[[g]] <- grid.list
        }
    }
    #vector of len, inc
    if(class(object)=="numeric"){
        if(length(object)==3){
            tar.grid.list <- vector(mode="list", length=object[3])
            for(g in 1:object[3]){
            grid.list <- Make_Grid(c(object[1], object[2]), grid.opts)
            tar.grid.list[[g]] <- grid.list
            }
        } else{
            tar.grid.list <- vector(mode="list", length=1)
            grid.list <- Make_Grid(c(object[1], object[2]), grid.opts)
            tar.grid.list[[1]] <- grid.list
        }
    }
    return(tar.grid.list)
}


