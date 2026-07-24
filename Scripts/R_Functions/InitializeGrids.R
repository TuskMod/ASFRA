InitializeGrids <- function(path, indv_ras_data,parameters0){
    # pulls user-supplied landscape attributes (if no raster present, these are used)
    pop_init_grid_opts <- parameters0$pop_init_grid_opts
    grid.opts <- parameters0$grid.opts
    len <- parameters0$len
    inc <- parameters0$inc
    km_len <- len*inc
    

    # previously in ReadLands.R
    # finds land tiles from path object and collects into sprc object
    # organizes lands by name (i.e. id number of the original tiles)
    
    if (grid.opts == 'ras'){
        fs <- list.files(path, full.names=TRUE)
        nm <- unlist(tstrsplit(fs, '/', keep=5))
        nm <- unlist(tstrsplit(nm, '[_.]', keep=2))
        plands_list <- vector(mode="list", length=length(fs))
        for(fi in 1:length(fs)){
            plands_list[[fi]] <- terra::rast(fs[fi])
            names(plands_list[[fi]]) <- nm[fi]
        }
        # stick lands into a sprc object
        plands_sprc <- terra::sprc(plands_list)
        names(plands_sprc) <- lapply(plands_list, names)
    }

    # previously in _targets.R
    # handle the variants -- homogeneous, random, or raster-based landscape
    if (pop_init_grid_opts == 'homogeneous'){ # i.e. initial sounders are distributed randomly
        if(grid.opts != 'ras'){ # if grid.opts is homogeneous or heterogeneous
            # make a grid either uniform or random with even initial pig locations
            land_grid_list <- InitializeGrids_sub(c(len, inc), grid.opts)
        } else if (grid.opts == 'indvras'){ # if there is an input raster
            inc <- terra::res(indv_ras_data)[1]/1000
            km_len <- dim(indv_ras_data)[1]*inc
            land_grid_list <- InitializeGrids_sub(indv_ras_data, grid.opts)
        }
    } else if (pop_init_grid_opts == 'heterogeneous'){ # i.e. initial sounders are distributed according to landscape preferences
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
        } else if (grid.opts == "indvras"){
          # random pig distribution with raster landscape
          inc <- terra::res(indv_ras_data[1])[1]/1000
          km_len <- dim(indv_ras_data[1])[1]*inc
          land_grid_list <- InitializeGrids_sub(indv_ras_data[1], grid.opts)
        }
    }
    return(list(land_grid_list, inc, km_len))
}

    # modified to include ReadLands.R in this targets step
    # not sure if this works for homogeneous, random, or individual land tiles anymore, but skips the need for geotargets so it can operate in a conda environment (hopefully)
InitializeGrids_sub <- function(object, grid.opts){
    if(class(object)=="SpatRaster"){ # using one raster
        tar.grid.list <- vector(mode="list", length=1)
        grid.list <- Make_Grid(object, grid.opts)
        tar.grid.list[[1]] <- grid.list
    }
    if(class(object)=="SpatRasterCollection"){ # using 2+ rasters
        tar.grid.list <- vector(mode="list", length=length(object))
        for(g in 1:length(object)){
            class(g)
            ras <- object[g]
            grid.list <- Make_Grid(ras, grid.opts)
            tar.grid.list[[g]] <- grid.list
        }
    }
    #vector of len, inc
    if(class(object)=="numeric"){ # using a generated landscape (homogeneous or random)
        if(length(object)==3){
            tar.grid.list <- vector(mode="list", length=object[3])
            for(g in 1:object[3]){ # random landscape
                grid.list <- Make_Grid(c(object[1], object[2]), grid.opts)
                tar.grid.list[[g]] <- grid.list
            }
        } else{ # just homogeneous
            tar.grid.list <- vector(mode="list", length=1)
            grid.list <- Make_Grid(c(object[1], object[2]), grid.opts)
            tar.grid.list[[1]] <- grid.list
        }
    }
    return(tar.grid.list)
}


