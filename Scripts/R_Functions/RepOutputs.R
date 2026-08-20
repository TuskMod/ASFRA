## repackage outputs after simulation/burn-in run

rep_outputs <- function(out.list, v, l, r, parameters, out.opts, prevrep.in = as.list(rep(NA, 5))){

    list2env(parameters, .GlobalEnv)

    # if it is the burnin, the r=0 below will overwrite NA defaults in this irrelevant
    # if it is after the burn-in, these will have existing values that new values will be added to
    print("starting repoutputs")
    tm.mat <- prevrep.in[[1]]
    summ.vals <- prevrep.in[[2]]
    incidence <- prevrep.in[[3]]
    detections <- prevrep.in[[4]]
    allzone <- prevrep.in[[5]]
    ## test these outputs if out.opts doesn't include them
    ## id by variable combination, landscape, rep, and timestep for one row per timestep data
    end.tm <- out.list$endtime[[1]]
    id.r <- c(v,l,r)
    tm.mat.r <- cbind(matrix(id.r, nrow=end.tm, ncol=3, byrow=TRUE), seq(end.tm))
    colnames(tm.mat.r) <- c("var","land","rep","timestep")
    #Handle effective removal rate (timestep output)
    tm.mat.r <- cbind(tm.mat.r, out.list$Ct[seq(end.tm)])
    colnames(tm.mat.r)[ncol(tm.mat.r)] <- "Ct"
    # Births
    tm.mat.r <- cbind(tm.mat.r, out.list$BB[seq(end.tm)])
    colnames(tm.mat.r)[ncol(tm.mat.r)] <- 'BB'

    #Handle sounderlocs (optional...)
    ## seems like other things were supposed to happen in sounderlocsSummarize, if we want to use those this will have to change
    if ("sounderlocs" %in% out.opts){
        print("about to save sounder locs!")
        solocs.r <- sounderlocsSummarize(out.list$sounderlocs, r)[[1]]
        solocs.all <- out.list$sounderlocs
        print("saved sounder locs!")
#         solocs.all <- cbind(v, l, r, out.list$sounderlocs)
        print(solocs.r)
        print(tm.mat.r)
        tm.mat.r <- cbind(tm.mat.r, solocs.r[3:8])
        print("column bind!")
        setDT(solocs.all)
        print("finished sounder locs!")
        
    }
    # detections (optional...)
    ## has a row for each timestep AND detection type, with timestep, code (1=live,0=dead), number of individuals detected, and position
    ## still need to test sample = 1
    if ('alldetections' %in% out.opts){
        print("starting detection saving!")
        print(out.list$alldetections)
       # detections.r <- as.matrix(out.list$alldetections)
        detections.r <- out.list$alldetections
       
        detections.r <- detections.r[detections.r$time <= end.tm,]
        
        replication_times <- rep(3,dim(tm.mat.r)[[1]])
       
        dup_tm <- tm.mat.r[rep(row.names(tm.mat.r), times = replication_times), ]
        #print(dup_tm)
        n.det <- nrow(detections.r)
        #print(n.det)
        #print("done")
        detections.r <- suppressWarnings(cbind(matrix(id.r, ncol=3, nrow=n.det, byrow=TRUE), detections.r)) # gave a warning if there were no detections; very annoying
        detections.r <- cbind(dup_tm,detections.r)
        #colnames(detections.r) <- c('var', 'land', 'rep', 'timestep', 'code', 'detected', 'loc')
        #allzone.r <- out.list$allzonecells
        #colnames(allzone.r) <- c('var', 'land', 'rep', 'timestep', 'loc')
        print("about to save detections!")
        setDT(as.data.frame(detections.r))
        
    }
    # incidence -- more rows than timesteps, separate output (optional...)
    if ('incidence' %in% out.opts){ # don't want this if it's a burn-in output
        incidence.r <- out.list$incidence
        n.inc = nrow(incidence.r)
        incidence.r <- suppressWarnings(cbind(matrix(id.r, ncol=3, nrow=n.inc, byrow=TRUE), incidence.r))
        colnames(incidence.r) <- c('var', 'land', 'rep', 'timestep', 'state', 'loc')
        setDT(incidence.r)
        
    }

    # single value per vlr combination outputs
    # Tinc # sumTculled # Mspread # IConDD # ICatDD # TincToDD # TincFromDD # DET (total detections) #
    summ.vals.r <- matrix(c(id.r, unlist(out.list[c(1, 2, 4:9, length(out.list))])), nrow=1)

    colnames(summ.vals.r) <- c('var','land','rep',names(out.list[c(1, 2, 4:9, length(out.list))]))
    setDT(tm.mat.r)
     #incidence.r[,state := as.numeric(factor(state, levels=c('exposed','infected','carcass')))]
    fwrite(tm.mat.r[,.(BB,S,E,I,R,C,Z)], paste0('./Output/tm.mat/tm.mat_r',r,'_l',l,'_v',v,'.gz'), compressLevel=4L)
#     fwrite(summ.vals.r, paste0('./Output/summ.vals/summ.vals_r',r,'_l',l,'_v',v,'.csv'))
     #fwrite(incidence.r[,.(timestep, state, loc)], paste0('./Output/incidence/incidence_r',r,'_l',l,'_v',v,'.gz'))
     fwrite(detections.r, paste0('./Output/detections/detections_r',r,'_l',l,'_v',v,'.csv'))
#     fwrite(allzone, paste0('./Output/allzone/allzone_r',r,'_l',l,'_v',v,'.csv'))
    fwrite(solocs.all[,.(time, cell, pref, S, E, I, R, C, Z)], paste0('./Output/solocs.all/solocs.all_r',r,'_l',l,'_v',v,'.gz'), compressLevel=9L)

#     return(list(tm.mat, summ.vals, incidence, detections, allzone, solocs.all))
    return(summ.vals.r)
}
