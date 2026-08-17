

combo.plans <- function(parameters, variables, reps, mv.parms){
    # defines all possible combinations of replicate, variables, and landscapes

    library(data.table)
    # if any output folders are missing, create them
    dirlist <- list.dirs()
    if ('./Output/tm.mat' %in% dirlist == FALSE){dir.create('./Output/tm.mat')}
    if ('./Output/summ.vals' %in% dirlist == FALSE){dir.create('./Output/summ.vals')}
    if ('./Output/incidence' %in% dirlist == FALSE){dir.create('./Output/incidence')}
        if ('./Output/detections' %in% dirlist == FALSE){dir.create('./Output/detections')}
#         if ('./Output/allzone' %in% dirlist == FALSE){dir.create('./Output/allzone')}
    if ('./Output/solocs.all' %in% dirlist == FALSE){dir.create('./Output/solocs.all')}

    setDT(mv.parms)
    # table of all desired simulation runs
    lvtable <- CJ(vars = seq(nrow(variables)), land = unique(mv.parms[,index]), rep = seq(reps))
    return(lvtable)
}

check.existing <- function(lvtable, out.repl){
    # checks if output files exist for existing simulation runs
    tm.mat.in <- list.files('./Output/tm.mat')
    detections.in <- list.files('./Output/detections')
    
#     summ.vals.in <- list.files('./Output/summ.vals')
#     incidence.in <- list.files('./Output/incidence')
    solocs.all.in <- list.files('./Output/solocs.all')
  #  detections <- list.files('./Output/detections')
    # if out.repl is not 1/TRUE, and there are things already in the output directories,
    # sort existing things into a table and find what is missing relative to lvtable
#     if (out.repl != TRUE & length(c(tm.mat.in, summ.vals.in, incidence.in, solocs.all.in)) != 0){
    if (out.repl == 0 & length(c(tm.mat.in, solocs.all.in,detections.in)) != 0){
        print("checking")
        splt.check <- function(nlst, nm){
            # slice and dice the names of existing files to get rep, var, and land id's
            instr <- as.data.table(tstrsplit(nlst, '_', keep=2:4))
            instr <- instr[,lapply(.SD, function(x) unlist(regmatches(x, gregexpr('(\\d+)', x))))]
            setnames(instr, c('r','l','v'))
            # indicate thet the combination has existing output files
            instr[,tmp := 1]
            setnames(instr, 'tmp', nm)
            return(instr)
        }
        # check for output files separately for each type of output file
        tm.tab <- unique(as.data.table(rbindlist(lapply(tm.mat.in, splt.check, nm = 'tm.mat'))))
#         summ.tab <- unique(as.data.table(rbindlist(lapply(summ.vals.in, splt.check, nm = 'summ.vals'))))
#         incid.tab <- unique(as.data.table(rbindlist(lapply(incidence.in, splt.check, nm = 'incidence'))))
        solocs.tab <- unique(as.data.table(rbindlist(lapply(solocs.all.in, splt.check, nm = 'solocs.all'))))
        detections.tab <- unique(as.data.table(rbindlist(lapply(detections.in, splt.check, nm = 'detection'))))
        print(solocs.tab)
        # connect all filetype lists together
#         bndtab <- merge(tm.tab, summ.tab, by=c('v','l','r'), all=TRUE)
        bndtab <- merge(tm.tab, solocs.tab, by=c('v','l','r'), all=TRUE)
#         bndtab <- merge(bndtab, incid.tab, by=c('v','l','r'), all=TRUE)
#         bndtab <- merge(bndtab, solocs.tab, by=c('v','l','r'), all=TRUE)
        # anything that has a missing file goes to 0
        bndtab[is.na(bndtab)] <- 0
        bndtab[,names(.SD) := lapply(.SD, as.numeric)]
        # connect existing file list to goal file list, and filter out lines that have all file types accounted for
        lvtable2 <- merge(lvtable, bndtab, by.x=c('vars', 'land', 'rep'), by.y=c('v', 'l', 'r'), all=TRUE)
        lvtable2 <- lvtable2[is.na(tm.mat) | is.na(solocs.all),]
#         lvtable2 <- lvtable2[is.na(tm.mat) | is.na(summ.vals) | is.na(incidence) | is.na(solocs.all),]
        lvtable2 <- lvtable2[,.(vars, land, rep)]

    } else {
        # just use the goal file as the list of simulations to run
        lvtable2 <- lvtable
    }

    print(lvtable2)
    return(lvtable2)
}
