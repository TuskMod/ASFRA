## est ridge with landscape-based cv
# make_rslt <- function(rslt, variables, mv.params){
make_rslt <- function(variables, mv.params){
#     library(targets)
#     library(glmnet)
#     library(data.table)

#     rslt <- as.data.table(t(fread('result.outputs.csv')))
    rslt <- fread('result.outputs.csv')

    setnames(rslt, c('v','l','r','est','edge.tm','max.dist','inf.area','inf.spd','sounder.weeks','prop.infd','max.inc','tm.esc'))
    rslt[est == 2, est := 1]
    setDT(variables)
    variables[,v := 1:.N]
    setkey(variables, v)
    setkey(mv.params, index)
    setkey(rslt, v, l, r)
    rslt <- rslt[variables, on=.(v)]
    rslt <- rslt[mv.params, on=.(l=index)]
    rslt[, names(.SD) := lapply(.SD, as.factor), .SDcols=c('est', 'contact', 'variant', 'density')]

    keep.names <- c("v", "l", "r", "est", "inf.area", "inf.spd", "sounder.weeks", "prop.infd", "max.inc", "tm.esc",
                    "contact", "variant", "density", "nnd_med", "nnd_range", "nnd_mean", "nnd_sd",
                    "moranI", "gearyC", "tc", "mast", "rgd", "rd1","rd2","rd3", "prcp", "tmin", "tmax",
                    "drt", "contag", "aggindex", "entropy", "simpindx", "nnd_cv", "mvmt.mean", "mvmt.cv")
    rslt <- rslt[,.SD, .SDcols=keep.names]

    # so that it can group by land ID for glmm
    rslt[,l := as.factor(l)]

    return(rslt)
}





glmm.lasso <- function(lasslist, dat, land.kf,
                          fix = 'as.factor(contact) + as.factor(variant) + as.factor(density) + nnd_med + nnd_range + nnd_mean + nnd_sd + moranI + gearyC + tc + mast + rgd + rd1 + rd2 + rd3 + prcp + tmin + tmax + drt + contag + aggindex + entropy + simpindx + nnd_cv + mvmt.mean + mvmt.cv',
                          rnd = 'l',
                          cv=TRUE
                          ){
    lambda <- lasslist[,lambda]
    kf <- lasslist[,kf.group]
    rsp.var <- lasslist[,rsp.var]

    overdisp <- FALSE
    if (rsp.var != 'est'){
        dat <- dat[est == 1,]
    }
    if (rsp.var %in% c('sounder.weeks','max.inc','tm.esc')){
        fam <- poisson() # apparently quasi-x families don't exist in glmmLasso
        overdisp <- TRUE
    } else if (rsp.var == 'prop.infd') {
        fam <- binomial()
        overdisp <- TRUE
    } else if (rsp.var == 'est') {
        fam <- binomial()
        dat[,est := as.numeric(as.character(est))]
    } else if (rsp.var == 'inf.spd') {fam <- gaussian()}
    if (rsp.var == 'tm.esc'){
        stop('this ain\'t right (fix this later)')
    }
    library(glmmLasso)
    in.lands <- land.kf[kf.group == kf,land]
    trn.dat <- dat[l %in% in.lands == FALSE, .SD, .SDcols=c(rsp.var, 'l', 'contact', 'variant', 'density', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd', 'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rd1','rd2','rd3', 'prcp', 'tmin', 'tmax', 'drt', 'contag', 'aggindex', 'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean','mvmt.cv')]
    tst.dat <- dat[l %in% in.lands, .SD, .SDcols=c(rsp.var, 'l', 'contact', 'variant', 'density', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd', 'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rd1','rd2','rd3', 'prcp', 'tmin', 'tmax', 'drt', 'contag', 'aggindex', 'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean','mvmt.cv')]
    fix2 <- formula(paste(rsp.var, '~', fix))
    rnd2 <- list(~1)
    names(rnd2) <- rnd
    lassoi <- glmmLasso(fix = fix2, rnd = rnd2, data=data.frame(trn.dat), lambda=lambda, family=fam, control=glmmLassoControl(overdispersion=overdisp), final.re=TRUE)
    if (cv == TRUE){
        y.pred <- predict(lassoi, newdata=tst.dat)
        rmse <- sqrt(sum((tst.dat[,..rsp.var] - y.pred)^2))
        outlist <- list(c(rsp.var, kf, lambda, rmse, lassoi$phi), lassoi)
    } else {
        outlist <- list(lasslist, lassoi)
    }
    return(outlist)
}

glmm.vanilla <- function(rsp.var, dat, #land.kf,
                          fix = '(contact : variant : density) * (nnd_med + nnd_range + nnd_mean + nnd_sd + moranI + gearyC + tc + mast + rgd + rd1 + rd2 + rd3 + prcp + tmin + tmax + drt + contag + aggindex + entropy + simpindx + nnd_cv + mvmt.mean + mvmt.cv)',
                          rnd = 'l'
                          ){
#     fix = 'nnd_med + nnd_range + nnd_mean + nnd_sd + moranI + gearyC + tc + mast + rgd + rd1 + rd2 + rd3 + prcp + tmin + tmax + drt + contag + aggindex + entropy + simpindx + nnd_cv + mvmt.mean + mvmt.cv + (1|contact) + (1|variant) + (1|density)'
#     rnd = 'l/(contact + variant + density)'
#     lambda <- lasslist[,lambda]
#     kf <- lasslist[,kf.group]
#     rsp.var <- lasslist[,rsp.var]
#     dat[,density := as.numeric(as.character(density))]

    zi_form <- ~0
#     zi_form <- ~.
    overdisp <- FALSE
    if (rsp.var != 'est'){
        dat <- dat[est == 1,]
#         fam <- gaussian()
    }
#     dat[,density := as.numeric(as.character(density))]
    print(rsp.var)
    if (rsp.var %in% c('sounder.weeks','max.inc','tm.esc')){
        fam <- nbinom2()
        overdisp <- TRUE
    } else if (rsp.var == 'prop.infd') {
        fam <- beta_family()
        overdisp <- TRUE
        dat[prop.infd == 1, prop.infd := 0.999999]
    } else if (rsp.var == 'est') {
        fam <- binomial()
        dat[,est := as.numeric(as.character(est))]
    } else if (rsp.var == 'inf.spd') {
        fam <- lognormal()
#         browser()
#         zi_form <- ~.
        dat[inf.spd <0.0001, inf.spd := 0.0001]
    }
    if (rsp.var == 'max.inc'){
        dat[is.infinite(max.inc), max.inc := 0]
#         browser()
    }
#     in.lands <- land.kf[kf.group == kf,land]
#     trn.dat <- dat[l %in% in.lands == FALSE, .SD, .SDcols=c(rsp.var, 'l', 'contact', 'variant', 'density', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd', 'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rd1','rd2','rd3', 'prcp', 'tmin', 'tmax', 'drt', 'contag', 'aggindex', 'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean','mvmt.cv')]
#     tst.dat <- dat[l %in% in.lands, .SD, .SDcols=c(rsp.var, 'l', 'contact', 'variant', 'density', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd', 'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rd1','rd2','rd3', 'prcp', 'tmin', 'tmax', 'drt', 'contag', 'aggindex', 'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean','mvmt.cv')]
    fix2 <- formula(paste(rsp.var, '~', fix, '+ (1|', rnd, ')'))
#     rnd2 <- list(~1)
#     names(rnd2) <- rnd
    lassoi <- glmmTMB(formula=fix2, data=data.frame(dat), family=fam, ziformula=zi_form)#, control=glmmLassoControl(overdispersion=overdisp), final.re=TRUE)
    browser()
#     lassoi <- glmmLasso(fix = fix2, rnd = rnd2, data=data.frame(trn.dat), lambda=lambda, family=fam, control=glmmLassoControl(overdispersion=overdisp), final.re=TRUE)
#     if (cv == TRUE){
#         y.pred <- predict(lassoi, newdata=tst.dat)
#         rmse <- sqrt(sum((tst.dat[,..rsp.var] - y.pred)^2))
#         outlist <- list(c(rsp.var, kf, lambda, rmse, lassoi$phi), lassoi)
#     } else {
#     outlist <- list(lasslist, lassoi)
#     }
    return(list(lassoi))
}


make_lkf <- function(rslt, kf.sz = 10){
    if (kf.sz == 1){
        print('leave one out CV')
        nlands <- length(unique(rslt[,l]))
        land.kf.group <- data.table(land = unique(rslt[,l]), kf.group=unique(rslt[,unclass(l)]))
    } else if (kf.sz == 0){
        print('no CV')
        land.kf.group <- data.table(land = 1, kf.group=kf.sz)
    } else {
        print(paste('CV of', kf.sz, 'folds'))
        land.kf.group <- data.table(land = unique(rslt[,l]), kf.group=sample(rep(seq(kf.sz), length.out=length(unique(rslt[,l]))), replace=FALSE))
    }
    return(land.kf.group)
}

make.lasslist <- function(rslt, rsp.var, lkf.table, lambdas=10, oom.lambda.min=-4, oom.lambda.max=5, cvlbd=NULL){

    if (lambdas != 1){
        lasslist <- CJ(lambda = 10^seq(oom.lambda.max, oom.lambda.min, length.out=lambdas), land=unique(rslt[,l]), rsp.var=rsp.var)
        lasslist <- lkf.table[lasslist, on='land']
        lasslist <- unique(lasslist[,-'land'])
    } else {
        cvlbd <- rbindlist(cvlbd)
        avg.rmse <- cvlbd[, mean(rmse), by=.(resp.var, lambda)]
        avg.rmse[,min.rmse := min(V1), by=.(resp.var)]
        min.rmse <- avg.rmse[min.rmse == V1,]
        min.rmse[,kf.group := 1]
        lasslist <- min.rmse[,.(kf.group, lambda, rsp.var)]
    }
    return(lasslist)
}



organize.lasso <- function(cv.out){
    tab.index <- cv.out[[1]]
    cv.table <- data.table(t(tab.index))
    setnames(cv.table, c('resp.var','kf','lambda','rmse','overdisp'))
    cv.table[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('lambda', 'rmse','overdisp')]
    return(cv.table)
}

cv.curves <- function(cv.table){

    cv.table <- rbindlist(cv.table)
    lambda_means <- cv.table[, mean(rmse), by=.(lambda, resp.var)]
    rv.list <- unique(cv.table[,resp.var])
    png('./Output/figures/cv_curves.png', width=1200, height=1100)
    par(mfrow=c(2, ceiling(length(rv.list)/2)))
    lapply(rv.list, function(x){
        plot(rmse ~ log(lambda), data=cv.table[resp.var==x,], col=as.factor(kf), main=x, cex=1.6)
        lines(V1 ~ log(lambda), data=lambda_means[resp.var==x,], lwd=2)
        abline(h=min(lambda_means[resp.var==x,V1]), col=2, lty=3)
        abline(v=log(lambda_means[V1==min(lambda_means[resp.var==x,V1]) & resp.var==x,lambda]), col=2, lty=3)
    })
    dev.off()
}



make_vector_lyr <- function(var.set, preds, grid_path){
    # really just need to tie the values to a vector object
    grid <- vect(grid_path, crs='EPSG:5070')[,c('ID')]
    # figure out which columns have response variables in them
#     resp.cols <- names(preds)[grep('*_pred$', names(preds))]
    # cast the table to wide
#     preds_cast <- dcast(preds, l +inside + edge~ var, value.var = resp.cols)
    # get rid of annoying naming leftovers
#     setnames(preds_cast, gsub('pred_','',names(preds_cast)))
    # not necessary, but makes the output cleaner
#     setnames(preds_cast, 'l', 'ID')
    # connect predictions to spatvector
    pred_geo <- merge(grid, preds[var == var.set,], by.x='ID', by.y='l')
    return(pred_geo)
}


oos.preds.glmm <- function(predmods, resp.var, tiledat){
#     resp.var <- as.character(predmods$fix[2])
# #     lambda.mins <- lasslist2[rsp.var == resp.var,lambda]
#     #     if(resp.var != lasslist2[,rsp.var]) stop('response variables are not matching with lambdas')
#     #     if (resp.var == 'exc') {
#     #         tiledat <- tiledat[CJ(l=unique(tiledat[,l]), tm=1:78), on=.(l), allow.cartesian=TRUE]
#     #     }
#     tiledat.trim <- tiledat[,.(l, contact, variant, density, nnd_med, nnd_range, nnd_mean, nnd_sd, moranI, gearyC, tc, mast, rgd, rds, prcp, tmin, tmax, drt, contag, aggindex, entropy, simpindx, nnd_cv, mvmt.mean, mvmt.cv)]
#     tiledat.trim[,l := as.factor(l)]
#     tiledat.trim[,as.character(resp.var) := 0]
#     pred <- predict(predmods, newdata=tiledat.trim)
#     coefs.vals <- coef(predmods)
#     coefs.used <- names(coefs.vals)[which(coefs.vals != 0)]
#     pred.table <- cbind(tiledat, resp.variable=resp.var, coefs.used=paste(coefs.used, collapse=' '), lambda.nom=lambda.mins, lambda.mod=predmods$lambda.max, pred)
# if (resp.var == 'inf.spd') browser()
    tiledat[,as.character(resp.var) := predict(predmods[[1]], newdata=tiledat, type='response', re.form=NA)]
    return(tiledat)
}


edge.correction <- function(td){
    tdi <- td[inside == 1,]
    nnd.names <- names(td)[grep('^nnd', names(td))]
    nnd.fncns <- tstrsplit(nnd.names, '_', keep=2)
    nnd.2si <- lapply(nnd.names, function(nms){
        mn <- mean(unlist(tdi[,..nms]))
        sd <- sd(unlist(tdi[,..nms]))
        up <- mn + 2*sd
        dn <- mn - 2*sd
        return(c(nms, mn, sd, up, dn))
    })
    nnd.2si <- data.table(matrix(unlist(nnd.2si), ncol=5, byrow=TRUE))
    nnd.2si[,names(.SD) := lapply(.SD, as.numeric), .SDcols=2:5]
    setnames(nnd.2si, c('metric','mean', 'sd','upper','lower'))
    td[inside==0,nnd_med := pmin(nnd_med, nnd.2si[metric=='nnd_med', upper])]
    td[inside==0,nnd_med := pmax(nnd_med, nnd.2si[metric=='nnd_med', lower])]
    td[inside==0,nnd_range := pmin(nnd_range, nnd.2si[metric=='nnd_range', upper])]
    td[inside==0,nnd_range := pmax(nnd_range, nnd.2si[metric=='nnd_range', lower])]
    td[inside==0,nnd_mean := pmin(nnd_mean, nnd.2si[metric=='nnd_mean', upper])]
    td[inside==0,nnd_mean := pmax(nnd_mean, nnd.2si[metric=='nnd_mean', lower])]
    td[inside==0,nnd_sd := pmin(nnd_sd, nnd.2si[metric=='nnd_sd', upper])]
    td[inside==0,nnd_sd := pmax(nnd_sd, nnd.2si[metric=='nnd_sd', lower])]
    td[inside==0,nnd_cv := pmin(nnd_cv, nnd.2si[metric=='nnd_cv', upper])]
    td[inside==0,nnd_cv := pmax(nnd_cv, nnd.2si[metric=='nnd_cv', lower])]
    return(td)
}
