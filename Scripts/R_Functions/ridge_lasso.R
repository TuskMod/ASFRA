## est ridge with landscape-based cv
make_rslt <- function(variables, mv.params){
    library(targets)
    library(glmnet)
    library(data.table)

    rslt <- fread('Output/result.outputs.csv')

    setnames(rslt, c('v','l','r','est','edge.tm','max.dist','inf.area','inf.spd','sounder.weeks','prop.infd','max.inc','tm.esc'))
    rslt[est == 2, est := 1]
    #variables <- tar_read(variables)
    #mv.params <- tar_read(mv.params)
    setDT(variables)
    setDT(mv.params)
    variables[,v := 1:.N]
    setkey(variables, v)
    setkey(mv.params, index)
    setkey(rslt, v, l, r)
    rslt <- rslt[variables, on=.(v)]
    rslt <- rslt[mv.params, on=.(l=index)]
    #     rslt[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('gamma.shape', 'gamma.scale', 'moranI', 'tc', 'mast', 'rgd', 'rds','dayl','prcp','drt','contag','aggindex','entropy','simpindx')]
    rslt[, names(.SD) := lapply(.SD, as.factor), .SDcols=c('est', 'contact', 'variant', 'density')]

    keep.names <- c("v", "l", "r", "est", "inf.area", "inf.spd", "sounder.weeks", "prop.infd", "max.inc", "tm.esc",
                    "contact", "variant", "density", "nnd_med", "nnd_range", "nnd_mean", "nnd_sd",
                    "moranI", "gearyC", "tc", "mast", "rgd", "rds", "dayl", "prcp", "tmin", "tmax",
                    "drt", "contag", "aggindex", "entropy", "simpindx", "nnd_cv", "mvmt.mean", "mvmt.cv")
    #                 "edge.tm", "max.dist",
    #                 "ss", "B1", "B2", "F1", "F2_int", "F2_B", "F2i_int", "F2i_B", "infpd", "incub",
    #                 "x", "y", "state", "gamma.shape", "gamma.scale",
    #                 "st.wt")
    rslt <- rslt[,.SD, .SDcols=keep.names]

    return(rslt)
}


ridge.lasso <- function(rslt, rslt1, tmtab3, variables, mv.params){

    # rslt <- tar_read(rslt)
    # rslt1 <- tar_read(rslt1)
    rslt1 <- rslt[est == 1,]
    rslt1[prop.infd == 1, prop.infd := 0.999]
    # tmtab3 <- tar_read(tmtab3)
    tmtab3 <- tm.tab.3(rslt1, variables, mv.params)
    keep.names <- c("v", "l", "r", "exc", "tm.esc","contact", "variant", "density", "nnd_med", "nnd_range",
                    "nnd_mean", "nnd_sd","moranI", "gearyC", "tc", "mast", "rgd", "rds", "dayl", "prcp",
                    "tmin", "tmax", "drt", "contag", "aggindex", "entropy", "simpindx", "nnd_cv", "mvmt.mean", "mvmt.cv")
    tmtab3 <- tmtab3[,.SD, .SDcols=keep.names]

    rslt  [,lno := unclass(as.factor(l))]
    rslt1 [,lno := unclass(as.factor(l))]
    tmtab3[,lno := unclass(as.factor(l))]
    tmtab3 <- tmtab3[exc==1,]

    glmdat.est           <- model.matrix(~.,rslt  [,.SD,      .SDcols=!c('v','l','r','est','sounder.weeks','prop.infd','max.inc','tm.esc','inf.spd','lno','inf.area')])[,-1]
    glmdat.sounder.weeks <- model.matrix(~.,rslt1 [,.SD,      .SDcols=!c('v','l','r','sounder.weeks','est','prop.infd','max.inc','tm.esc','inf.spd','lno','inf.area')])[,-1]
    glmdat.prop.infd     <- model.matrix(~.,rslt1 [,.SD,      .SDcols=!c('v','l','r','prop.infd','sounder.weeks','est','max.inc','tm.esc','inf.spd','lno','inf.area')])[,-1]
    glmdat.max.inc       <- model.matrix(~.,rslt1 [,.SD,      .SDcols=!c('v','l','r','max.inc','sounder.weeks','prop.infd','est','tm.esc','inf.spd','lno','inf.area')])[,-1]
    glmdat.tm.esc        <- model.matrix(~.,tmtab3[exc==1,.SD,.SDcols=!c('v','l','r','tm.esc','lno','exc')])[,-1]
    glmdat.inf.spd       <- model.matrix(~.,rslt1 [,.SD,      .SDcols=!c('v','l','r','inf.spd','sounder.weeks','prop.infd','max.inc','tm.esc','est','lno','inf.area')])[,-1]

    # lambda grid values
    grid <- 10^seq(10, -2, length=100)

    # overall ridge model
#     ridge.mod.est           <-    glmnet(glmdat.est,           as.numeric(rslt[,est]),            alpha=0, lambda=grid, family='binomial')
    lasso.mod.est           <-    glmnet(glmdat.est,           as.numeric(rslt[,est]),            alpha=1, lambda=grid, family='binomial', intercept=FALSE)
#     ridge.mod.sounder.weeks <-    glmnet(glmdat.sounder.weeks, as.numeric(rslt1[,sounder.weeks]), alpha=0, lambda=grid, family=quasipoisson(), control=list(maxit=999999))
    lasso.mod.sounder.weeks <-    glmnet(glmdat.sounder.weeks, as.numeric(rslt1[,sounder.weeks]), alpha=1, lambda=grid, family=quasipoisson(), control=list(maxit=999999), intercept=FALSE)
#     ridge.mod.prop.infd     <-    glmnet(glmdat.prop.infd,     as.numeric(rslt1[,prop.infd]),     alpha=0, lambda=grid, family=quasibinomial())
    lasso.mod.prop.infd     <-    glmnet(glmdat.prop.infd,     as.numeric(rslt1[,prop.infd]),     alpha=1, lambda=grid, family=quasibinomial(), intercept=FALSE)
#     ridge.mod.max.inc       <-    glmnet(glmdat.max.inc,       as.numeric(rslt1[,max.inc]),       alpha=0, lambda=grid, family=quasipoisson())
    lasso.mod.max.inc       <-    glmnet(glmdat.max.inc,       as.numeric(rslt1[,max.inc]),       alpha=1, lambda=grid, family=quasipoisson(), intercept=FALSE)
#     ridge.mod.tm.esc        <-    glmnet(glmdat.tm.esc,        as.numeric(tmtab3[,tm.esc]), alpha=0, lambda=grid, family=quasipoisson())
    lasso.mod.tm.esc        <-    glmnet(glmdat.tm.esc,        as.numeric(tmtab3[,tm.esc]), alpha=1, lambda=grid, family=quasipoisson(), intercept=FALSE)
#     ridge.mod.inf.spd       <-    glmnet(glmdat.inf.spd,       as.numeric(rslt1[,inf.spd]),       alpha=0, lambda=grid)
    lasso.mod.inf.spd       <-    glmnet(glmdat.inf.spd,       as.numeric(rslt1[,inf.spd]),       alpha=1, lambda=grid, intercept=FALSE)

#     cv.est.ridge            <- cv.glmnet(glmdat.est,           as.numeric(rslt[,est]),            alpha=0, lambda=grid, foldid=rslt[,lno],   family='binomial')
    cv.est.lasso            <- cv.glmnet(glmdat.est,           as.numeric(rslt[,est]),            alpha=1, lambda=grid, foldid=rslt[,lno],   family='binomial', intercept=FALSE)
#     cv.sounder.weeks.ridge  <- cv.glmnet(glmdat.sounder.weeks, as.numeric(rslt1[,sounder.weeks]), alpha=0, lambda=grid, foldid=rslt1[,lno],  family=quasipoisson(), control=list(maxit=999999))
    cv.sounder.weeks.lasso  <- cv.glmnet(glmdat.sounder.weeks, as.numeric(rslt1[,sounder.weeks]), alpha=1, lambda=grid, foldid=rslt1[,lno],  family=quasipoisson(), control=list(maxit=999999), intercept=FALSE)
#     cv.prop.infd.ridge      <- cv.glmnet(glmdat.prop.infd,     as.numeric(rslt1[,prop.infd]),     alpha=0, lambda=grid, foldid=rslt1[,lno],  family=quasibinomial())
    cv.prop.infd.lasso      <- cv.glmnet(glmdat.prop.infd,     as.numeric(rslt1[,prop.infd]),     alpha=1, lambda=grid, foldid=rslt1[,lno],  family=quasibinomial(), intercept=FALSE)
#     cv.max.inc.ridge        <- cv.glmnet(glmdat.max.inc,       as.numeric(rslt1[,max.inc]),       alpha=0, lambda=grid, foldid=rslt1[,lno],  family=quasipoisson())
    cv.max.inc.lasso        <- cv.glmnet(glmdat.max.inc,       as.numeric(rslt1[,max.inc]),       alpha=1, lambda=grid, foldid=rslt1[,lno],  family=quasipoisson(), intercept=FALSE)
#     cv.tm.esc.ridge         <- cv.glmnet(glmdat.tm.esc,        as.numeric(tmtab3[,tm.esc]),       alpha=0, lambda=grid, foldid=tmtab3[,lno], family=quasipoisson())
    cv.tm.esc.lasso         <- cv.glmnet(glmdat.tm.esc,        as.numeric(tmtab3[,tm.esc]),       alpha=1, lambda=grid, foldid=tmtab3[,lno], family=quasipoisson(), intercept=FALSE)
#     cv.inf.spd.ridge        <- cv.glmnet(glmdat.inf.spd,       as.numeric(rslt1[,inf.spd]),       alpha=0, lambda=grid, foldid=rslt1[,lno])
    cv.inf.spd.lasso        <- cv.glmnet(glmdat.inf.spd,       as.numeric(rslt1[,inf.spd]),       alpha=1, lambda=grid, foldid=rslt1[,lno], intercept=FALSE)

    png('ridgelasso.png', width=2000, height=1600)
    par(mfrow=c(2,6))
#     par(mfrow=c(4,6))
#     plot(cv.est.ridge,           main='ridge est')
#     plot(cv.sounder.weeks.ridge, main='ridge sounder.weeks')
#     plot(cv.prop.infd.ridge,     main='ridge prop.infd')
#     plot(cv.max.inc.ridge,       main='ridge max.inc')
#     plot(cv.tm.esc.ridge,        main='ridge tm.esc')
#     plot(cv.inf.spd.ridge,       main='ridge inf.spd')
#
#     plot(cv.est.ridge$glmnet.fit,           main='ridge est')
#     abline(v=-log(cv.est.ridge$lambda.min),           lty=3)
#     abline(v=-log(cv.est.ridge$lambda.1se),           lty=3)
#     plot(cv.sounder.weeks.ridge$glmnet.fit, main='ridge sounder.weeks')
#     abline(v=-log(cv.sounder.weeks.ridge$lambda.min), lty=3)
#     abline(v=-log(cv.sounder.weeks.ridge$lambda.1se), lty=3)
#     plot(cv.prop.infd.ridge$glmnet.fit,     main='ridge prop.infd')
#     abline(v=-log(cv.prop.infd.ridge$lambda.min),     lty=3)
#     abline(v=-log(cv.prop.infd.ridge$lambda.1se),     lty=3)
#     plot(cv.max.inc.ridge$glmnet.fit,       main='ridge max.inc')
#     abline(v=-log(cv.max.inc.ridge$lambda.min),       lty=3)
#     abline(v=-log(cv.max.inc.ridge$lambda.1se),       lty=3)
#     plot(cv.tm.esc.ridge$glmnet.fit,        main='ridge tm.esc')
#     abline(v=-log(cv.tm.esc.ridge$lambda.min),        lty=3)
#     abline(v=-log(cv.tm.esc.ridge$lambda.1se),        lty=3)
#     plot(cv.inf.spd.ridge$glmnet.fit,       main='ridge inf.spd')
#     abline(v=-log(cv.inf.spd.ridge$lambda.min),       lty=3)
#     abline(v=-log(cv.inf.spd.ridge$lambda.1se),       lty=3)

    plot(cv.est.lasso,           main='lasso est')
    plot(cv.sounder.weeks.lasso, main='lasso sounder.weeks')
    plot(cv.prop.infd.lasso,     main='lasso prop.infd')
    plot(cv.max.inc.lasso,       main='lasso max.inc')
    plot(cv.tm.esc.lasso,        main='lasso tm.esc')
    plot(cv.inf.spd.lasso,       main='lasso inf.spd')

    plot(cv.est.lasso$glmnet.fit,           main='lasso est')
    abline(v=-log(cv.est.lasso$lambda.min),           lty=3)
    abline(v=-log(cv.est.lasso$lambda.1se),           lty=3)
    plot(cv.sounder.weeks.lasso$glmnet.fit, main='lasso sounder.weeks')
    abline(v=-log(cv.sounder.weeks.lasso$lambda.min), lty=3)
    abline(v=-log(cv.sounder.weeks.lasso$lambda.1se), lty=3)
    plot(cv.prop.infd.lasso$glmnet.fit,     main='lasso prop.infd')
    abline(v=-log(cv.prop.infd.lasso$lambda.min),     lty=3)
    abline(v=-log(cv.prop.infd.lasso$lambda.1se),     lty=3)
    plot(cv.max.inc.lasso$glmnet.fit,       main='lasso max.inc')
    abline(v=-log(cv.max.inc.lasso$lambda.min),       lty=3)
    abline(v=-log(cv.max.inc.lasso$lambda.1se),       lty=3)
    plot(cv.tm.esc.lasso$glmnet.fit,        main='lasso tm.esc')
    abline(v=-log(cv.tm.esc.lasso$lambda.min),        lty=3)
    abline(v=-log(cv.tm.esc.lasso$lambda.1se),        lty=3)
    plot(cv.inf.spd.lasso$glmnet.fit,       main='lasso inf.spd')
    abline(v=-log(cv.inf.spd.lasso$lambda.min),       lty=3)
    abline(v=-log(cv.inf.spd.lasso$lambda.1se),       lty=3)
    dev.off()

#     coef.est.ridge           <- predict(ridge.mod.est,           s=cv.est.ridge$lambda.min,           type='coefficients', exact=TRUE)
    coef.est.lasso           <- predict(lasso.mod.est,           s=cv.est.lasso$lambda.min,           type='coefficients', exact=TRUE)
#     coef.sounder.weeks.ridge <- predict(ridge.mod.sounder.weeks, s=cv.sounder.weeks.ridge$lambda.min, type='coefficients', exact=TRUE)
    coef.sounder.weeks.lasso <- predict(lasso.mod.sounder.weeks, s=cv.sounder.weeks.lasso$lambda.min, type='coefficients', exact=TRUE)
#     coef.prop.infd.ridge     <- predict(ridge.mod.prop.infd,     s=cv.prop.infd.ridge$lambda.min,     type='coefficients', exact=TRUE)
    coef.prop.infd.lasso     <- predict(lasso.mod.prop.infd,     s=cv.prop.infd.lasso$lambda.min,     type='coefficients', exact=TRUE)
#     coef.max.inc.ridge       <- predict(ridge.mod.max.inc,       s=cv.max.inc.ridge$lambda.min,       type='coefficients', exact=TRUE)
    coef.max.inc.lasso       <- predict(lasso.mod.max.inc,       s=cv.max.inc.lasso$lambda.min,       type='coefficients', exact=TRUE)
#     coef.tm.esc.ridge        <- predict(ridge.mod.tm.esc,        s=cv.tm.esc.ridge$lambda.min,        type='coefficients', exact=TRUE)
    coef.tm.esc.lasso        <- predict(lasso.mod.tm.esc,        s=cv.tm.esc.lasso$lambda.min,        type='coefficients', exact=TRUE)
#     coef.inf.spd.ridge       <- predict(ridge.mod.inf.spd,       s=cv.inf.spd.ridge$lambda.min,       type='coefficients', exact=TRUE)
    coef.inf.spd.lasso       <- predict(lasso.mod.inf.spd,       s=cv.inf.spd.lasso$lambda.min,       type='coefficients', exact=TRUE)

    browser()
    coefs.out <- cbind(#coef.est.ridge,
 coef.est.lasso,
# coef.sounder.weeks.ridge,
 coef.sounder.weeks.lasso ,
# coef.prop.infd.ridge,
 coef.prop.infd.lasso ,
# coef.max.inc.ridge,
 coef.max.inc.lasso ,
# coef.tm.esc.ridge,
 coef.tm.esc.lasso ,
# coef.inf.spd.ridge,
 coef.inf.spd.lasso)
    coefs.out <- data.table(as.matrix(coefs.out))
    s.vals <- unlist(lapply(colnames(coefs.out), function(x) tstrsplit(x, '=', keep=2)))
    coefs.out[, rwnames := rownames(coef.est.lasso)]
    coefs.out <- rbind(coefs.out, t((c(s.vals, 's'))), use.names=FALSE)
    coefs.out <- rbind(coefs.out, t(c(
#         cv.est.ridge$lambda.min,
        cv.est.lasso$lambda.min,
#         cv.sounder.weeks.ridge$lambda.min,
        cv.sounder.weeks.lasso$lambda.min,
#         cv.prop.infd.ridge$lambda.min,
        cv.prop.infd.lasso$lambda.min,
#         cv.max.inc.ridge$lambda.min,
        cv.max.inc.lasso$lambda.min,
#         cv.tm.esc.ridge$lambda.min,
        cv.tm.esc.lasso$lambda.min,
#         cv.inf.spd.ridge$lambda.min,
        cv.inf.spd.lasso$lambda.min,
        'lambda.min'
    )), use.names=FALSE)
    coefs.out <- rbind(coefs.out, t(c(
#         cv.est.ridge$lambda.1se,
        cv.est.lasso$lambda.1se,
#         cv.sounder.weeks.ridge$lambda.1se,
        cv.sounder.weeks.lasso$lambda.1se,
#         cv.prop.infd.ridge$lambda.1se,
        cv.prop.infd.lasso$lambda.1se,
#         cv.max.inc.ridge$lambda.1se,
        cv.max.inc.lasso$lambda.1se,
#         cv.tm.esc.ridge$lambda.1se,
        cv.tm.esc.lasso$lambda.1se,
#         cv.inf.spd.ridge$lambda.1se,
        cv.inf.spd.lasso$lambda.1se,
        'lambda.1se'
    )), use.names=FALSE)
    coefs.out <- rbind(coefs.out, t(c(
#         cv.est.ridge$cvm[cv.est.ridge$index[1]],
        cv.est.lasso$cvm[cv.est.lasso$index[1]],
#         cv.sounder.weeks.ridge$cvm[cv.sounder.weeks.ridge$index[1]],
        cv.sounder.weeks.lasso$cvm[cv.sounder.weeks.lasso$index[1]],
#         cv.prop.infd.ridge$cvm[cv.prop.infd.ridge$index[1]],
        cv.prop.infd.lasso$cvm[cv.prop.infd.lasso$index[1]],
#         cv.max.inc.ridge$cvm[cv.max.inc.ridge$index[1]],
        cv.max.inc.lasso$cvm[cv.max.inc.lasso$index[1]],
#         cv.tm.esc.ridge$cvm[cv.tm.esc.ridge$index[1]],
        cv.tm.esc.lasso$cvm[cv.tm.esc.lasso$index[1]],
#         cv.inf.spd.ridge$cvm[cv.inf.spd.ridge$index[1]],
        cv.inf.spd.lasso$cvm[cv.inf.spd.lasso$index[1]],
        'mse.min'
    )), use.names=FALSE)
    coefs.out <- rbind(coefs.out, t(c(
#         cv.est.ridge$cvm[cv.est.ridge$index[2]],
        cv.est.lasso$cvm[cv.est.lasso$index[2]],
#         cv.sounder.weeks.ridge$cvm[cv.sounder.weeks.ridge$index[2]],
        cv.sounder.weeks.lasso$cvm[cv.sounder.weeks.lasso$index[2]],
#         cv.prop.infd.ridge$cvm[cv.prop.infd.ridge$index[2]],
        cv.prop.infd.lasso$cvm[cv.prop.infd.lasso$index[2]],
#         cv.max.inc.ridge$cvm[cv.max.inc.ridge$index[2]],
        cv.max.inc.lasso$cvm[cv.max.inc.lasso$index[2]],
#         cv.tm.esc.ridge$cvm[cv.tm.esc.ridge$index[2]],
        cv.tm.esc.lasso$cvm[cv.tm.esc.lasso$index[2]],
#         cv.inf.spd.ridge$cvm[cv.inf.spd.ridge$index[2]],
        cv.inf.spd.lasso$cvm[cv.inf.spd.lasso$index[2]],
        'mse.1se'
    )), use.names=FALSE)
#     coefs.out <- rbind(coefs.out, unlist(list(colnames(coefs.out), 's'))
    setnames(coefs.out, c(#'est.ridge',
'est.lasso',
#                           'sounder.weeks.ridge',
'sounder.weeks.lasso',
#                           'prop.infd.ridge',
'prop.infd.lasso',
#                           'max.inc.ridge',
'max.inc.lasso',
#                           'tm.esc.ridge',
'tm.esc.lasso',
#                           'inf.spd.ridge',
'inf.spd.lasso',
'rwnames'))


#     coefs.out <- rbind(coefs.out, c(
#                     cv.est.ridge$cvm[cv.est.ridge$index[1]],
#                     cv.est.lasso$cvm[cv.est.lasso$index[1]],
#                     cv.sounder.weeks.ridge$cvm[cv.sounder.weeks.ridge$index[1]],
#                     cv.sounder.weeks.lasso$cvm[cv.sounder.weeks.lasso$index[1]],
#                     cv.prop.infd.ridge$cvm[cv.prop.infd.ridge$index[1]],
#                     cv.prop.infd.lasso$cvm[cv.prop.infd.lasso$index[1]],
#                     cv.max.inc.ridge$cvm[cv.max.inc.ridge$index[1]],
#                     cv.max.inc.lasso$cvm[cv.max.inc.lasso$index[1]],
#                     cv.tm.esc.ridge$cvm[cv.tm.esc.ridge$index[1]],
#                     cv.tm.esc.lasso$cvm[cv.tm.esc.lasso$index[1]],
#                     cv.inf.spd.ridge$cvm[cv.inf.spd.ridge$index[1]],
#                     cv.inf.spd.lasso$cvm[cv.inf.spd.lasso$index[1]]))
#

#     c(cv.tm.esc.ridge$cvm[cv.tm.esc.ridge$index[1]], cv.tm.esc.lasso$cvm[cv.tm.esc.lasso$index[1]])

#     setnames()

    fwrite(coefs.out, 'coefsout.csv')
#     coefs.out <- coefs.out[, c(2, 4, 6, 8, 10, 11)]

    modlist <- list(#ridge.mod.est,
                    lasso.mod.est,
#                     ridge.mod.sounder.weeks,
                    lasso.mod.sounder.weeks,
#                     ridge.mod.prop.infd,
                    lasso.mod.prop.infd,
#                     ridge.mod.max.inc,
                    lasso.mod.max.inc,
#                     ridge.mod.tm.esc,
                    lasso.mod.tm.esc,
#                     ridge.mod.inf.spd,
                    lasso.mod.inf.spd)

#     sval.list <- c(cv.est.ridge$lambda.min, cv.est.lasso$lambda.min, cv.sounder.weeks.ridge$lambda.min, cv.sounder.weeks.lasso$lambda.min,
#                    cv.prop.infd.ridge$lambda.min, cv.prop.infd.lasso$lambda.min, cv.max.inc.ridge$lambda.min, cv.max.inc.lasso$lambda.min,
#                    cv.tm.esc.ridge$lambda.min, cv.tm.esc.lasso$lambda.min, cv.inf.spd.ridge$lambda.min, cv.inf.spd.lasso$lambda.min)

    return(list(coefs.out, modlist))# sval.list))

}



