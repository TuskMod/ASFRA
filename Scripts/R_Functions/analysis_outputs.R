
# Analyze the outputs from data.looper (i.e. stats tests, etc.)

# make_rslt <- function(result.outputs, variables, mv.params){
# #     rslt <- data.table(matrix(unlist(result.outputs), ncol=length(result.outputs[[1]]), byrow=TRUE))
# rslt <- fread('result.outputs.csv')
#     setnames(rslt, c('v','l','r','est','edge.tm','max.dist','inf.area','inf.spd','sounder.weeks','prop.infd','max.inc','tm.esc'))
#     rslt[est==2, est:=1]
#     setDT(variables)
#     variables[,v := 1:.N]
#     setkey(variables, v)
#     setkey(mv.params, index)
#     setkey(rslt, v, l, r)
#     rslt <- rslt[variables, on=.(v)]
#     rslt <- rslt[mv.params, on=.(l=index)]
# #     rslt[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('gamma.shape', 'gamma.scale', 'moranI', 'tc', 'mast', 'rgd', 'rds','dayl','prcp','drt','contag','aggindex','entropy','simpindx')]
#     rslt[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('est', 'contact', 'variant', 'density')]
#
#     ## cleared out colinear predictors based on car::vif(), eliminating redundant or semi-redundant predictors one at a time (e.g. contag and aggindex, or entropy and simpindx, reflect similar processes)
#     rslt <- rslt[,.(v,l,r # identifying replicates
#                     # response variables
#                         , est             # establishment (0, 1)
# #                         , edge.tm         # time until first exposure in edge cell (<= 78 weeks) or end of simulation
# #                         , max.dist        # max distance attained by exposure (up to 50.25, distance to closest edge, if reached)
# #                         , inf.area        # total number of cells with positive exposures
#                         , inf.spd         # speed of infection wave (average change in average distance of sounders with E|I|C > 0 by week)
#                         , sounder.weeks   # sounders with E|I>0 across time (i.e. sum all at t=1, t=2, t=3, etc.)
# #                         , infpot.cells    # count of all cells within the maximum distance of disease spread (infected and clear)
#                         , prop.infd       # proportion of the potentially infected cells (dist <= max.dist) that have actually had E|I>0
#                         , max.inc         # maximum incidence (E_t - E_t-1)
#                         , tm.esc          # time to where 5% of the sounders are outside the 10km radius
#                     # epi/population variables
#                         , contact         # (b1, b2, f1, f2_int, f2_b, f2i_int, f2i_b),
#                         , variant         # DR, Pol
#                         , density         # 1.5, 5 pigs/km^2
#                     # landscape variables
#                         , nnd_med         # median nearest neighbor value of sounder dispersal on landscape according to habitat preference
# #                         , nnd_range       # range of nnd values as above
#                         , nnd_mean
# #                         , nnd_sd
#                         , nnd_cv
# #                         , gamma.shape     # movement gamma distribution shape parameter; these should be treated as one parameter (maybe shape*scale to give mean?)
# #                         , gamma.scale     # movement gamma distribution scale parameter; could also do coefficient of variation?
#                         , mvmt.mean
#                         , mvmt.cv
#                         , moranI          # Moran's I statistic of spatial autocorrelation of habitat preference value
#                         , tc              # tree cover index
#                         , mast            # masting species
#                         , rgd             # ruggedness measure (?)
#                         , rds             # roads index
#                         , dayl            # daylight
# #                         , prcp            # precipitation
# #                         , tmin
#                         , tmax
# #                         , drt             # drought index
# #                         , contag          # contagion of land cover classes
# #                         , aggindex        # aggregation index of land cover classes
# #                         , entropy         # marginal entropy of land cover classes
#                         , simpindx        # simpson's diversity index of land cover classes
#                         )]
#     return(rslt)
# }
#
#
#
# analysis.output <- function(rslt, variables, mv.params){

#     rslt <- rbindlist(rslt) # for good measure
# library(data.table)
# library(targets)
# library(car)
# library(lme4)
# library(glmmTMB)

# library(BART)
# rslt <- fread('result.outputs.csv')
#     rslt[est==2, est:=1]
# variables <- tar_read(variables)
# mv.params <- tar_read(mv.params)
#
#     setDT(variables)
#     variables[,v := 1:.N]
#     setkey(variables, v)
#     setkey(mv.params, index)
#     setkey(rslt, v, l, r)

#     rslt <- rslt[variables, on=.(v)]
#     rslt <- rslt[mv.params, on=.(l=index)]
#     rslt[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('gamma.shape', 'gamma.scale', 'moranI', 'tc', 'mast', 'rgd', 'rds','dayl','prcp','drt','contag','aggindex','entropy','simpindx')]
#     rslt[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('est', 'contact', 'variant', 'density')]
    ## check colinearlity
#     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd_med + nnd_range + nnd_mean + nnd_sd + gamma.shape + gamma.scale + moranI + gearyC + tc + mast + rgd + rds + dayl + prcp + tmin + tmax + drt + contag + aggindex + entropy + simpindx + nnd_cv + mvmt.mean + mvmt.cv + st.wt + (1|l), data=rslt, family = poisson, control=glmerControl(autoscale=TRUE)))
    ## to document the parameter colinearity reduction
#     >     car::vif(glmer(inf.area ~ contact + variant + density + nnd_med + nnd_ra$
#                          contact      variant      density      nnd_med    nnd_range     nnd_mean       nnd_sd  gamma.shape  gamma.scale       moranI       gearyC           tc         mast          rgd          rds         dayl         prcp         tmin         tmax          drt
#                          1.000000     1.000000     1.000000    12.121197     3.404680   103.525497   115.952679   152.098523   303.126406 49966.019562 49820.535170    14.682812    24.059521     7.364861     5.780511     1.677864    22.009235    59.374986    68.629495    12.075931
#                          contag     aggindex      entropy     simpindx       nnd_cv    mvmt.mean      mvmt.cv        st.wt
#                          2350.564705   119.044290  1371.856335    55.872202    84.237482   163.513747   399.533759     7.484814
#                          > sort(    car::vif(glmer(inf.area ~ contact + variant + density + nnd_med + n$
#
#                                                    > vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd_med$
#                                                                                  > vif.preds
#                                                                                  contact      variant      density      nnd_med    nnd_range     nnd_mean       nnd_sd  gamma.shape  gamma.scale       moranI       gearyC           tc         mast          rgd          rds         dayl         prcp         tmin         tmax          drt
#                                                                                  1.000000     1.000000     1.000000    12.121197     3.404680   103.525497   115.952679   152.098523   303.126406 49966.019562 49820.535170    14.682812    24.059521     7.364861     5.780511     1.677864    22.009235    59.374986    68.629495    12.075931
#                                                                                  contag     aggindex      entropy     simpindx       nnd_cv    mvmt.mean      mvmt.cv        st.wt
#                                                                                  2350.564705   119.044290  1371.856335    55.872202    84.237482   163.513747   399.533759     7.484814
#                                                                                  > sort(vif.preds)
#                                                                                  contact      variant      density         dayl    nnd_range          rds          rgd        st.wt          drt      nnd_med           tc         prcp         mast     simpindx         tmin         tmax       nnd_cv     nnd_mean       nnd_sd     aggindex
#                                                                                  1.000000     1.000000     1.000000     1.677864     3.404680     5.780511     7.364861     7.484814    12.075931    12.121197    14.682812    22.009235    24.059521    55.872202    59.374986    68.629495    84.237482   103.525497   115.952679   119.044290
#                                                                                  gamma.shape    mvmt.mean  gamma.scale      mvmt.cv      entropy       contag       gearyC       moranI
#                                                                                  152.098523   163.513747   303.126406   399.533759  1371.856335  2350.564705 49820.535170 49966.019562
#                                                                                  >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                    > sort(vif.preds)
#                                                                                                                    contact     variant     density        dayl   nnd_range         rds       st.wt      gearyC         rgd     nnd_med         drt          tc        prcp        mast    simpindx        tmin        tmax      nnd_cv    nnd_mean    aggindex      nnd_sd
#                                                                                                                    1.000000    1.000000    1.000000    1.645454    3.354041    4.426941    5.785051    6.107677    7.008249   12.022526   12.058202   14.532550   21.848368   23.779556   55.822585   59.049248   65.910857   81.267949  101.305794  109.207038  113.151316
#                                                                                                                    gamma.shape   mvmt.mean gamma.scale     mvmt.cv     entropy      contag
#                                                                                                                    150.880418  159.997288  290.647891  398.895259 1212.711581 2079.450800
#                                                                                                                    >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                      > vif.preds
#                                                                                                                                                      contact    variant    density    nnd_med  nnd_range   nnd_mean     nnd_sd     moranI         tc       mast        rgd        rds       dayl       prcp       tmin       tmax        drt     contag   simpindx     nnd_cv  mvmt.mean    mvmt.cv      st.wt
#                                                                                                                                                      1.000000   1.000000   1.000000   9.359084   3.021223  96.956603 109.297635   3.791586  10.804425  16.390466   5.452375   3.876169   1.570789  18.966716  40.505177  53.300962  10.315643  46.510157  42.342030  78.512529   8.012636   3.673606   4.866551
#                                                                                                                                                      > sort(vif.preds)
#                                                                                                                                                      contact    variant    density       dayl  nnd_range    mvmt.cv     moranI        rds      st.wt        rgd  mvmt.mean    nnd_med        drt         tc       mast       prcp       tmin   simpindx     contag       tmax     nnd_cv   nnd_mean     nnd_sd
#                                                                                                                                                      1.000000   1.000000   1.000000   1.570789   3.021223   3.673606   3.791586   3.876169   4.866551   5.452375   8.012636   9.359084  10.315643  10.804425  16.390466  18.966716  40.505177  42.342030  46.510157  53.300962  78.512529  96.956603 109.297635
#                                                                                                                                                      >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                                                        > sort(vif.preds)
#                                                                                                                                                                                        contact   variant   density      dayl nnd_range    nnd_cv   mvmt.cv    moranI       rds     st.wt       rgd mvmt.mean       drt   nnd_med  nnd_mean        tc      mast      prcp      tmin  simpindx    contag      tmax
#                                                                                                                                                                                        1.000000  1.000000  1.000000  1.568812  2.831891  3.089148  3.647844  3.787138  3.868417  4.866193  5.146342  7.711748  8.800868  9.328942 10.043508 10.747726 16.355950 18.643617 35.410553 41.325208 43.964622 47.373656
#                                                                                                                                                                                        >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                                                                                          > vif.preds
#                                                                                                                                                                                                                          contact   variant   density   nnd_med nnd_range  nnd_mean    moranI        tc      mast       rgd       rds      dayl      prcp      tmin      tmax       drt  simpindx    nnd_cv mvmt.mean   mvmt.cv     st.wt
#                                                                                                                                                                                                                          1.000000  1.000000  1.000000  9.248048  2.752404 10.003631  3.468402 10.712805 15.180907  5.118402  3.756139  1.568309 18.512439 33.047269 44.751016  8.390026  4.720711  2.920490  5.639246  3.624411  4.666068
#                                                                                                                                                                                                                          > sort(vif.preds)
#                                                                                                                                                                                                                          contact   variant   density      dayl nnd_range    nnd_cv    moranI   mvmt.cv       rds     st.wt  simpindx       rgd mvmt.mean       drt   nnd_med  nnd_mean        tc      mast      prcp      tmin      tmax
#                                                                                                                                                                                                                          1.000000  1.000000  1.000000  1.568309  2.752404  2.920490  3.468402  3.624411  3.756139  4.666068  4.720711  5.118402  5.639246  8.390026  9.248048 10.003631 10.712805 15.180907 18.512439 33.047269 44.751016
#                                                                                                                                                                                                                          >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                                                                                                                            > vif.preds
#                                                                                                                                                                                                                                                            contact   variant   density   nnd_med  nnd_mean    moranI        tc      mast       rgd       rds      dayl      prcp      tmax       drt  simpindx    nnd_cv mvmt.mean   mvmt.cv     st.wt
#                                                                                                                                                                                                                                                            1.000000  1.000000  1.000000  8.613450  7.401882  3.245992  8.179982 13.681745  3.506283  3.671831  1.481985 12.914078  6.297713  7.786109  4.672815  1.847182  5.595587  3.438926  3.407144
#                                                                                                                                                                                                                                                            > sort(vif.preds)
#                                                                                                                                                                                                                                                            contact   variant   density      dayl    nnd_cv    moranI     st.wt   mvmt.cv       rgd       rds  simpindx mvmt.mean      tmax  nnd_mean       drt        tc   nnd_med      prcp      mast
#                                                                                                                                                                                                                                                            1.000000  1.000000  1.000000  1.481985  1.847182  3.245992  3.407144  3.438926  3.506283  3.671831  4.672815  5.595587  6.297713  7.401882  7.786109  8.179982  8.613450 12.914078 13.681745
#                                                                                                                                                                                                                                                            >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                                                                                                                                 > sort(vif.preds)
#                                                                                                                                                                                                                                                                 contact   variant   density      dayl    nnd_cv    moranI     st.wt   mvmt.cv       rds       rgd  simpindx        tc mvmt.mean      tmax  nnd_mean       drt   nnd_med      mast
#                                                                                                                                                                                                                                                                 1.000000  1.000000  1.000000  1.481985  1.846112  2.672872  3.205743  3.334964  3.390232  3.505754  3.767774  5.468735  5.521903  6.191158  6.880466  7.449204  8.310320 13.556838
#                                                                                                                                                                                                                                                                 >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                                                                                                                                 > sort(vif.preds)
#                                                                                                                                                                                                                                                                 contact   variant   density      dayl    nnd_cv       rds     st.wt   mvmt.cv    moranI       rgd  simpindx mvmt.mean      tmax  nnd_mean        tc   nnd_med      mast      prcp
#                                                                                                                                                                                                                                                                 1.000000  1.000000  1.000000  1.453563  1.834521  2.857555  2.914208  3.121008  3.245027  3.467743  4.256479  4.605902  5.693778  7.067865  7.884551  8.517221  9.063288 12.355284
#                                                                                                                                                                                                                                                                 >     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd$
#                                                                                                                                                                                                                                                                 > sort(vif.preds)
#                                                                                                                                                                                                                                                                 contact   variant   density      dayl    nnd_cv    moranI       rds     st.wt   mvmt.cv       rgd  simpindx mvmt.mean        tc      tmax  nnd_mean   nnd_med      mast
#                                                                                                                                                                                                                                                                 1.000000  1.000000  1.000000  1.452262  1.834431  2.656171  2.734619  2.818417  2.918911  3.463481  3.558602  4.601552  5.426497  5.444779  6.689233  8.270296  8.393573
#     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd_med + nnd_mean + moranI + tc + mast + rgd + rds + dayl + tmax + simpindx + nnd_cv + mvmt.mean + mvmt.cv + (1|l), data=rslt, family = poisson, control=glmerControl(autoscale=TRUE)))

#     rslt <- rslt[,.(v,l,r # identifying replicates
#                     # response variables
#                         , est             # establishment (0, 1)
# #                         , edge.tm         # time until first exposure in edge cell (<= 78 weeks) or end of simulation
# #                         , max.dist        # max distance attained by exposure (up to 50.25, distance to closest edge, if reached)
# #                         , inf.area        # total number of cells with positive exposures
#                         , inf.spd         # speed of infection wave (average change in average distance of sounders with E|I|C > 0 by week)
#                         , sounder.weeks   # sounders with E|I>0 across time (i.e. sum all at t=1, t=2, t=3, etc.)
# #                         , infpot.cells    # count of all cells within the maximum distance of disease spread (infected and clear)
#                         , prop.infd       # proportion of the potentially infected cells (dist <= max.dist) that have actually had E|I>0
#                         , max.inc         # maximum incidence (E_t - E_t-1)
#                         , tm.esc          # time to where 5% of the sounders are outside the 10km radius
#                     # epi/population variables
#                         , contact         # (b1, b2, f1, f2_int, f2_b, f2i_int, f2i_b),
#                         , variant         # DR, Pol
#                         , density         # 1.5, 5 pigs/km^2
#                     # landscape variables
#                         , nnd_med         # median nearest neighbor value of sounder dispersal on landscape according to habitat preference
# #                         , nnd_range       # range of nnd values as above
#                         , nnd_mean
# #                         , nnd_sd
#                         , nnd_cv
# #                         , gamma.shape     # movement gamma distribution shape parameter; these should be treated as one parameter (maybe shape*scale to give mean?)
# #                         , gamma.scale     # movement gamma distribution scale parameter; could also do coefficient of variation?
#                         , mvmt.mean
#                         , mvmt.cv
#                         , moranI          # Moran's I statistic of spatial autocorrelation of habitat preference value
#                         , tc              # tree cover index
#                         , mast            # masting species
#                         , rgd             # ruggedness measure (?)
#                         , rds             # roads index
#                         , dayl            # daylight
# #                         , prcp            # precipitation
# #                         , tmin
#                         , tmax
# #                         , drt             # drought index
# #                         , contag          # contagion of land cover classes
# #                         , aggindex        # aggregation index of land cover classes
# #                         , entropy         # marginal entropy of land cover classes
#                         , simpindx        # simpson's diversity index of land cover classes
#                         )]
#
#
#
#
#
#     ## use only the data with est==1 (fizzler) or est==2 (established)
#     rslt <- rslt[est==2, est:=1]
#     rslt1 <- rslt[est==1,]


#     rslt2 <- rslt[est==2,]
#     rslt2.preds <- rslt2[,.(v, l, r, contact, variant, density, nnd_med, nnd_range, gamma.shape, gamma.scale, moranI, tc, mast, rgd, rds, dayl, prcp, drt)]


    ## INFERENCE
    # ESTABLISHMENT -- definitely this one
        # if max infections <= 1 at each time step, and total exposures overall are <= 10, est = 0 (non established)
        # else if time to edge is <78 and maximum distance is <50, est = 1 (established, but stopped before reaching edge/end time, "fizzled out")
        # else est = 2 (reched the edge before 78 weeks, OR lasted for the entire 78 weeks)
    # establishment is probably the result of contact, density, and maybe variant (for longevity allowing the first few transmissions)
    # establishmet probably doesn't have anything to do with habitat fracturing (because the threshold of establishment is so low), roads, masting,
    # tc, and probably ruggedness except at the most extreme cases (very bad environments will have minimal establishment, but those are probably not in our datasets)
    # it may have more to do with seasonality, if we included that (carcass degradation, movement, daylight hours, temperatures, etc.)
    # I think we can test this with the cvd core and nnd parameters, as well as movement parameters
#     contact is the intermediary between density and variant -- variant controls duration of opportunity, density controls frequency of opportunity, but contact controsl success of each opportunity
    est.glmer.dredge <- function(rslt){

        est.glmer <- glmmTMB(est ~
                contact + variant + density +    # epi/population parameters (Hi/Lo)
                contact:density + contact:variant +
                nnd_med + nnd_mean + nnd_cv +    # sounder distribution parameters
                mvmt.mean + mvmt.cv +            # sounder movement (mean, coefficint of variance)
                moranI +                         # landscape preference autocorrelation
                rgd +                            # landscape ruggedness
                tc + mast +                      # vegetation (tree cover, masting species)
                rds +                            # roads
#                 dayl +                         # turns out average year-long daylight averages have a range of about 20 seconds-- should remove that one from everything...
                tmax +
                simpindx +             # land cover classification fragmentation metrics
                (1|l),                           # random effects -- landscape ID
            data = rslt # includes non-established
            , family=binomial
            , na.action = na.fail
#             , control = glmerControl(autoscale=TRUE)
            , control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000), profile=TRUE)
    )
    # what should drive establishment? obviously contact, density, variant, but also processes upstream from that --
    # - how far away sounders tend to be (nnd_med)
    # - how much sounders move (gamma shape and scale)
    # -
    library(MuMIn)
    est.dredge <- dredge(est.glmer, trace=2,
                                evaluate=FALSE,
                         fixed=c('cond(contact)','cond(variant)','cond(density)','cond(contact:variant)','cond(contact:density)'))
    return(est.dredge)
    }
#     est.avg <- model.avg(est.dredge)
#     est.avg <- model.avg(subset(est.dredge, cumsum(weight) <= 0.95))


    # INFECTION SPEED -- super important... use that one
    # what would affect infection speed? most things... things related to contact, pig movement, including variances/cv... probably not climate
    # connectivity nnd, sounder movement, autocorrelation of habitat preference? maybe!
    # ruggedness probably yeah, tree cover I suspect has to do with how clumpy the sounders are on the environment
    inf.spd.glmer.dredge <- function(rslt1){
        inf.spd.lmer <- glmmTMB(inf.spd ~
                contact + variant + density +    # epi/population parameters (Hi/Lo)
                contact:density + contact:variant +
                nnd_med + nnd_mean + nnd_cv +    # sounder distribution parameters
                mvmt.mean + mvmt.cv +            # sounder movement (mean, coefficint of variance)
                moranI +                         # landscape preference autocorrelation
                rgd +                            # landscape ruggedness
                tc + mast +                      # vegetation (tree cover, masting species)
                rds +                            # roads
#                 dayl +                         # turns out average year-long daylight averages have a range of about 20 seconds-- should remove that one from everything...
                tmax +
                simpindx +             # land cover classification fragmentation metrics
                (1|l),                           # random effects -- landscape ID
            data = rslt1 # includes those who established and didn't necessarily make it to the 50km mark
            ,family=gaussian
#             , control = lmerControl(autoscale=TRUE)
            , control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000), profile=TRUE)
            , na.action = na.fail
        )
        inf.spd.dredge <- dredge(inf.spd.lmer, trace=2,
                                evaluate=FALSE,
                         fixed=c('cond(contact)','cond(variant)','cond(density)','cond(contact:variant)','cond(contact:density)'))
        return(inf.spd.dredge)
    }
#     inf.spd.avg <- model.avg(subset(inf.spd.dredge, cumsum(weight) <= 0.95), fit=TRUE)

    # SOUNDER WEEKS
    # should be a factor of spread, longevity, ... distribution on the landscape
    sounder.weeks.glmer.dredge <- function(rslt1){
        sounder.weeks.glmer <- glmmTMB(sounder.weeks ~
                contact + variant + density +    # epi/population parameters (Hi/Lo)
                contact:density + contact:variant +
                nnd_med + nnd_mean + nnd_cv +    # sounder distribution parameters
                mvmt.mean + mvmt.cv +            # sounder movement (mean, coefficint of variance)
                moranI +                         # landscape preference autocorrelation
                rgd +                            # landscape ruggedness
                tc + mast +                      # vegetation (tree cover, masting species)
                rds +                            # roads
#                 dayl +                         # turns out average year-long daylight averages have a range of about 20 seconds-- should remove that one from everything...
                tmax +
                simpindx +             # land cover classification fragmentation metrics
                (1|l),                           # random effects -- landscape ID
                family=nbinom2,
            data = rslt1 # includes those who didn't make it to the 50km mark
            , na.action = na.fail
            , control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000), profile=TRUE)
        )
        snd.weeks.dredge <- dredge(sounder.weeks.glmer, trace=2,
                                evaluate=FALSE,
                                fixed = c('cond(contact)','cond(variant)','cond(density)','cond(contact:variant)','cond(contact:density)'))
        return(snd.weeks.dredge)
    }
#     snd.weeks.avg <- model.avg(subset(snd.weeks.dredge, cumsum(weight) <= 0.95), fit=TRUE)

    # PROPORTION INFECTED CELLS -- makes sense as a beta
    # infected proportion of cells would have to do with connectivity (nnd, and fragmentation), movement and contact
    prop.infd.beta.dredge <- function(rslt1){
        rslt1[prop.infd == 1, prop.infd := 0.999]
        prop.infd.beta <- glmmTMB::glmmTMB(prop.infd ~
                contact + variant + density +    # epi/population parameters (Hi/Lo)
                contact:density + contact:variant +
                nnd_med + nnd_mean + nnd_cv +    # sounder distribution parameters
                mvmt.mean + mvmt.cv +            # sounder movement (mean, coefficint of variance)
                moranI +                         # landscape preference autocorrelation
                rgd +                            # landscape ruggedness
                tc + mast +                      # vegetation (tree cover, masting species)
                rds +                            # roads
#                 dayl +                         # turns out average year-long daylight averages have a range of about 20 seconds-- should remove that one from everything...
                tmax +
                simpindx +             # land cover classification fragmentation metrics
                (1|l),                           # random effects -- landscape ID
                family=glmmTMB::beta_family(link='logit'),
            data = rslt1 # includes those who didn't make it to the 50km mark
            , na.action = na.fail
            , control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000), profile=TRUE)
        )
        prop.infd.dredge <- dredge(prop.infd.beta, trace=2,
                                evaluate=FALSE,
                                fixed = c('cond(contact)','cond(variant)','cond(density)','cond(contact:variant)','cond(contact:density)'))
        return(prop.infd.dredge)
    }
#     prop.infd.avg <- model.avg(subset(prop.infd.dredge, cumsum(weight) <= 0.95), fit=TRUE)
#     summary(prop.infd.beta)

    # MAXIMUM INCIDENCE
    # should be density, movement, contact, continuity, variant, etc.
    max.inc.glmer.dredge <- function(rslt1){
        max.inc.glmer.nbinom <- glmmTMB(max.inc ~
                contact + variant + density +    # epi/population parameters (Hi/Lo)
                contact:density + contact:variant +
                nnd_med + nnd_mean + nnd_cv +    # sounder distribution parameters
                mvmt.mean + mvmt.cv +            # sounder movement (mean, coefficint of variance)
                moranI +                         # landscape preference autocorrelation
                rgd +                            # landscape ruggedness
                tc + mast +                      # vegetation (tree cover, masting species)
                rds +                            # roads
#                 dayl +                         # turns out average year-long daylight averages have a range of about 20 seconds-- should remove that one from everything...
                tmax +
                simpindx +             # land cover classification fragmentation metrics
                (1|l),                           # random effects -- landscape ID
            family=nbinom2,
            data = rslt1 # includes those who didn't make it to the 50km mark
            , na.action = na.fail
            , control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000), profile=TRUE)
        )
        max.inc.dredge <- dredge(max.inc.glmer.nbinom, trace=2,
                                evaluate=FALSE,
                                fixed = c('cond(contact)','cond(variant)','cond(density)','cond(contact:variant)','cond(contact:density)'))
        return(max.inc.dredge)
    }
#     max.inc.avg <- model.avg(subset(max.inc.dredge, cumsum(weight) <= 0.95), fit=TRUE)

    # TIME TO ESCAPE
    # contact, density, connectivity,
    # reset data to event status at each time step
    tm.tab.3 <- function(rslt1, variables, mv.params){
        library(data.table)
        setDT(variables)
        variables[,v := 1:.N]
        setDT(mv.params)
        tmtab <- CJ(v=1:8, l=unique(rslt1[,l]), r=1:100, tm=1:78)
        tmtab2 <- tmtab[rslt1[,.(v,l,r,tm.esc)], on=.(v,l,r)]
        tmtab3 <- tmtab2[tm <= tm.esc | is.na(tm.esc),]
        tmtab3[, exc:=0]
        tmtab3[tm == tm.esc, exc := 1]
        tmtab3 <- tmtab3[variables[,.(v, contact, variant, density)], on=.(v)]
        tmtab3 <- tmtab3[mv.params, on=.(l=index)]
        tmtab3[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact','variant','density')]
    }

    tm.esc.binom.dredge <- function(tmtab3, variables, mv.params){
        tm.esc.bin <- glmmTMB(exc ~ tm +
                contact + variant + density +    # epi/population parameters (Hi/Lo)
                contact:density + contact:variant +
                nnd_med + nnd_mean + nnd_cv +    # sounder distribution parameters
                mvmt.mean + mvmt.cv +            # sounder movement (mean, coefficint of variance)
                moranI +                         # landscape preference autocorrelation
                rgd +                            # landscape ruggedness
                tc + mast +                      # vegetation (tree cover, masting species)
                rds +                            # roads
#                 dayl +                         # turns out average year-long daylight averages have a range of about 20 seconds-- should remove that one from everything...
                tmax +
                simpindx +             # land cover classification fragmentation metrics
                (1|l),                           # random effects -- landscape ID
            family=binomial(link = 'logit'),
            data = tmtab3 # includes those who didn't make it to the 50km mark
            , na.action = na.fail
    #          , control = glmerControl(autoscale=TRUE)
            , control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000), profile=TRUE)
        )
        tm.esc.dredge <- dredge(tm.esc.bin, trace=2,
                                evaluate=FALSE,
                                fixed = c('cond(contact)','cond(variant)','cond(density)','cond(contact:variant)','cond(contact:density)','cond(tm)'))
        return(tm.esc.dredge)
    }
#     tm.esc.avg <- model.avg(subset(tm.esc.dredge, cumsum(weight) <= 0.95), fit=TRUE)

    model.eval <- function(dredge.list, rslt) {
        # depending on which dredge.list is coming in (from which response variable) any of rslt, rslt1, or tmtab3 may be required
        mod <- eval(dredge.list[[1]])
        return(mod)
    }

    model.eval1 <- function(dredge.list, rslt1) {
        # depending on which dredge.list is coming in (from which response variable) any of rslt, rslt1, or tmtab3 may be required
        rslt1[prop.infd == 1, prop.infd := 0.999]
        mod <- eval(dredge.list[[1]])
        return(mod)
    }

    model.eval.esc <- function(dredge.list, tmtab3) {
        # depending on which dredge.list is coming in (from which response variable) any of rslt, rslt1, or tmtab3 may be required
        mod <- eval(dredge.list[[1]])
        return(mod)
    }

    ref.to.list <- function(modoutput, is.esc){
        if (is.esc == 1){
            modoutput <- modoutput[[1]]
        }
        formula.resp <- paste(formula(modoutput)[2])
        formula.pred <- paste(formula(modoutput)[3])
        nvars <- extractAIC(modoutput)[1]
        AIC <- extractAIC(modoutput)[2]
        nll <- logLik(modoutput)[1]
        coefs <- apply(coef(modoutput)[[1]]$l, 2, mean)
        return(c(list(resp = formula.resp, preds = formula.pred, nvars = nvars, AIC = AIC, nll = nll), coefs))
    }

    list.to.table <- function(modlist){
        dat <- rbindlist(modlist, fill=TRUE)
        dat[,id := 1:.N]
        setorder(dat, 'AIC')
        return(dat[, rank := 1:.N][,dAIC := AIC - min(AIC)])
    }

    table.filter <- function(dat, daicval, modelmin){
        dat.out <- dat[dAIC<=daicval,]
        if (nrow(dat.out) <= modelmin){
            dat.out <- dat[order(dAIC)][1:modelmin,]
        }
        return(dat.out)
    }


    model.avgs <- function(x.table, x.dredge){
        modids <- x.table[,id]
        avgmod <- model.avg(x.dredge, fit=TRUE)
        coefs <- data.table(t(coef(avgmod)))
        nms <- lapply(names(coefs), function(x) gsub('cond','', x))
        nms <- lapply(nms, function(x) gsub('\\(','', x))
        nms <- lapply(nms, function(x) gsub('\\)','', x))
        names(coefs) <- unlist(nms)
        setnames(coefs, 'Int','(Intercept)')
        outtab <- rbind(x.table, coefs, fill=TRUE)
        return(outtab)
    }

    combine.models <- function(est.model.avg, inf.spd.model.avg, snd.wk.model.avg, prop.infd.model.avg, max.inc.model.avg, tm.esc.model.avg){
        est.model.avg[is.na(resp),resp := 'est.avg']
        inf.spd.model.avg[is.na(resp),resp := 'inf.spd.avg']
        snd.wk.model.avg[is.na(resp),resp := 'sounder.weeks.avg']
        prop.infd.model.avg[is.na(resp),resp := 'prop.infd.avg']
        max.inc.model.avg[is.na(resp),resp := 'max.inc.avg']
        tm.esc.model.avg[is.na(resp),resp := 'exc.avg']
        try(setnames(tm.esc.model.avg, 'contactLo:density', 'contactLo:density5'), silent=TRUE)
        mod.avgs <- rbindlist(list(est.model.avg, inf.spd.model.avg, snd.wk.model.avg, prop.infd.model.avg, max.inc.model.avg, tm.esc.model.avg), fill=TRUE)
        #??? add names of models and column for loo.land
    }

    lolo.table <- function(mv.params, allmods, rslt){
        lolo.table <- CJ(lindx=unique(rslt[,l]), resp.var=c('est','inf.spd','sounder.weeks','prop.infd','max.inc','exc'), fam='none')#, est.resp=0, exc.resp=0)
        ## connect this table with output of combine.models; add exp.var based on allmods results and lnd based on the names of lands
        lolo.table[resp.var == 'est', fam := 'binomial']
#         lolo.table[resp.var == 'est', est.resp := 1]
        lolo.table[resp.var == 'inf.spd', fam := 'gaussian']
        lolo.table[resp.var == 'sounder.weeks', fam := 'nbinom2']
        lolo.table[resp.var == 'prop.infd', fam := 'beta_family']
        lolo.table[resp.var == 'max.inc', fam := 'nbinom2']
        lolo.table[resp.var == 'exc', fam := 'binomial']
#         lolo.table[resp.var == 'exc', exc.resp := 1]
        lolo.table[,resp.avg := paste0(resp.var, '.avg')]
#         lolo.table <- allmods[lolo.table, on=.(resp = resp.avg)]
        lolo.table <- allmods[lolo.table, on=.(resp = resp.var), allow.cartesian=TRUE]
        return(lolo.table)
    }
    # resp, preds, nvars, AIC, nll, [coefficients], id, rank, dAIC, lindx, fam, resp.avg

    lolo.predict <- function(lolo.table.pred, rslt, rslt1, tmtab3){
        lindx <- lolo.table.pred[,lindx]
        resp.var <- lolo.table.pred[,resp]
        fam <- lolo.table.pred[,fam]
#         coefs <- 'contact + variant + density + contact:density + contact:variant + nnd_med + nnd_mean + nnd_cv + mvmt.mean + mvmt.cv + moranI + rgd + tc + mast + rds + tmax + simpindx + (1|l)'
        coefs <- lolo.table.pred[,preds]
#         tm.in <- ''
        dat <- rslt
        if (resp.var == 'exc'){
            dat <- tmtab3
        } else if (resp.var != 'est') {
            dat <- rslt1
        }
        if (resp.var == 'prop.infd') {
            rslt1[prop.infd == 1, prop.infd := 0.999]
            dat <- rslt1
        }
        trn <- dat[l != lindx,]
        tst <- dat[l == lindx,]
        ## will need to change how models are run in glmmTMB based on structure of mod.list
        m <- glmmTMB(formula(paste0(resp.var, '~', coefs)),
                     data=trn, family=fam,
                     control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))#, profile=TRUE))
        tst[,yhat := predict(m, newdata=tst, re.form=NA, type='response', allow.new.levels = TRUE)]
        tst[,resp := resp.var][,mod := coefs]
        return(tst)
    }
    # v, l, r, est, inf.spd, sounder.weeks, prop.infd, max.inc, tm.esc, contact, variant, density, [exp.vars], resp.var

    calc.rmse <- function(test.fit.tab, lolo.combos){
        #might need to add model numbers to the tables before they are combined
#         tab[,mod := .GRP, by=.(l, resp.var, mod)] # labels groups by number according to the by= statement
        test.fit.tab <- test.fit.tab[unique(lolo.combos[,.(resp, fam)]), on=.NATURAL]
        test.fit.tab[,ydiff.sq := (as.numeric(.SD[[.BY[[1]]]]) - yhat)^2, by=resp] # pulls the values from resp.var column and uses them to choose which columns (named like the values) to use for .SD
        rmse.tab <- unique(test.fit.tab[,rmse := sqrt(sum(ydiff.sq)/.N), by=.(l, resp, mod)][,.(l, resp, mod, rmse, fam)])
        return(rmse.tab)
    }

    get.pred.mods <- function(rmse.tab, allmods, rslt, rslt1, tmtab3){
        avg.mods <- rmse.tab[,mean(rmse), by=.(resp, mod, fam)]
        best.mods <- avg.mods[order(V1)][,.SD[1], by=resp]
        best.mods.tab <- allmods[best.mods, on=.(resp,preds=mod)]
        setnames(best.mods.tab, 'V1', 'rmse.avg')
        mod.list <- lapply(1:nrow(best.mods.tab), function(x){
            resp.var <- best.mods.tab[x, resp]
            fam <- best.mods.tab[x, fam]
            coefs <- best.mods.tab[x, preds]
            dat <- rslt
            if (resp.var == 'exc'){
                dat <- tmtab3
            } else if (resp.var != 'est') {
                dat <- rslt1
            }
            if (resp.var == 'prop.infd') {
                rslt1[prop.infd == 1, prop.infd := 0.999]
                dat <- rslt1
            }
            m <- glmmTMB(formula(paste0(resp.var, '~', coefs)),
                     data=dat, family=fam,
                     control=glmmTMBControl(optCtrl=list(iter.max=1000,eval.max=1000)))
            return(m)
        })
        return(mod.list)
    }

oos.tiles <- function(tiledat, rslt, variables){
    # the rest of the tiles' data for predictions
    tiledat <- fread(tiledat)
#     tiledat <- fread('./Landscape_Setup/NND_Lands/all_tile_attribs.csv')
    setnames(tiledat, c('index', 'nnd_med', 'nnd_range', 'nnd_mean', 'nnd_sd', 'x', 'y', 'state', 'gamma.shape', 'gamma.scale',
                        'moranI', 'gearyC', 'tc', 'mast', 'rgd', 'rds', 'dayl', 'prcp', 'tmin', 'tmax', 'drt', 'contag', 'aggindex',
                        'entropy', 'simpindx', 'nnd_cv', 'mvmt.mean', 'mvmt.cv'))
#     tiledat <- tiledat[index %in% rslt[,l] == FALSE, ]
    setDT(variables)
    variables[,v := 1:.N]
    tvar <- CJ(l = unique(tiledat[,index]), var = variables[,v])
    tiledat <- tvar[tiledat, on=.(l=index)][variables, on=.(var=v)]

#     tiledat <- tiledat[,.(l, var, contact, variant, density, nnd_med, nnd_range, nnd_mean, nnd_sd, nnd_cv, gamma.shape, gamma.scale, moranI, tc, mast, rgd, rds, dayl, prcp, drt, contag, aggindex, entropy, simpindx, mvmt.mean, mvmt.cv, tmax)]
    tiledat[,names(.SD) := lapply(.SD, as.numeric), .SDcols=c('gamma.shape', 'gamma.scale', 'moranI', 'tc', 'mast', 'rgd', 'rds','dayl','prcp','drt','contag','aggindex','entropy','simpindx')]
    tiledat[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('contact', 'variant', 'density')]

    # non-correlated covariates from above
#     vif.preds <- car::vif(glmer(inf.area ~ contact + variant + density + nnd_med + nnd_mean + moranI + tc + mast + rgd + rds + dayl + tmax + simpindx + nnd_cv + mvmt.mean + mvmt.cv + (1|l), data=rslt, family = poisson, control=glmerControl(autoscale=TRUE)))
#     tiledat <- tiledat[,.(l, var, contact, variant, density, nnd_med, nnd_mean, nnd_cv, moranI, tc, mast, rgd, rds, dayl, tmax, simpindx, mvmt.mean, mvmt.cv)]
    return(tiledat)
}


oos.preds <- function(predmods, lambda.mins, tiledat){
    resp.var <- str_trim(unlist(tstrsplit(unlist(tstrsplit(paste(predmods[[1]]$call[3]), ',', keep=2)), ']', keep=1)))
    if (resp.var == 'exc') {
        tiledat <- tiledat[CJ(l=unique(tiledat[,l]), tm=1:78), on=.(l), allow.cartesian=TRUE]
    }
    tiledat.trim <- tiledat[,.(contact, variant, density, nnd_med, nnd_range, nnd_mean, nnd_sd, moranI, gearyC, tc, mast, rgd, rds, dayl, prcp, tmin, tmax, drt, contag, aggindex, entropy, simpindx, nnd_cv, mvmt.mean, mvmt.cv)]
    pred <- predict(predmods[[1]],
                    newx=model.matrix(~., tiledat.trim)[,-1],
                    s=lambda.mins)#exact=TRUE)
    coefs.vals <- coef(predmods[[1]], s=lambda.mins)
    coefs.used <- rownames(coefs.vals)[which(coefs.vals != 0)]
    # add some kind of identifier to the output for the response variable, maybe the model too
    dimnames(pred) <- NULL
    pred.table <- cbind(tiledat, resp.variable=resp.var, coefs.used=paste(coefs.used, collapse=' '), lambda=lambda.mins, pred)
    setnames(pred.table, 'V1', 'pred')
#     setnames(pred.table, c('V2','V3'), c('resp.variable','preds'))
#     if (paste(formula(predmods[[1]])[2]) == 'exc') {
#         pred.table <- pred.table[pred >= 0.5, .SD[1], by=.(l, var)]
#     }
    return(pred.table)
}

preds.compile <- function(preds, preds_edge, tiledat, tiledat_edge){
    preds <- rbindlist(preds, fill=TRUE, use.names=TRUE)
    preds_edge <- rbindlist(preds_edge, fill=TRUE, use.names=TRUE)
    setnames(preds_edge, 'l', 'index')
    # dcast to wide-ish format
#     preds.out <- dcast(preds[is.na(tm) | tm == 52,], l + var ~ resp.variable, value.var=c('pred'), fun.aggregate = mean)
    preds.out <- dcast(preds, l + var ~ resp.variable, value.var=c('pred'), fun.aggregate = mean)
    preds.out.edge <- dcast(preds_edge, index + var ~ resp.variable, value.var=c('pred'), fun.aggregate = mean)
#     preds.out.edge <- dcast(preds_edge[is.na(tm) | tm == 52,], index + var ~ resp.variable, value.var=c('pred'), fun.aggregate = mean)
    preds.out <- preds.out[tiledat, on=.(l, var)]
#     coef.table <- unique(preds[,.(resp.variable, coefs.used)])
    preds.out[,edg := 0]
    preds.out.edge <- preds.out.edge[tiledat_edge, on=.(index=l, var)]
    preds.out.edge[,edg := 1]
    preds.out <- rbind(preds.out, preds.out.edge, fill=TRUE)
    return(preds.out)
#     return(list(preds.out, coef.table))
}




maps.plot <- function(preds.table, rslt, rslt1, tmtab3, variables){
    # national maps of response variables
    # not looped by targets
    library(sf)
    library(raster)
    library(terra)

    # grab the selected lands shapefile if it exists (not the first run of this function)
    sel.files <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)
    if (file.exists('./Input/lands_bound/sel_bounds.shp')){
        selpts.box <- st_read('./Input/lands_bound/sel_bounds.shp')
#         same.tiles <- all(selpts.box$l %in% unique(rslt[,l]))
    } else {
        # creates the polygons outlining the borders of the landscape tiles
        # this is slow, so only runs if necessary
        selpts.box <- bind(lapply(sel.files, function(tif){
            box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
            return(box)
        }))
        # get the names of the selected files
        name <- tstrsplit(sel.files, '_')
        name <- name[length(name)]
        name <- unlist(tstrsplit(unlist(name), '\\.', keep=1))
        name <- as.numeric(name)
        selpts.box <- SpatialPolygonsDataFrame(selpts.box, data=data.table(l=name))
        # create the output shapefile
        raster::shapefile(x=selpts.box, file='./Input/lands_bound/sel_bounds.shp')
    }

    # same as above, but for all tiles
    if (file.exists('./Input/lands_bound/all_bounds.shp')){
        pot.pts.box <- st_read('./Input/lands_bound/all_bounds.shp')
    } else {
        tile.potential <- list.files('./Landscape_Setup/Pipeline_SSF_Weekly/4_Output/indiv_plands', full.names=TRUE)
        pot.pts.box <- bind(lapply(tile.potential, function(tif){
            name <- tstrsplit(tif, '_')
            name <- name[length(name)]
            name <- paste0('tile', unlist(tstrsplit(unlist(name), '\\.', keep=1)))
            box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
            names(box) <- name
            return(box)
        }))
        name <- tstrsplit(tile.potential, '_')
        name <- name[length(name)]
        name <- unlist(tstrsplit(unlist(name), '\\.', keep=1))
        name <- as.numeric(name)
        pot.pts.box <- SpatialPolygonsDataFrame(pot.pts.box, data=data.table(l=name))
        raster::shapefile(x=pot.pts.box, file='./Input/lands_bound/all_bounds.shp')
    }

    # same as above but for edge (non-candidate) tiles
    if (file.exists('./Input/lands_bound/edge_bounds.shp')){
        edgepts.box <- st_read('./Input/lands_bound/edge_bounds.shp')
#         same.tiles <- all(selpts.box$l %in% unique(rslt[,l]) == FALSE)
    } else {
        edge.files <- list.files('./Landscape_Setup/Pipeline_SSF_Weekly/4_Output/edge_plands', full.names=TRUE)
        edgepts.box <- bind(lapply(edge.files, function(tif){
            box <- rasterToPolygons(reclassify(raster(tif), matrix(c(0,1,1), ncol=3)), dissolve=TRUE)[1]
            return(box)
        }))
        name <- tstrsplit(edge.files, '_')
        name <- name[length(name)-1]
        name <- as.numeric(unlist(name))
        edgepts.box <- SpatialPolygonsDataFrame(edgepts.box, data=data.table(l=name))
        raster::shapefile(x=edgepts.box, file='./Input/lands_bound/edge_bounds.shp')
    }

    # get the contiguous us state border shapefile
    contus <- st_read('./Input/contus_shp/') # for windows must be directory, might be something else for linux?
    # grab the geographic data from a landscape tile to make sure everything matches
    tile.one <- list.files('./Landscape_Setup/NND_Lands/4_Output/sel_plands', full.names=TRUE)[1]
    rast.crs <- crs(raster(tile.one))
    contus.transform <- st_transform(contus[1], rast.crs)

    # create bounding box of contiguous US
    bounds <- st_bbox(contus.transform)
    xrng <- bounds[c(1,3)]
    yrng <- bounds[c(2,4)]

    # loop through response variables
    lapply(c('est','inf.spd','sounder.weeks','prop.infd','max.inc','exc'), function(pred.var){
        # get the values predicted by the models to set the ranges necessary for the map colors
        pred.var.vals <- values(vect(merge(pot.pts.box, preds.table[,.SD,.SDcols=c('l','var',pred.var)], by='l')))
        pred.var.vals$edge <- 0
        pred.var.vals.edge <- values(vect(merge(edgepts.box, preds.table[,.SD,.SDcols=c('l','var',pred.var)], by='l')))
        pred.var.vals.edge$edge <- 1
        pred.var.vals.all <- rbind(pred.var.vals, pred.var.vals.edge)
        max.val <- max(pred.var.vals.all[, 3])
        min.val <- min(pred.var.vals.all[, 3])
        # scaled response variables from 0 to 1
        pixels <- (pred.var.vals.all[,3]-min.val)/(max.val-min.val)
        pred.var.vals.all$pixels <- pixels
        # define color gradient
        color.grads <- colorRampPalette(c('blue','purple','red','orange','yellow'))
        # match colors to response variable
        rank <- as.factor( as.numeric( cut(pixels, 150)))
        pred.var.vals.all$rank <- as.numeric(as.character(rank))
        pred.var.vals <- pred.var.vals.all[pred.var.vals.all$edge == 0,]
        pred.var.vals.edge <- pred.var.vals.all[pred.var.vals.all$edge == 1,]
        setDT(variables)
        variables[,v:=1:.N]
        # plot
        png(paste0('./Output/figures/lm_tile_map_', pred.var, '.png'), width=2000, height=1000)
        par(oma=c(0.2,0.2,3,0.2))
        layout(matrix(c(seq(1,8),rep(9,4)), nrow=3, byrow=TRUE), widths=c(1,1,1,1), heights=c(1,1,0.3))
        lapply(c(3,1,4,2,7,5,8,6), function(x){
            variable.line <- unlist(variables[v==x, 1:3])
            #     connect tiledat and prediction data to specific tiles
            selpts.box.sub <- vect(merge(selpts.box, preds.table[var==x, .SD, .SDcols=c('l', pred.var)], by='l'))
            # escape time needs its own form, mostly just flipping the colors because lower = worse case
            if (pred.var %in% c('tm.esc')) {
                terra::plot(vect(merge(pot.pts.box, pred.var.vals[pred.var.vals$var == x,], by='l')), col=rev(color.grads(150))[pred.var.vals[pred.var.vals$var == x,'rank']], pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n', border=NA)
                terra::plot(vect(merge(edgepts.box, pred.var.vals.edge[pred.var.vals.edge$var==x,], by='l')), col=rev(color.grads(150))[pred.var.vals.edge[pred.var.vals.edge$var == x,'rank']], pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n', border=NA, add=TRUE)
#                 mtext(paste(pred.var, 'vars', x), 3, 0, cex=2)
                mtext(paste('contact', variable.line[1], '| variant', variable.line[2], '| density', variable.line[3]), 3, -1.5, cex=1.5)
                plot(contus.transform, col=NA, lwd=1.5, alpha=0, fill=NA, border='darkgrey', add=TRUE)
            } else {
                terra::plot(vect(merge(pot.pts.box, pred.var.vals[pred.var.vals$var == x,], by='l')), col=color.grads(150)[pred.var.vals[pred.var.vals$var == x,'rank']], pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n', border=NA)
                terra::plot(vect(merge(edgepts.box, pred.var.vals.edge[pred.var.vals.edge$var == x,], by='l')), col=color.grads(150)[pred.var.vals.edge[pred.var.vals.edge$var == x, 'rank']], pch=1, cex=1.3, xlim=xrng, ylim=yrng, ann=FALSE, asp=1, xaxt='n', yaxt='n', bty='n', border=NA, add=TRUE)
                mtext(paste('contact', variable.line[1], '| variant', variable.line[2], '| density', variable.line[3]), 3, -1.5, cex=1.5)
                plot(contus.transform, col=NA, lwd=1.5, alpha=0, fill=NA, border='darkgrey', add=TRUE)
            }
        })
        # add a legend
        mtext(resp.translate(pred.var), 3, 0, cex=2, outer=TRUE)
        legendimg <- as.raster(matrix(color.grads(150), nrow=1))
        plot(c(0,1), c(0,1), type='n', axes=F, xlab = '', ylab='', main=resp.translate(pred.var), cex.main=1.5)
        min.val2 <- floor(log10(abs(min.val)))
        mtext(round(seq(min.val, max.val, length.out=5), -min.val2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
        rasterImage(legendimg, 0, 0, 1, 1)
        dev.off()
    })
}


## heatmaps
    heatmap <- function(mod, oos.tile.dat, oos.tile.dat.edge, lambda.mins){
        # looped by targets inputs on mod and lambda.mins
        # use the 2 most important landscape attributes, vary them for the axes, and for each combination of contact, variant, and density make a colored heatmap of response value predictions
        # mark points where actual landscapes were
        # use averages of other landscape values for model inputs

        # pull in model coefficients
        coefs <- coef(mod[[1]], s=lambda.mins)
        # keep only the landscape-associated values
        coefs <- coefs[(rownames(coefs) %in% c('(Intercept)', 'contactLo', 'density5', 'variantPol')==FALSE), 1]
        if (length(coefs[coefs!=0]) >=2){
            # grab the two coefficients with the greatest effect sizes for the plot axes
            ## make sure these leave negative coefficients as negative
            coefs <- sort(abs(unlist(coefs)), decreasing=TRUE)[1:2]
            # grab coefficient names
            cnames <- names(coefs)
            cf1nm <- cnames[1]
            cf2nm <- cnames[2]
            # get the tile data for the non-simulated landscapes
            tiledats <- rbind(oos.tile.dat, oos.tile.dat.edge, fill=TRUE, use.names=TRUE)
            # get ranges of coefficients for each
            coef1.range <- range(tiledats[,..cf1nm])
            coef2.range <- range(tiledats[,..cf2nm])
            # create a table with all combinations of the active coefficients and the non-landscape variables
            ht.table <- CJ(
                cf1 = seq(coef1.range[1], coef1.range[2], length.out=25)
                , cf2 = seq(coef2.range[1], coef2.range[2], length.out=25)
                , density = c(1.5,5)
                , contact = c('Hi','Lo')
                , variant = c('DR','Pol'))
            # keep the landscape variables from the tile data and average them
            landval.avgs <- apply(oos.tile.dat[,c(3:6,12:29)], 2, mean)
            # filter out landscape variables that are included in the plot axes
            landval.avgs <- landval.avgs[names(landval.avgs) %in% names(coefs) == FALSE]
            # get the name of the response variable
            resp.var <- str_trim(unlist(tstrsplit(unlist(tstrsplit(paste(mod[[1]]$call[3]), ',', keep=2)), ']', keep=1)))
            setnames(ht.table, c('cf1','cf2'), c(cnames[1], cnames[2]))
            # connect coefficient combination with average landscape variables
            ht.table <- cbind(ht.table, as.data.table(t(landval.avgs)))
            browser()
            # set population/epidemiology variables as factors
            ht.table[,names(.SD) := lapply(.SD, as.factor), .SDcols=c('density','contact','variant')]
            # predict the response for each variable value combination
            ht.table[,pred := predict(mod[[1]], newx=model.matrix(~., ht.table)[,-1], s=lambda.mins)]#, type='response', allow.new.levels=TRUE)]
            # make table of variaable combinations to feed into mapply function below (to generate 8 heatmaps)
            quick.vars <- unique(ht.table[,.(density, contact, variant)])
            # define color ramp
            color.grads <- colorRampPalette(c('blue','purple','red','orange','yellow'))
            # match the color gradient with response variable values
            ht.table[,rank := as.factor( as.numeric( cut(pred, 150)))]
            ht.table[,colr := color.grads(150)[as.numeric(as.character(rank))]]
            # plot
            png(paste0('./Output/figures/heat_test_', cf1nm, '_', cf2nm,'_', resp.var, '.png'), width=1700, height=1050)
            par(mfrow=c(2,4), oma=c(1,1,4,1))
            layout(matrix(c(1:9,9,9,9), nrow=3, ncol=4, byrow=TRUE), widths=c(1,1,1,1), heights=c(1,1,0.3))
            mapply(function(dens, cont, varnt){
                rows <- ht.table[,which(density == dens & contact == cont & variant == varnt)]
    #             sub.rank <- rank[rows]
                ht.sub <- ht.table[rows, ]
                plot(unlist(ht.sub[,1]), unlist(ht.sub[,2]), col=ht.sub[,colr], pch=15, cex=4,ann=FALSE, pty='s')
                mtext(paste('density:', dens,'| contact:', cont, '| variant:', varnt), 3, 1.2, cex=1.4)
                mtext(exp.translate(cf1nm), 1, 2)
                mtext(exp.translate(cf2nm), 2, 2)
                tile.pts <- tiledats[variant == varnt & contact == cont & density == dens,.SD, .SDcols=c(cf1nm, cf2nm)]
                points(tile.pts)
            }, dens=quick.vars[,density], cont=quick.vars[,contact], var=quick.vars[,variant])
            mtext(resp.translate(resp.var), 3, 2, outer=TRUE, cex=2)
            # create gradient legend object
            legendimg <- as.raster(matrix(color.grads(150), nrow=1))
            plot(c(0,1), c(0,1), type='n', axes=F, xlab = '', ylab='', main=resp.translate(resp.var), cex.main=1.5)
            min.val2 <- floor(log10(abs(min(ht.table[,pred]))))
            mtext(round(seq(range(ht.table[,pred])[1], range(ht.table[,pred])[2], length.out=5), -min.val2), 1, 1, at= c(0,0.25,0.5,0.75,1), cex=1.4)
            rasterImage(legendimg, 0, 0, 1, 1)
            dev.off()
        } else {
            # if there aren't two or more landscape attributes that are kept by the lasso regression, no heatmap is possible
            print('this isn\'t going to work')
        }
    return(NA)
    }



resp.translate <- function(rv){
    # translates response variable codes into text for labeling figures
    print(rv)
    if(rv == 'est') rv.out <- 'Epidemic Establishment Proportion'
    if(rv == 'inf.spd') rv.out <- 'Epidemic Wave Max Speed (km/wk)'
    if(rv == 'sounder.weeks') rv.out <- 'Epidemic Intensity (Sounder-Weeks)'
    if(rv == 'prop.infd') rv.out <- 'Proportion of Cells Infected'
    if(rv == 'max.inc') rv.out <- 'Maximum Incidence'
    if(rv == 'exc') rv.out <- 'Proportion of Epidemics Escaped by 52 Weeks'
    if(rv == 'tm.esc') rv.out <- 'Proportion of Epidemics Escaped by 52 Weeks'
    return(rv.out)
}

exp.translate <- function(exp.var){
    # translates explanatory variable codes into text for labeling figures
    if(exp.var == 'nnd_med') ev.out <- 'Median Nearest Neighbor Distance'
    if(exp.var == 'nnd_range') ev.out <- 'Range Nearest Neighbor Distance'
    if(exp.var == 'nnd_mean') ev.out <- 'Mean Nearest Neighbor Distance'
    if(exp.var == 'nnd_cv') ev.out <- 'Nearest Neighbor Coefficient of Variation'
    if(exp.var == 'nnd_sd') ev.out <- 'Nearest Neighbor Standard Deviation'
    if(exp.var == 'mvmt.mean') ev.out <- 'Mean Sounder Movement'
    if(exp.var == 'mvmt.cv') ev.out <- 'Sounder Movement Coefficient of Variation'
    if(exp.var == 'moranI') ev.out <- 'Landscape Preference Autocorrelation (Moran\'s I)'
    if(exp.var == 'gearyC') ev.out <- 'Landscape Preference Autocorrelation (Geary\'s C)'
    if(exp.var == 'tc') ev.out <- 'Tree Cover'
    if(exp.var == 'mast') ev.out <- 'No. Masting Species'
    if(exp.var == 'rgd') ev.out <- 'Landscape Ruggedness'
    if(exp.var == 'rds') ev.out <- 'Roads Index'
    if(exp.var == 'dayl') ev.out <- 'Daylight'
    if(exp.var == 'prcp') ev.out <- 'Precip'
    if(exp.var == 'tmax') ev.out <- 'Mean Annual Maximum Daily Temperature'
    if(exp.var == 'tmin') ev.out <- 'Mean Annual Minimum Daily Temperature'
    if(exp.var == 'simpindx') ev.out <- 'Landscape Cover Simpson Diversity Index'
    if(exp.var == 'drt') ev.out <- 'Drought'
    if(exp.var == 'contag') ev.out <- 'Contagion'
    if(exp.var == 'aggindex') ev.out <- 'aggregaiton index'
    if(exp.var == 'entropy') ev.out <- 'entropy'
    return(ev.out)
}
