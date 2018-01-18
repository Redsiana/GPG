rm(list=ls())

library(compiler)



source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/master_asexuals_ded_171120.R")

master <- cmpfun(master)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/initialization_asexuals_171120.R")

genesis <- cmpfun(genesis)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/competition_170103.R")

dogeatdog <- cmpfun(dogeatdog)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/faircompet_170103.R")

faircompetition <- cmpfun(faircompetition)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/diaspora_customkernel_asexuals_170103.R")

diaspora <- cmpfun(diaspora)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/meetic_asexuals_171120.R")

meetic <- cmpfun(meetic)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/stillborns_asexuals_171120.R")

stillborn <- cmpfun(stillborn)

source("D:/Kokkonuts/Where sex/modelling/nouvelle geographie/plotdensity_asexuals.R")

plotdensity <- cmpfun(plotdensity)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/plottruc.R")

plottruc <- cmpfun(plottruc)

source("D:/Kokkonuts/Where sex/modelling/hexagonal grid/plotGRbabies.R")

plotGR <- cmpfun(plotGR)


competition <- dogeatdog
# competition <- faircompetition
 