rm(list=ls())

library(compiler)



source("master.R")

master <- cmpfun(master)

source("initialization.R")

genesis <- cmpfun(genesis)

source("competition")

competition <- cmpfun(competition)

source("diaspora.R")

diaspora <- cmpfun(diaspora)

source("meetic.R")

meetic <- cmpfun(meetic)

source("stillborns.R")

stillborn <- cmpfun(stillborn)

source("plotdensity_asexuals.R")

plotdensity <- cmpfun(plotdensity)

source("plottruc.R")

plottruc <- cmpfun(plottruc)

source("plotGRbabies.R")

plotGR <- cmpfun(plotGR)