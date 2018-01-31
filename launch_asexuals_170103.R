rm(list=ls())

library(compiler)
library(ggplot2)
library(RColorBrewer)


source("master.R")

master <- cmpfun(master)

source("genesis.R")

genesis <- cmpfun(genesis)

source("competition.R")

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