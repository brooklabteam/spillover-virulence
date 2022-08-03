# This script is similar to the one which accompanies the brief -
# But it produces Fig 3 of the main textThis is the script to accompany the brief


rm(list=ls())

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(RColorBrewer)

homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))
