
rm(list=ls())



homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))

predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)
names(predict.dat)


# And rescale alpha in both vectors

# Constant
predict.dat$alpha_star_human_constant[!is.na(predict.dat$alpha_star_human_constant)] <- scales::rescale(x =predict.dat$alpha_star_human_constant[!is.na(predict.dat$alpha_star_human_constant)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_constant_lci[!is.na(predict.dat$alpha_star_human_constant_lci)] <- scales::rescale(x =predict.dat$alpha_star_human_constant_lci[!is.na(predict.dat$alpha_star_human_constant_lci)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_constant_uci[!is.na(predict.dat$alpha_star_human_constant_uci)] <- scales::rescale(x =predict.dat$alpha_star_human_constant_uci[!is.na(predict.dat$alpha_star_human_constant_uci)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 


# Complete
predict.dat$alpha_star_human_complete[!is.na(predict.dat$alpha_star_human_complete)] <- scales::rescale(x =predict.dat$alpha_star_human_complete[!is.na(predict.dat$alpha_star_human_complete)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_complete_lci[!is.na(predict.dat$alpha_star_human_complete_lci)] <- scales::rescale(x =predict.dat$alpha_star_human_complete_lci[!is.na(predict.dat$alpha_star_human_complete_lci)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_complete_uci[!is.na(predict.dat$alpha_star_human_complete_uci)] <- scales::rescale(x =predict.dat$alpha_star_human_complete_uci[!is.na(predict.dat$alpha_star_human_complete_uci)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 


#

# And make table S2
TableS2 <- dplyr::select(predict.dat, -(N_Tw), -(N_cumulative), -(phylo_dist), -(N_g0), (N_mu))
head(TableS2)

TableS2[,2:ncol(TableS2)]<- signif(TableS2[,2:ncol(TableS2)], 3)
TableS2$Tw_constant <- paste0(TableS2$Tw_constant, " [", TableS2$Tw_constant_lci, "-", TableS2$Tw_constant_uci, "]")
TableS2$Tw_complete <- paste0(TableS2$Tw_complete, " [", TableS2$Tw_complete_lci, "-", TableS2$Tw_complete_uci, "]")
TableS2$Tv_constant <- TableS2$Tv_human_constant
TableS2$Tv_complete <- TableS2$Tv_human_complete
TableS2$mu <- paste0(TableS2$mu, " [", TableS2$mu_lci, "-", TableS2$mu_uci, "]")
TableS2$g0 <- paste0(TableS2$g0, " [", TableS2$g0_lci, "-", TableS2$g0_uci, "]")
TableS2$rstar_constant <- paste0(TableS2$rstar_constant, " [", TableS2$rstar_constant_lci, "-", TableS2$rstar_constant_uci, "]")
TableS2$rstar_complete <- paste0(TableS2$rstar_complete, " [", TableS2$rstar_complete_lci, "-", TableS2$rstar_complete_uci, "]")
TableS2$alphastar_human_constant <- paste0(TableS2$alpha_star_human_constant, " [", TableS2$alpha_star_human_constant_lci, "-", TableS2$alpha_star_human_constant_uci, "]")
TableS2$alphastar_human_complete <- paste0(TableS2$alpha_star_human_complete, " [", TableS2$alpha_star_human_complete_lci, "-", TableS2$alpha_star_human_complete_uci, "]")

TableS2 <- dplyr::select(TableS2, order, mu,  Tw_constant, Tw_complete, g0, Tv_constant, Tv_complete, rstar_constant, rstar_complete, alphastar_human_constant, alphastar_human_complete)
head(TableS2)

TableS2$order <- as.character(TableS2$order)
TableS2 <- arrange(TableS2, order)
TableS2$order <- as.factor(TableS2$order)

head(TableS2)

#TableS2 <- dcast(melt(TableS2, id.vars = "order"), formula = variable~order)

# And save
write.csv(TableS2, file = paste0(homewd,"supp-figs/TableS2.csv"), row.names = F)
