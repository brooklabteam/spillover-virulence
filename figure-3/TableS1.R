
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

head(predict.dat)
# express mu in scientific notation
predict.dat$mu <- scales:: scientific(predict.dat$mu, digits = 3)
predict.dat$mu_lci <- scales:: scientific(predict.dat$mu_lci, digits = 3)
predict.dat$mu_uci <- scales:: scientific(predict.dat$mu_uci, digits = 3)

# And make table S1
TableS1 <- dplyr::select(predict.dat, -(N_Tw), -(N_cumulative), -(phylo_dist), -(N_g0), (N_mu))
head(TableS1)

TableS1[,5:ncol(TableS1)]<- signif(TableS1[,5:ncol(TableS1)], 3)
TableS1$Tw_constant <- paste0(TableS1$Tw_constant, " [", TableS1$Tw_constant_lci, "-", TableS1$Tw_constant_uci, "]")
TableS1$Tw_complete <- paste0(TableS1$Tw_complete, " [", TableS1$Tw_complete_lci, "-", TableS1$Tw_complete_uci, "]")
TableS1$Tv_constant <- TableS1$Tv_human_constant
TableS1$Tv_complete <- TableS1$Tv_human_complete
TableS1$mu <- paste0(TableS1$mu, " [", TableS1$mu_lci, "-", TableS1$mu_uci, "]")
TableS1$g0 <- paste0(TableS1$g0, " [", TableS1$g0_lci, "-", TableS1$g0_uci, "]")
TableS1$rstar_constant <- paste0(TableS1$rstar_constant, " [", TableS1$rstar_constant_lci, "-", TableS1$rstar_constant_uci, "]")
TableS1$rstar_complete <- paste0(TableS1$rstar_complete, " [", TableS1$rstar_complete_lci, "-", TableS1$rstar_complete_uci, "]")
TableS1$alphastar_human_constant <- paste0(TableS1$alpha_star_human_constant, " [", TableS1$alpha_star_human_constant_lci, "-", TableS1$alpha_star_human_constant_uci, "]")
TableS1$alphastar_human_complete <- paste0(TableS1$alpha_star_human_complete, " [", TableS1$alpha_star_human_complete_lci, "-", TableS1$alpha_star_human_complete_uci, "]")

TableS1 <- dplyr::select(TableS1, order, mu,  Tw_constant, Tw_complete, g0, Tv_constant, Tv_complete, rstar_constant, rstar_complete, alphastar_human_constant, alphastar_human_complete)
head(TableS1)

TableS1$order <- as.character(TableS1$order)

TableS1 <- dplyr::arrange(TableS1, order)
TableS1$order <- as.factor(TableS1$order)

head(TableS1)

#TableS1 <- dcast(melt(TableS1, id.vars = "order"), formula = variable~order)

# And save
write.csv(TableS1, file = paste0(homewd,"supp-figs/TableS1.csv"), row.names = F)
