
# This script produces Figure 3 of the main text.

rm(list=ls())

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(ggimage)

homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))


######################################################
######################################################

#load parameters from GAM analysis
predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)

# Now the literature to compare against
load(paste0(homewd,"source/gam.dat.Guth.et.al.2021.Rdata"))

# run this script to assemble them into a comparison

# assemble a function that takes in the three parameters of interest 
# (host mortality rate, tolerance of immunopathology, extent of constitutive
# immunity) and produces an estimate of relative virulence


b= .2
q = .0002
c=.5
m=1/(21)
g1=.9
#g0 = .3
zeta=.2
v=1
w=1
Tv=.005
# Tw=.005
# mu = 1/(20*365))

dat1 = predict.dat[1,]
dat1 <- dplyr::select(dat1, order, mu, mu_lci, mu_uci, N_mu, g0, g0_lci, g0_uci,  N_g0,Tw_constant, Tw_complete,  phylo_dist, Tv_human_constant, Tv_human_complete)
dat1 <- dplyr::select(predict.dat, order, mu, mu_lci, mu_uci,  N_mu, g0,  g0_lci, g0_uci, N_g0,Tw_constant, Tw_complete,  phylo_dist, Tv_human_constant, Tv_human_complete)
head(dat1)

Tw_guess = dat1$Tw_constant
Tw_joke = rep(1.01, length(Tw_guess))
Tw_high = rep(100, length(Tw_guess))

tolerance_type = "constant"


dat <- dplyr::select(dat1, order, mu,  mu_lci, mu_uci, N_mu, g0, g0_lci, g0_uci,  N_g0, phylo_dist, Tv_human_constant, Tv_human_complete)

dat

#Tw_guess is a vector of guesses in the same order as the orders of the data

make.alpha.prediction <- function(Tw, b, q, c, m, g1, zeta, v, w, Tv, dat, tolerance_type){
  
 
  vir.par <- list(b= b, q = q, c=c,  m=m,  g1=g1,  zeta=zeta,
                  v=v, w=w,  Tv=Tv)
  
  dat$Tw = Tw
  

  
  # Constant tolerance rstar predictions using order-specific g0, Tw, and mu
  # and default values for all other parameters
  if(tolerance_type=="constant"){
    dat$rstar <- (vir.par$c*dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g1*dat$g0*dat$mu*dat$Tw*vir.par$Tv*(vir.par$v*dat$Tw+vir.par$g1*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*dat$Tw) +(vir.par$g1*vir.par$w*(vir.par$m^2)*vir.par$Tv))  

    } else if(tolerance_type=="complete"){
    dat$rstar <- (vir.par$c*dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g1*dat$g0*dat$mu)/(sqrt((vir.par$m^2)*(vir.par$c^2)*dat$mu*vir.par$g1*dat$g0*(vir.par$g1*vir.par$w+vir.par$v-vir.par$g1*dat$Tw-vir.par$Tv)))

    }
  

  
  ######################################################
  ######################################################
  ## 5. Now, estimate alpha* in spillover human host
  
  # Virulence in the human host is a function of phylogenetic distance,
  # represented in order-specific Tvs (shown here as Tv_human_constant)
  # Tolerance of immunopathology is now human, so we hold constant across 
  # all reservoir host orders.
  
  # First, calculate Vmax in the human host: (this is the same whether constant or complete tolerance)
  dat$Vs_max <- dat$rstar/(vir.par$g1*vir.par$c) - dat$rstar/(2*vir.par$g1*vir.par$c) + 1 - 1/(vir.par$g1) + vir.par$c/(2*dat$rstar*vir.par$g1)

  # Next, take that viral load and calculate spillover virulence, here for the constant tolerance assumption:
  if(tolerance_type=="constant"){
  vir.par$Tw = 1
  dat$alpha_star_human <- ((dat$rstar*vir.par$v)/dat$Tv_human_constant  + (vir.par$g1*vir.par$w*dat$rstar)/vir.par$Tw)*dat$Vs_max

  } else if(tolerance_type=="complete"){
    vir.par$Tw = 0
    dat$alpha_star_human <- (dat$rstar*(vir.par$v-dat$Tv_human_complete) + dat$rstar_complete*vir.par$g1*(vir.par$w-vir.par$Tw))*dat$Vs_max

  }
  
  #return, and compare with the data
  return(dat)
}
predict.compare.RSS <- function(par, b, q, c, m, g1, zeta, v, w, Tv, dat, gam.dat, tolerance_type){
  
  #split dataset by order and extract the corresponding human virulence prediction appropriate 
  dat.split <- dlply(dat,.(order))
  

  print(as.list(exp(par)))
  #apply prediction function to data
  out.df <- mapply(FUN = make.alpha.prediction, Tw = as.list(exp(par)), dat=dat.split, MoreArgs = list(b= b, q = q, c=c,  m=m,  g1=g1, zeta=zeta, v=v, w=w,  Tv=Tv, tolerance_type=tolerance_type), SIMPLIFY = FALSE)
  
  #and get the  lci and uci
  dat.lci <- dat.uci <- dat
  dat.lci$mu <- dat.lci$mu_lci
  dat.uci$mu <- dat.uci$mu_uci
  
  dat.lci$g0 <- dat.lci$g0_lci
  dat.uci$g0 <- dat.uci$g0_uci
  
  dat.split.lci <- dlply(dat.lci,.(order))
  dat.split.uci <- dlply(dat.uci,.(order))
  
  #out.df.lci <- mapply(FUN = make.alpha.prediction, Tw = as.list(exp(par)), dat=dat.split.lci, MoreArgs = list(b= b, q = q, c=c,  m=m,  g1=g1, zeta=zeta, v=v, w=w,  Tv=Tv, tolerance_type=tolerance_type), SIMPLIFY = FALSE)
  #out.df.uci <- mapply(FUN = make.alpha.prediction, Tw = as.list(exp(par)), dat=dat.split.uci, MoreArgs = list(b= b, q = q, c=c,  m=m,  g1=g1, zeta=zeta, v=v, w=w,  Tv=Tv, tolerance_type=tolerance_type), SIMPLIFY = FALSE)
  
  out.df <- data.table::rbindlist(out.df)
  #out.df.lci <- data.table::rbindlist(out.df.lci)
  #out.df.uci <- data.table::rbindlist(out.df.uci)
  #head(out.df)
  
  #add uci/lci predictions
  #out.df.lci <- dplyr::select(out.df.lci, order, alpha_star_human)
  #out.df.uci <- dplyr::select(out.df.uci, order, alpha_star_human)
  #names(out.df.lci)[names(out.df.lci)=="alpha_star_human"] <- "alpha_star_human_lci"
  #names(out.df.uci)[names(out.df.uci)=="alpha_star_human"] <- "alpha_star_human_uci"
  
  #merge
  #out.df <- merge(out.df, out.df.lci, by="order")
  #out.df <- merge(out.df, out.df.uci, by="order")
  #head(out.df)
  #now compare with the data from the literature
  
  
  
  #rank by descending alpha, then confidence
  gam.dat <- arrange(gam.dat, desc(alpha), desc(Nobs))
  
  
  # Merge our nested model predictions with those from Guth et al. 2022
  
  # Make columns match
  #head(gam.dat)
  names(gam.dat)[1] <- "order"
  names(gam.dat)[3] <- "N"
  
  
  # And rescale alpha 
  #out.df$alpha_star_human[!is.na(out.df$alpha_star_human)] <- scales::rescale(x =out.df$alpha_star_human[!is.na(out.df$alpha_star_human)], from=c(min(out.df$alpha_star_human_lci, na.rm = T), max(out.df$alpha_star_human_uci, na.rm = T)), to =c(0,1)) 
  out.df$alpha_star_human[!is.na(out.df$alpha_star_human)] <- scales::rescale(x =out.df$alpha_star_human[!is.na(out.df$alpha_star_human)], from=c(min(out.df$alpha_star_human, na.rm = T), max(out.df$alpha_star_human, na.rm = T)), to =c(0,1)) 
  
  # And do the same for gam.dat
  #gam.dat$alpha <- scales::rescale(x =gam.dat$alpha, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
  gam.dat$alpha <- scales::rescale(x =gam.dat$alpha, from=c(min(gam.dat$alpha), max(gam.dat$alpha)), to =c(0,1)) 
  gam.dat <- select(gam.dat, -(alpha_lci), -(alpha_uci))
  
  
  # And merge the data
  gam.dat <-  dplyr::select(gam.dat, order, alpha)
  share.dat <- dplyr::select(out.df,order, alpha_star_human)
  names(share.dat)[names(share.dat)=="alpha_star_human"] <- "alpha_from_nested_model"
  
  merge.dat <- merge(gam.dat, share.dat, by="order")
  #with(merge.dat, plot(alpha, alpha_from_nested_model))
  
  merge.dat$sq_diff <- (merge.dat$alpha-merge.dat$alpha_from_nested_model)^2
  
  RSS = sum(merge.dat$sq_diff)
  
  return(RSS)  
}
  


wrap.fit <- function(Tw_guess, b, q, c, m, g1, zeta, v, w, Tv, dat, tolerance_type, gam.dat){
  
  #load your dataset and optimize
  predict.compare.RSS(par=log(Tw_guess), b= b, q = q, c=c,  m=m,  g1=g1,  zeta=zeta,
                      v=v, w=w,  Tv=Tv, dat=dat, tolerance_type=tolerance_type, gam.dat=gam.dat)
  
  predict.compare.RSS(par =log(Tw_joke), b= b, q = q, c=c,  m=m,  g1=g1,  zeta=zeta,
                      v=v, w=w,  Tv=Tv, dat=dat, tolerance_type=tolerance_type, gam.dat=gam.dat)

  
  predict.compare.RSS(par =log(Tw_high), b= b, q = q, c=c,  m=m,  g1=g1,  zeta=zeta,
                      v=v, w=w,  Tv=Tv, dat=dat, tolerance_type=tolerance_type, gam.dat=gam.dat)
  
  out.optim <- optim(par = log(Tw_guess), fn=predict.compare.RSS, method = "Nelder-Mead", b= b, q = q, c=c,  m=m,  g1=g1,  zeta=zeta,
                  v=v, w=w,  Tv=Tv, dat=dat, tolerance_type=tolerance_type, gam.dat=gam.dat, hessian = T)
  
  
  
  
}
  


assemble.dat <- function(predict.dat, gam.dat){
  
  #rank by descending alpha, then confidence
  gam.dat <- arrange(gam.dat, desc(alpha), desc(Nobs))
  
  
  # Merge our nested model predictions with those from Guth et al. 2022
  
  # Make columns match
  head(gam.dat)
  names(gam.dat)[1] <- "order"
  names(gam.dat)[3] <- "N"
  names(predict.dat)[names(predict.dat)=="N_cumulative"] <- "N"
  head(predict.dat)
  
  # And rescale alpha in both vectors
  
  # Constant
  predict.dat$alpha_star_human_constant[!is.na(predict.dat$alpha_star_human_constant)] <- scales::rescale(x =predict.dat$alpha_star_human_constant[!is.na(predict.dat$alpha_star_human_constant)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
  predict.dat$alpha_star_human_constant_lci[!is.na(predict.dat$alpha_star_human_constant_lci)] <- scales::rescale(x =predict.dat$alpha_star_human_constant_lci[!is.na(predict.dat$alpha_star_human_constant_lci)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
  predict.dat$alpha_star_human_constant_uci[!is.na(predict.dat$alpha_star_human_constant_uci)] <- scales::rescale(x =predict.dat$alpha_star_human_constant_uci[!is.na(predict.dat$alpha_star_human_constant_uci)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
  
  
  # Complete
  predict.dat$alpha_star_human_complete[!is.na(predict.dat$alpha_star_human_complete)] <- scales::rescale(x =predict.dat$alpha_star_human_complete[!is.na(predict.dat$alpha_star_human_complete)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
  predict.dat$alpha_star_human_complete_lci[!is.na(predict.dat$alpha_star_human_complete_lci)] <- scales::rescale(x =predict.dat$alpha_star_human_complete_lci[!is.na(predict.dat$alpha_star_human_complete_lci)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
  predict.dat$alpha_star_human_complete_uci[!is.na(predict.dat$alpha_star_human_complete_uci)] <- scales::rescale(x =predict.dat$alpha_star_human_complete_uci[!is.na(predict.dat$alpha_star_human_complete_uci)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
  
  
  # And do the same for gam.dat
  gam.dat$alpha <- scales::rescale(x =gam.dat$alpha, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
  gam.dat$alpha_lci <- scales::rescale(x =gam.dat$alpha_lci, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
  gam.dat$alpha_uci <- scales::rescale(x =gam.dat$alpha_uci, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
  
  
  # And merge the data
  gam.plot.dat <-  dplyr::select(gam.dat, order, N, alpha, alpha_lci, alpha_uci)
  share.dat.constant <- dplyr::select(predict.dat,order, N, alpha_star_human_constant, alpha_star_human_constant_lci, alpha_star_human_constant_uci)
  share.dat.complete <- dplyr::select(predict.dat,order, N, alpha_star_human_complete, alpha_star_human_complete_lci, alpha_star_human_complete_uci)
  names(share.dat.complete) <- names(share.dat.constant) <- names(gam.plot.dat)
  
  share.dat.constant <- arrange(share.dat.constant, desc(alpha))
  share.dat.complete <- arrange(share.dat.complete, desc(alpha))
  
  
  gam.plot.dat[,3:5] <- -1*gam.plot.dat[,3:5] 
  share.dat.complete$tolerance = "complete"
  share.dat.constant$tolerance = "constant"
  gam.plot.dat$tolerance = "natural"
  share.dat.complete$source <- share.dat.constant$source <-  "predicted from\nnested model"
  gam.plot.dat$source <- "predicted from zoonoses"
  plot.dat <- rbind(gam.plot.dat, share.dat.complete, share.dat.constant)
  plot.dat$source <- factor(plot.dat$source, levels=c("predicted from zoonoses", "predicted from\nnested model"))
  
  # Reorder, ranked by the constant data
  plot.dat$order <- factor(plot.dat$order, levels=unique(arrange(share.dat.constant, desc(alpha))$order))
  
  head(plot.dat)
  
  # And take only the complete data
  plot.dat <- plot.dat[complete.cases(plot.dat),]
  
  
  
  #and now compare predictions
  unique(plot.dat$tolerance)
  unique(plot.dat$source[plot.dat$tolerance=="natural"])
  nat.dat = subset(plot.dat, tolerance=="natural")
  nat.dat <- dplyr::select(nat.dat, order, N, alpha, alpha_lci, alpha_uci)
  names(nat.dat) <- c("order", "lit_N", "lit_alpha", "lit_alpha_lci", "lit_alpha_uci")
  
  dat.compare=subset(plot.dat, tolerance!="natural")
  unique(dat.compare$source)
  head(dat.compare)
  
  dat.compare <- merge(dat.compare, nat.dat, by=c("order"), all.x = T)
  dat.compare$lit_alpha <- dat.compare$lit_alpha*-1 
  
  p1 <- ggplot(data = dat.compare) + theme_bw() +
    geom_point(aes(x=lit_alpha, y=alpha, color=order), size=3) +
    facet_grid(~tolerance)
  
  m1 <- lm(alpha~lit_alpha, dat=subset(dat.compare, tolerance=="complete"))
  summary(m1) #Multiple R-squared:  0.6208. Adjusted R-squared:  0.5577.  p-value: 0.02021
  
  m2 <- lm(alpha~lit_alpha, dat=subset(dat.compare, tolerance=="constant"))
  summary(m2) #Multiple R-squared:  0.5696,	Adjusted R-squared:  0.4979. p-value: 0.03044
  
  #and add to the plot
  dat.compare$predicted_alpha <- NA
  
  dat.compare$predicted_alpha[dat.compare$tolerance=="complete" & !is.na(dat.compare$lit_alpha)] <- predict(m1)
  dat.compare$predicted_alpha[dat.compare$tolerance=="constant" & !is.na(dat.compare$lit_alpha)] <- predict(m2)
  
  droplevels(dat.compare$order)
  
  dat.compare$residuals = dat.compare$alpha-dat.compare$predicted_alpha
  
  return(dat.compare)
  
}

dat.compare <-  assemble.dat(predict.dat = predict.dat, gam.dat = gam.dat)

head(dat.compare)
