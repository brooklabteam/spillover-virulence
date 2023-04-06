
# Goal: hold all parameters constant and fit just the one of interest

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
dat <- read.csv(file=paste0(homewd, "figure-3/compare_predictions.csv"), header = T, stringsAsFactors = F)
dat


# assemble a function that takes in the three parameters of interest 
# (host mortality rate, tolerance of immunopathology, extent of constitutive
# immunity) and produces an estimate of relative virulence

#here are your default parameter values. you can override them as needed for fitting
par.list <- list(b= .2,
                 q = .0002,
                 c=.5,
                 m=1/(21),
                 g1=.9,
                 g0 = .3,
                 zeta=.2,
                 v=1,
                 w=1,
                 Tv=1.005,# .005 for complete tolerance; 1.005 for constant tolerance
                 Tw=1.005,
                 Tv_human=1.005)

# we will practice fitting Tv, Tw, and g0
par.fit = "Tv_human"
#you have to write over the dat vector with the appropriate guess vector for your parameter of interest

# 
# # 
# # dat1 = predict.dat[1,]
# # dat1 <- dplyr::select(dat1, order, mu, mu_lci, mu_uci, N_mu, g0, g0_lci, g0_uci,  N_g0,Tw_constant, Tw_complete,  phylo_dist, Tv_human_constant, Tv_human_complete)
# # dat1 <- dplyr::select(predict.dat, order, mu, mu_lci, mu_uci,  N_mu, g0,  g0_lci, g0_uci, N_g0,Tw_constant, Tw_complete,  phylo_dist, Tv_human_constant, Tv_human_complete)
# # head(dat1)
# 
# Tw_guess = dat1$Tw_constant
# Tw_joke = rep(1.01, length(Tw_guess))
# Tw_high = rep(100, length(Tw_guess))

tolerance_type = "constant"
alt.par.hold=TRUE

# dat <- dplyr::select(dat1, order, mu,  mu_lci, mu_uci, N_mu, g0, g0_lci, g0_uci,  N_g0, phylo_dist, Tv_human_constant, Tv_human_complete)
# 
# dat

#Tw_guess is a vector of guesses in the same order as the orders of the data

make.alpha.prediction <- function(dat, vir.par, alt.par.hold, tolerance_type){
  
  
  if(alt.par.hold==TRUE){
    #g0, Tw, Tv_human are listed in the dat$ data

  
  # Constant tolerance rstar predictions using order-specific g0, Tw, and mu
  # and default values for all other parameters
  if(tolerance_type=="constant"){
    dat$rstar <- (vir.par$c*dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g1*dat$g0*dat$mu*dat$Tw*vir.par$Tv*(vir.par$v*dat$Tw+vir.par$g1*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*dat$Tw) +(vir.par$g1*vir.par$w*(vir.par$m^2)*vir.par$Tv))  

    } else if(tolerance_type=="complete"){
    dat$rstar <- (vir.par$c*dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g1*dat$g0*dat$mu)/(sqrt((vir.par$m^2)*(vir.par$c^2)*dat$mu*vir.par$g1*dat$g0*(vir.par$g1*vir.par$w+vir.par$v-vir.par$g1*dat$Tw-vir.par$Tv)))

    }
  
  
  ## Now, estimate alpha* in spillover human host
  
  # Virulence in the human host is a function of phylogenetic distance,
  # represented in order-specific Tvs (shown here as Tv_human_constant)
  # Tolerance of immunopathology is now human, so we hold constant across 
  # all reservoir host orders.
  
  # First, calculate Vmax in the human host: (this is the same whether constant or complete tolerance)
  dat$Vs_max <- dat$rstar/(vir.par$g1*vir.par$c) - dat$rstar/(2*vir.par$g1*vir.par$c) + 1 - 1/(vir.par$g1) + vir.par$c/(2*dat$rstar*vir.par$g1)

  # Next, take that viral load and calculate spillover virulence, here for the constant tolerance assumption:
  if(tolerance_type=="constant"){
  vir.par$Tw = 1
  dat$alpha_star_human <- ((dat$rstar*vir.par$v)/dat$Tv_human  + (vir.par$g1*vir.par$w*dat$rstar)/vir.par$Tw)*dat$Vs_max

  } else if(tolerance_type=="complete"){
    vir.par$Tw = 0
    dat$alpha_star_human <- (dat$rstar*(vir.par$v-dat$Tv_human) + dat$rstar_complete*vir.par$g1*(vir.par$w-vir.par$Tw))*dat$Vs_max

  }
  } else if(alt.par.hold==FALSE){
    #g0, Tw, Tv_human are listed according to vir.par
    #with only the one that is being optimized changing
    
    
    # Constant tolerance rstar predictions using order-specific g0, Tw, and mu
    # and default values for all other parameters
    if(tolerance_type=="constant"){
      dat$rstar <- (vir.par$c*vir.par$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g1*vir.par$g0*dat$mu*vir.par$Tw*vir.par$Tv*(vir.par$v*vir.par$Tw+vir.par$g1*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*vir.par$Tw) +(vir.par$g1*vir.par$w*(vir.par$m^2)*vir.par$Tv))  
      
    } else if(tolerance_type=="complete"){
      dat$rstar <- (vir.par$c*vir.par$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g1*vir.par$g0*dat$mu)/(sqrt((vir.par$m^2)*(vir.par$c^2)*dat$mu*vir.par$g1*vir.par$g0*(vir.par$g1*vir.par$w+vir.par$v-vir.par$g1*vir.par$Tw-vir.par$Tv)))
      
    }
    
    
    ## Now, estimate alpha* in spillover human host
    
    # Virulence in the human host is a function of phylogenetic distance,
    # represented in order-specific Tvs (shown here as Tv_human_constant)
    # Tolerance of immunopathology is now human, so we hold constant across 
    # all reservoir host orders.
    
    # First, calculate Vmax in the human host: (this is the same whether constant or complete tolerance)
    dat$Vs_max <- dat$rstar/(vir.par$g1*vir.par$c) - dat$rstar/(2*vir.par$g1*vir.par$c) + 1 - 1/(vir.par$g1) + vir.par$c/(2*dat$rstar*vir.par$g1)
    
    # Next, take that viral load and calculate spillover virulence, here for the constant tolerance assumption:
    if(tolerance_type=="constant"){
      vir.par$Tw = 1
      dat$alpha_star_human <- ((dat$rstar*vir.par$v)/dat$Tv_human  + (vir.par$g1*vir.par$w*dat$rstar)/vir.par$Tw)*dat$Vs_max
      
    } else if(tolerance_type=="complete"){
      vir.par$Tw = 0
      dat$alpha_star_human <- (dat$rstar*(vir.par$v-dat$Tv_human) + dat$rstar_complete*vir.par$g1*(vir.par$w-vir.par$Tw))*dat$Vs_max
      
    }
  }
  
  
  #return, and compare with the data
  return(dat)
}
predict.compare.RSS <- function(par, par.list, par.fit, dat, alt.par.hold, tolerance_type){
  
  
  #first, replace the correct element of your parameter set with the ones being optimize
  #first, break your par list down by species
  
  par.df <- cbind.data.frame(lapply(cbind.data.frame(par.list), rep, nrow(dat)))
  par.df$order <- dat$order
  
  dat.merge <- cbind.data.frame(order=dat$order, par.fit=exp(par))

  par.df <- dplyr::select(par.df, -all_of(par.fit))
  par.df <- merge(par.df, dat.merge, by="order")
  
  
  #split dataset by order and extract the corresponding human virulence prediction appropriate 
  dat.split <- dlply(dat,.(order))
  par.split <- dlply(par.df, .(order))
  

  #apply prediction function to data
  out.df <- mapply(FUN= make.alpha.prediction, dat=dat.split, vir.par=par.split,  MoreArgs = list(alt.par.hold=alt.par.hold, tolerance_type=tolerance_type), SIMPLIFY = FALSE)
  
  
  #and get the  lci and uci
  dat.lci <- dat.uci <- dat
  dat.lci$mu <- dat.lci$mu_lci
  dat.uci$mu <- dat.uci$mu_uci
  
  dat.lci$g0 <- dat.lci$g0_lci
  dat.uci$g0 <- dat.uci$g0_uci
  
  dat.split.lci <- dlply(dat.lci,.(order))
  dat.split.uci <- dlply(dat.uci,.(order))
  
  out.df.lci <- mapply(FUN= make.alpha.prediction, dat=dat.split.lci, vir.par=par.split,  MoreArgs = list(alt.par.hold=alt.par.hold, tolerance_type=tolerance_type), SIMPLIFY = FALSE)
  out.df.uci <- mapply(FUN= make.alpha.prediction, dat=dat.split.uci, vir.par=par.split,  MoreArgs = list(alt.par.hold=alt.par.hold, tolerance_type=tolerance_type), SIMPLIFY = FALSE)
  
  
  out.df <- data.table::rbindlist(out.df)
  out.df.lci <- data.table::rbindlist(out.df.lci)
  out.df.uci <- data.table::rbindlist(out.df.uci)
  #head(out.df)
  
  #add uci/lci predictions
  out.df.lci <- dplyr::select(out.df.lci, order, alpha_star_human)
  out.df.uci <- dplyr::select(out.df.uci, order, alpha_star_human)
  names(out.df.lci)[names(out.df.lci)=="alpha_star_human"] <- "alpha_star_human_lci"
  names(out.df.uci)[names(out.df.uci)=="alpha_star_human"] <- "alpha_star_human_uci"
  
  #merge
  out.df <- merge(out.df, out.df.lci, by="order")
  out.df <- merge(out.df, out.df.uci, by="order")
  head(out.df)
  
  
  # And rescale alpha 
  out.df$alpha_star_human[!is.na(out.df$alpha_star_human)] <- scales::rescale(x =out.df$alpha_star_human[!is.na(out.df$alpha_star_human)], from=c(min(out.df$alpha_star_human_lci, na.rm = T), max(out.df$alpha_star_human_uci, na.rm = T)), to =c(0,1)) 
  
  #with(out.df, plot(lit_alpha, alpha_star_human))
  
  out.df$sq_diff <- (out.df$lit_alpha-out.df$alpha_star_human)^2
  
  RSS = sum(out.df$sq_diff)
  
  return(RSS)  
}
wrap.fit <- function(par.fit, par.list, dat, alt.par.hold, tolerance_type){
  
  dat = subset(dat, tolerance==tolerance_type)
  
  par.guess = c(unlist(dat[par.fit]))
  
  #run optim with this parameter list and with the par.fit option
  out.optim <- optim(par = log(par.guess), 
                     fn=predict.compare.RSS, method = "Nelder-Mead", 
                     par.fit=par.fit,
                     par.list=par.list,
                     dat=dat,
                     alt.par.hold=alt.par.hold,
                     tolerance_type=tolerance_type,
                     hessian = T)
  
  dat[par.fit] <- exp(out.optim$par)
  
  dat$convergence=out.optim$convergence
  dat$RSS=out.optim$value
  
return(dat)  
}
  
