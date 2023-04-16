
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
dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
head(dat)
stat.par <- dplyr::select(dat, order, tolerance, g0, Tw, Tv_human, alpha_star_human)

dat <- dplyr::select(dat, -(alpha_star_human))

max(dat$Tw[dat$tolerance=="constant"])#1.9
min(dat$Tw[dat$tolerance=="constant"])#1.3

max(dat$Tw[dat$tolerance=="complete"])# 0.9694559
min(dat$Tw[dat$tolerance=="complete"])# 0.3223903


#when you desire to fix pars all the same, write over such that
# dat$g0 = .3
# dat$Tv_human = 1.005 #constant
# dat$Tv_human = 0.005 #complete
# dat$Tw = 1.005 #constant
# dat$Tw = 0.005 #complete

#don't edit mu
make.alpha.prediction <- function(par.dat,  tolerance_type){
  
  
    par.old = par.dat
    # Constant tolerance rstar predictions using order-specific g0, Tw, and mu
    # and default values for all other parameters
    if(tolerance_type=="constant"){
      par.dat$rstar <- (par.dat$c*par.dat$g0/par.dat$m) + (sqrt((par.dat$m^2)*(par.dat$c^2)*par.dat$g1*par.dat$g0*par.dat$mu*par.dat$Tw*par.dat$Tv_bat*(par.dat$v*par.dat$Tw+par.dat$g1*par.dat$w*par.dat$Tv_bat)))/((par.dat$v*(par.dat$m^2)*par.dat$Tw) +(par.dat$g1*par.dat$w*(par.dat$m^2)*par.dat$Tv_bat))  
      
    } else if(tolerance_type=="complete"){
      par.dat$rstar <- (par.dat$c*par.dat$g0/par.dat$m) + ((par.dat$c^2)*par.dat$g1*par.dat$g0*par.dat$mu)/(sqrt((par.dat$m^2)*(par.dat$c^2)*par.dat$mu*par.dat$g1*par.dat$g0*(par.dat$g1*par.dat$w+par.dat$v-par.dat$g1*par.dat$Tw-par.dat$Tv_bat)))
      
    }
    
    
    ## Now, estimate alpha* in spillover human host
    
    # Virulence in the human host is a function of phylogenetic distance,
    # represented in order-specific Tvs (shown here as Tv_human_constant)
    # Tolerance of immunopathology is now human, so we hold constant across 
    # all reservoir host orders.
    
    # First, calculate Vmax in the human host: (this is the same whether constant or complete tolerance)
    par.dat$Vs_max <- par.dat$rstar/(par.dat$g1*par.dat$c) - par.dat$rstar/(2*par.dat$g1*par.dat$c) + 1 - 1/(par.dat$g1) + par.dat$c/(2*par.dat$rstar*par.dat$g1)
    
    # Next, take that viral load and calculate spillover virulence, here for the constant tolerance assumption:
    if(tolerance_type=="constant"){
      par.dat$Tw = 1 #humans have no tolerance
      par.dat$alpha_star_human <- ((par.dat$rstar*par.dat$v)/par.dat$Tv_human  + (par.dat$g1*par.dat$w*par.dat$rstar)/par.dat$Tw)*par.dat$Vs_max
      
    } else if(tolerance_type=="complete"){
      par.dat$Tw = 0 #humans have no tolerance
      par.dat$alpha_star_human <- (par.dat$rstar*(par.dat$v-par.dat$Tv_human) + par.dat$rstar*par.dat$g1*(par.dat$w-par.dat$Tw))*par.dat$Vs_max
      
    }
  
  par.dat$Tw = par.old$Tw
  #return, and compare with the data
  return(par.dat)
} #just finished editing this... 
#fix the rest of functions
predict.compare.llik <- function(par, par.fit, dat, tolerance_type, error_rate){
  
  if(par.fit=="all"){
    
    #first, replace the correct element of your parameter set with the ones being optimized
    n.par = length(par)/3
    dat$g0 <- exp(par[1:n.par])
    dat$Tw <- exp(par[(n.par+1):(2*n.par)])
    dat$Tv_human <- exp(par[(2*n.par+1):(n.par*3)])
  }else{
    
  
  
  #first, replace the correct element of your parameter set with the ones being optimized
  dat[par.fit] <- exp(par)
  
  }
  
  #split dataset by order and extract the corresponding human virulence prediction appropriate 
  dat.split <- dlply(dat,.(order))
  
  

  #apply prediction function to data
  out.df <- lapply(dat.split,  make.alpha.prediction, tolerance_type=tolerance_type)
  
  
  out.df <- data.table::rbindlist(out.df)
 
  #with(out.df, plot(lit_alpha_scale, alpha_star_human))
  
  #and get likelihoods of model fit to dat - this should be a perfect 1:1
  likelihoods <- dnorm(out.df$alpha_star_human, mean = out.df$lit_alpha_scale, sd = error_rate)
  
  #correct if 0
  
  likelihoods[likelihoods==0] <- .00000000000001
  
  #how likely is the data, given your model? try a linear regression
  #m1 <- lm(alpha_star_human~lit_alpha, dat=out.df)
  #summary(m1)
  
  #and minimize neg llik *2
  
  neg_llik = -2*log(sum(likelihoods))
  
  # #add constraints
  # if((par.fit=="Tw" & tolerance_type=="constant" & length(par[par>2])>0) | (par.fit=="Tw" & tolerance_type=="constant" & length(par[par<1])>0)){
  #      
  # neg_llik=1000000
  # 
  # }
  # 
  #out.df$sq_diff <- (out.df$lit_alpha-out.df$alpha_star_human)^2
  
  #RSS = sum(out.df$sq_diff)
  
  return(neg_llik)  
}
predict.compare.llik.single <- function(par, par.fit, dat, tolerance_type, error_rate){
  
  
    
    
    
    #first, replace the correct element of your parameter set with the ones being optimized
    dat[par.fit] <- par
    
  
  
  #split dataset by order and extract the corresponding human virulence prediction appropriate 
  dat.split <- dlply(dat,.(order))
  
  
  
  #apply prediction function to data
  out.df <- lapply(dat.split,  make.alpha.prediction, tolerance_type=tolerance_type)
  
  
  out.df <- data.table::rbindlist(out.df)
  
  #with(out.df, plot(lit_alpha_scale, alpha_star_human))
  
  #and get likelihoods of model fit to dat - this should be a perfect 1:1
  likelihoods <- dnorm(out.df$alpha_star_human, mean = out.df$lit_alpha_scale, sd = error_rate)
  
  #correct if 0
  
  likelihoods[likelihoods==0] <- .00000000000001
  
  #how likely is the data, given your model? try a linear regression
  #m1 <- lm(alpha_star_human~lit_alpha, dat=out.df)
  #summary(m1)
  
  #and minimize neg llik *2
  
  neg_llik = -2*log(sum(likelihoods))
  
  # #add constraints
  # if((par.fit=="Tw" & tolerance_type=="constant" & length(par[par>2])>0) | (par.fit=="Tw" & tolerance_type=="constant" & length(par[par<1])>0)){
  #      
  # neg_llik=1000000
  # 
  # }
  # 
  #out.df$sq_diff <- (out.df$lit_alpha-out.df$alpha_star_human)^2
  
  #RSS = sum(out.df$sq_diff)
  
  return(neg_llik)  
}
optimize.single.par <- function(par, dat, lower.par, upper.par, par.fit, tolerance_type, error_rate){
  
  par.range = seq(lower.par, upper.par, length=100)
  print(dat$order)
  print(par.fit)
  
  
  # #run optim with this parameter list and with the par.fit option
  
  
  llik = list()
  
  for (i in 1:length(par.range)){
    
  par.tmp = par.range[i]  
  llik[[i]] = predict.compare.llik.single(par =par.tmp, dat=dat, par.fit=par.fit, tolerance_type=tolerance_type, error_rate = error_rate)  
  
  }
  
  #and combine
  par.fit.df <- cbind.data.frame(par_df = par.range, llik=unlist(llik))
  #with(par.fit.df, plot(par_df, llik, type="l"))
  
  par.best <- par.fit.df$par_df[par.fit.df$llik==min(par.fit.df$llik)]
  
  if(length(par.best)>1){
    print("likelihood landscape flat: used optim")
     out.optim <- optim(par = par, 
                        fn=predict.compare.llik.single, 
                        method = "Brent", 
                        lower=lower.par,
                        upper = upper.par,
                        dat=dat,
                        par.fit=par.fit,
                        tolerance_type=tolerance_type,
                        error_rate = error_rate,
                        hessian = T)
     
     print(out.optim$par)
     
     dat[par.fit] <- out.optim$par
     
     dat.new =  make.alpha.prediction(par.dat = dat, tolerance_type = tolerance_type)
     names(dat.new)[names(dat.new)=="alpha_star_human"] <- "alpha_human_fitted"
     
     dat.new$llik = out.optim$value
    
  }else{
    print("successfully profiled likelihood surface")
    print(par.best)
    
    
    dat[par.fit] <- par.best
    
    #and what virulence results from this?
    dat.new =  make.alpha.prediction(par.dat = dat, tolerance_type = tolerance_type)
    names(dat.new)[names(dat.new)=="alpha_star_human"] <- "alpha_human_fitted"
    #dat.new$convergence = out.optim$convergence
    dat.new$llik = par.fit.df$llik[par.fit.df$par_df==par.best]
    
  }
  
  print(dat.new$alpha_human_fitted)
  
  return(dat.new)
  
}

wrap.fit <- function(par.fit, dat, single.par, lower.par, upper.par, tolerance_type, error_rate, stat.par){
  
  if(par.fit!="all"){
    dat = subset(dat, tolerance==tolerance_type)
    stat.par = subset(stat.par, tolerance==tolerance_type)
    
    if(single.par==FALSE){
      
      par.guess = c(unlist(dat[par.fit]))
      
      #run optim with this parameter list and with the par.fit option
      out.optim <- optim(par = log(par.guess), 
                         fn=predict.compare.llik, 
                         method = "Nelder-Mead", 
                         dat=dat,
                         par.fit=par.fit,
                         tolerance_type=tolerance_type,
                         error_rate = error_rate,
                         hessian = T)
      
      
      #plot(unlist(dat[par.fit]), exp(out.optim$par))
      dat[par.fit] <- exp(out.optim$par)
      
      
      dat$convergence=out.optim$convergence
      dat$llik=out.optim$value
      
      #add new alpha produced from this estimate
      #split dataset by order and extract the corresponding human virulence prediction appropriate 
      dat.split <- dlply(dat,.(order))
      
      
      #apply prediction function to data
      new.out  <- lapply(dat.split,  make.alpha.prediction, tolerance_type=tolerance_type)
      new.df <- data.table::rbindlist(new.out)
      new.df <- dplyr::select(new.df, order, alpha_star_human)
      names(new.df) <- c("order", "alpha_human_fitted")
      
      dat <- merge(dat, new.df, by="order")
      #with(new.df, plot(lit_alpha_scale, alpha_star_human)) 
      
      #return data with new parameters intact
      
      dat$fitted_par=par.fit
      
      #and merge in the old stat.par
      stat.merge <- dplyr::select(stat.par, order, par.fit)
      names(stat.merge)[length(names(stat.merge))] <- "lit_fit" #paste0(par.fit,"_lit")
      dat <- merge(dat, stat.merge)
      
      return(dat)  
      
    }else if (single.par==TRUE){
      
    
      
      par.guess.list = as.list(c(unlist(dat[par.fit])))
      
      dat.list <- dlply(dat, .(order))
      
      #apply fitting script over all
      fitted.list <- mapply(optimize.single.par, par=par.guess.list, dat=dat.list, 
                          MoreArgs= list(lower.par=lower.par, upper.par=upper.par, 
                                         error_rate=error_rate,
                          par.fit=par.fit, tolerance_type= tolerance_type), SIMPLIFY = FALSE)
      
      
      fit.df <- data.table::rbindlist(fitted.list)
      #ggplot(fit.df) + geom_point(aes(x=lit_alpha_scale,  y=alpha_human_fitted, color=order))
      
      #and subset to just those of interest
      
      #return()
      
      #and merge in the old stat.par
      stat.merge <- dplyr::select(stat.par, order, par.fit, alpha_star_human)
      names(stat.merge)[(length(names(stat.merge))-1)] <- "lit_fit" #paste0(par.fit,"_lit")
      names(stat.merge)[length(names(stat.merge))] <- "alpha_human_fitted_orig_stat"
      fit.df <- merge(fit.df, stat.merge)
      
      head(fit.df)
      fit.df <- dplyr::select(fit.df, order, tolerance, par.fit, lit_fit, lit_alpha_scale, alpha_human_fitted, alpha_human_fitted_orig_stat)
      names(fit.df) <- c("order", "tolerance", "fitted_par", "lit_fit", "lit_alpha_scale", "alpha_human_fitted", "alpha_human_fitted_orig_stat")
      
      fit.df$fit_par_ID <- par.fit
      
      return(fit.df)  
      
      
    }
     
    
  }else if(par.fit=="all"){
    
    dat = subset(dat, tolerance==tolerance_type)
    stat.par = subset(stat.par, tolerance==tolerance_type)
    
    par.guess = c(unlist(dat$g0), unlist(dat$Tw), unlist(dat$Tv_human))
    
    #run optim with this parameter list and with the par.fit option
    out.optim <- optim(par = log(par.guess), 
                       fn=predict.compare.llik, 
                       method = "Nelder-Mead", 
                       dat=dat,
                       par.fit=par.fit,
                       tolerance_type=tolerance_type,
                       error_rate = error_rate,
                       hessian = T)
    
    n.par = length(par.guess)/3
    dat$g0 <- exp(out.optim$par[1:n.par])
    dat$Tw <- exp(out.optim$par[(n.par+1):(2*n.par)])
    dat$Tv_human <- exp(out.optim$par[(2*n.par+1):(n.par*3)])
    
    
    
    dat$convergence=out.optim$convergence
    dat$llik=out.optim$value
    
    #add new alpha produced from this estimate
    #split dataset by order and extract the corresponding human virulence prediction appropriate 
    dat.split <- dlply(dat,.(order))
    
    
    #apply prediction function to data
    new.out  <- lapply(dat.split,  make.alpha.prediction, tolerance_type=tolerance_type)
    new.df <- data.table::rbindlist(new.out)
    new.df <- dplyr::select(new.df, order, alpha_star_human)
    names(new.df) <- c("order", "alpha_human_fitted")
    
    dat <- merge(dat, new.df, by="order")
    #with(dat, plot(lit_alpha_scale, alpha_human_fitted)) 
    
    #return data with new parameters intact
    
    dat$fitted_par=par.fit
    
    #and merge in the old stat.par
    stat.merge <- dplyr::select(stat.par, order, g0, Tw, Tv_human)
    names(stat.merge)[2:length(names(stat.merge))] <- paste0(names(stat.merge)[2:length(names(stat.merge))],"_lit")
    dat <- merge(dat, stat.merge)
    
    return(dat)    
    
    
  }
  
}

#first, fit Tw individually for each order, holding everything else constant

out.fit.1 <- wrap.fit(par.fit = "Tw",  dat=dat, single.par=TRUE,  lower.par = 1.00001, 
                      upper.par = 1.99999, tolerance_type = "constant", 
                      error_rate = .5, stat.par = stat.par)

#and g0
out.fit.2 <- wrap.fit(par.fit = "g0",  dat=dat, single.par=TRUE,  lower.par = .00001, 
                      upper.par = 0.99999, tolerance_type = "constant", 
                      error_rate = .5, stat.par = stat.par)
#and Tv_human
out.fit.3 <- wrap.fit(par.fit = "Tv_human",  dat=dat, single.par=TRUE,  lower.par = 1.00001, 
                      upper.par = 1.99999, tolerance_type = "constant", 
                      error_rate = .5, stat.par = stat.par)

#and repeat at complete
out.fit.4 <- wrap.fit(par.fit = "Tw",  dat=dat, single.par=TRUE,  lower.par = 0.00001, 
                      upper.par = 0.99999, tolerance_type = "complete", 
                      error_rate = .5, stat.par = stat.par)

#and g0
out.fit.5 <- wrap.fit(par.fit = "g0",  dat=dat, single.par=TRUE,  lower.par = .00001, 
                      upper.par = 0.99999, tolerance_type = "complete", 
                      error_rate = .5, stat.par = stat.par)
#and Tv_human
out.fit.6 <- wrap.fit(par.fit = "Tv_human",  dat=dat, single.par=TRUE,  lower.par = 0.00001, 
                      upper.par = 0.99999, tolerance_type = "complete", 
                      error_rate = .5, stat.par = stat.par)


#and join
out.fit <- rbind(out.fit.1, out.fit.2, out.fit.3, out.fit.4, out.fit.5, out.fit.6)
head(out.fit)

out.fit$tolerance <- factor(out.fit$tolerance, levels=c("constant", "complete"))

# First plot a comparison with the alpha fits
# this is with the optimized value for each of the three parameters,
# holding the others as they already vary in the literature

p1 <- ggplot(data=out.fit) + theme_bw() + 
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted, color=order), size=3, shape=16) +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted_orig_stat, color=order), size=3,shape=17) +
  facet_grid(tolerance~fit_par_ID) + geom_abline()

p1 

#original fits predicted too high of virulence for Eulipotyphla. 
#we can rescue this by maniupating g0, or to a lesser extent, Tv_human



#and versus the old fits

p2 <- ggplot(data=out.fit) + theme_bw() +
  geom_point(aes( x= lit_fit,y=fitted_par, color=order), size=3, shape=16) +
  facet_wrap(tolerance~fit_par_ID, scales = "free") + geom_abline()

p2 
#g0 has to be much lower for Eulipotyphla to rescue relationshio - or Tv_human has to be much higher
#Tw does not have much pull on the parameter interactions


# # then, repeat when everything is held still for all others

dat$Tw[dat$tolerance=="constant"] <- 1.5
dat$Tw[dat$tolerance=="complete"] <- 0.5

dat$Tv_human[dat$tolerance=="constant"] <- 1.5
dat$Tv_human[dat$tolerance=="complete"] <- 0.5

dat$g0 <- .05


out.fit.1 <- wrap.fit(par.fit = "Tw",  dat=dat, single.par=TRUE,  lower.par = 1.00001, 
                      upper.par = 1000000, tolerance_type = "constant", 
                      error_rate = .5, stat.par = stat.par)

#and g0
out.fit.2 <- wrap.fit(par.fit = "g0",  dat=dat, single.par=TRUE,  lower.par = .00001, 
                      upper.par = 0.99999, tolerance_type = "constant", 
                      error_rate = .5, stat.par = stat.par)
#and Tv_human
out.fit.3 <- wrap.fit(par.fit = "Tv_human",  dat=dat, single.par=TRUE,  lower.par = 1.00001, 
                      upper.par = 1.99999, tolerance_type = "constant", 
                      error_rate = .5, stat.par = stat.par)

#and repeat at complete
out.fit.4 <- wrap.fit(par.fit = "Tw",  dat=dat, single.par=TRUE,  lower.par = 0.00001, 
                      upper.par = 0.99999, tolerance_type = "complete", 
                      error_rate = .5, stat.par = stat.par)

#and g0
out.fit.5 <- wrap.fit(par.fit = "g0",  dat=dat, single.par=TRUE,  lower.par = .00001, 
                      upper.par = 0.99999, tolerance_type = "complete", 
                      error_rate = .5, stat.par = stat.par)
#and Tv_human
out.fit.6 <- wrap.fit(par.fit = "Tv_human",  dat=dat, single.par=TRUE,  lower.par = 0.00001, 
                      upper.par = 0.99999, tolerance_type = "complete", 
                      error_rate = .5, stat.par = stat.par)


#and join
out.fit <- rbind(out.fit.1, out.fit.2, out.fit.3, out.fit.4, out.fit.5, out.fit.6)
head(out.fit)



# First plot a comparison with the alpha fits
p3 <- ggplot(data=out.fit) + theme_bw() +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted, color=order), size=3, shape=16) +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted_orig_stat, color=order), size=3,shape=17) +
  facet_grid(tolerance~fit_par_ID) + geom_abline()

p3

#We can't really get te virulence we need for bats without allowing for variation in g0 AND Tv_human
#though g0 is stronger. When g0 is fixed at low, we can correct for Eulipotyphla over estimates

#and versus the old fits

p4 <- ggplot(data=out.fit) + theme_bw() +
  geom_point(aes(y=fitted_par, x= lit_fit, color=order), size=3, shape=16) +
  facet_wrap(tolerance~fit_par_ID, scales = "free") + geom_abline()

p4 

# Prior work overestimated g0 for Eulipotyphla and overestimated Tv_human for bats
# Fitting suggests that natural populations, excepting bats and carnivores, 
# should be totally intolerant of the viruses they evolve
# 

