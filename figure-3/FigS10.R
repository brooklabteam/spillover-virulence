
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

#load parameters from GlM analysis
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
predict.compare.llik.single.alpha <- function(par, par.fit, dat, tolerance_type, error_rate){
  
  
  
  
  
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
  
  #and return both the parameter and the alpha prediction
  out.df$neg_llik <- neg_llik
  
  return(out.df)  
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
profile.single.par <- function(par, dat, lower.par, upper.par, par.fit, tolerance_type, error_rate){
  
  par.range = seq(lower.par, upper.par, length=100)
  #print(dat$order)
  #print(par.fit)
  
  
  # #run optim with this parameter list and with the par.fit option
  
  
  out.df = list()
  
  for (i in 1:length(par.range)){
    
    par.tmp = par.range[i]  
    out.df[[i]] = predict.compare.llik.single.alpha(par =par.tmp, dat=dat, par.fit=par.fit, tolerance_type=tolerance_type, error_rate = error_rate)  
    
  }
  
  out.df <- data.table::rbindlist(out.df)
  #with(out.df, plot(g0, alpha_star_human))
  
  out.df$par_fit <- par.fit
  return(out.df)
  
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

#first, fit each of the 4 modulated parameters individually for each order,
#holding everything else constant

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
out.fit.A <- rbind(out.fit.1, out.fit.2, out.fit.3, out.fit.4, out.fit.5, out.fit.6)
head(out.fit.A)

out.fit.A$tolerance <- factor(out.fit.A$tolerance, levels=c("constant", "complete"))

# First plot a comparison with the alpha fits
# this is with the optimized value for each of the three parameters,
# holding the others as they already vary in the literature

p1 <- ggplot(data=out.fit.A) + theme_bw() + 
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted, color=order), size=3, shape=16) +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted_orig_stat, color=order), size=3,shape=17) +
  facet_grid(tolerance~fit_par_ID) + geom_abline()

p1 

#original fits predicted too high of virulence for Eulipotyphla. 
#we can rescue this by lowering the g0 estimate or increasing the Tv_human
# value. Modulating Tw has minimal impact




#and versus the old fits

p2 <- ggplot(data=out.fit.A) + theme_bw() +
  geom_point(aes( x= lit_fit,y=fitted_par, color=order), size=3, shape=16) +
  facet_wrap(tolerance~fit_par_ID, scales = "free") + geom_abline()

p2 

#g0 has to be much lower for Eulipotyphla to rescue relationshio -
#or Tv_human has to be much higher. 
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
out.fit.B <- rbind(out.fit.1, out.fit.2, out.fit.3, out.fit.4, out.fit.5, out.fit.6)
head(out.fit.B)



# First plot a comparison with the alpha fits
p3 <- ggplot(data=out.fit.B) + theme_bw() +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted, color=order), size=3, shape=16) +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted_orig_stat, color=order), size=3,shape=17) +
  facet_grid(tolerance~fit_par_ID) + geom_abline()

p3

# We can't really get te virulence we need for bats without allowing for 
# variation in g0 AND Tv_human though g0 does drive the realtionship even more
# When g0 is fixed at low, we can correct for Eulipotyphla over estimates

#and versus the old fits

p4 <- ggplot(data=out.fit.B) + theme_bw() +
  geom_point(aes(y=fitted_par, x= lit_fit, color=order), size=3, shape=16) +
  facet_wrap(tolerance~fit_par_ID, scales = "free") + geom_abline()

p4 


#now join and make plot


#melt to modulate alpha - we have an original fitted alpha and 
#a new estimate for each
alpha.melt.A <- melt(out.fit.A, id.vars = c("order", "tolerance", "fitted_par", "lit_fit", "fit_par_ID", "lit_alpha_scale"))

names(alpha.melt.A)[names(alpha.melt.A)=="variable"] <- "alpha_type"
names(alpha.melt.A)[names(alpha.melt.A)=="value"] <- "alpha"
alpha.melt.orig = subset(alpha.melt.A, alpha_type=="alpha_human_fitted_orig_stat")
#alpha.melt.orig <- dplyr::select(alpha.melt.orig, order, tolerance, lit_alpha_scale, alpha_type, alpha)
alpha.melt.orig <- alpha.melt.orig[!duplicated(alpha.melt.orig),]
alpha.melt.orig$alpha_type <- "glm fit for\nall parameters"

alpha.melt.A = subset(alpha.melt.A, alpha_type!="alpha_human_fitted_orig_stat")
alpha.melt.A$alpha_type <- "one parameter profiled;\nother parameters\nfrom glm fitting"

alpha.melt.B <- melt(out.fit.B, id.vars = c("order", "tolerance", "fitted_par", "lit_fit", "fit_par_ID", "lit_alpha_scale"))

names(alpha.melt.B)[names(alpha.melt.B)=="variable"] <- "alpha_type"
names(alpha.melt.B)[names(alpha.melt.B)=="value"] <- "alpha"

alpha.melt.B = subset(alpha.melt.B, alpha_type!="alpha_human_fitted_orig_stat")
alpha.melt.B$alpha_type <- "one parameter profiled;\nother parameters\nheld constant"


alpha.fit <- rbind(alpha.melt.A, alpha.melt.B, alpha.melt.orig)

head(alpha.fit)


#and plot 
shapez = c('glm fit for\nall parameters' = 21, 'one parameter profiled;\nother parameters\nfrom glm fitting' = 24, 'one parameter profiled;\nother parameters\nheld constant'=22)
colz2 = c('glm fit for\nall parameters' = "red4", 'one parameter profiled;\nother parameters\nfrom glm fitting' = "black", 'one parameter profiled;\nother parameters\nheld constant'="black")

load(paste0(homewd, "/figure-3/color-bar.Rdata"))

alpha.fit$tolerance <- as.character(alpha.fit$tolerance)
alpha.fit$tolerance[alpha.fit$tolerance=="constant"] <- "'constant tolerance'"
alpha.fit$tolerance[alpha.fit$tolerance=="complete"] <- "'complete tolerance'"
alpha.fit$tolerance <- factor(alpha.fit$tolerance, levels=c("'constant tolerance'", "'complete tolerance'"))

# and add labels for strip labels
alpha.fit$label= NA

alpha.fit$label[alpha.fit$fit_par_ID=="Tw"] <- "atop(T[w]~', reservoir tolerance', 'of immunopathology')"
alpha.fit$label[alpha.fit$fit_par_ID=="Tv_human"] <- "atop(T[v]~', spillover host tolerance', 'of direct virus pathology')"
alpha.fit$label[alpha.fit$fit_par_ID=="g0"] <- "atop(g[0]~', constitutive ', 'immunity in reservoir')"




alpha.fit$label = factor(alpha.fit$label, levels = c("atop(T[w]~', reservoir tolerance', 'of immunopathology')",
                                                     "atop(T[v]~', spillover host tolerance', 'of direct virus pathology')",
                                                     "atop(g[0]~', constitutive ', 'immunity in reservoir')"))


label_parseall <- function(variable, value) {
  plyr::llply(value, function(x) parse(text = paste("", 
                                                    x, sep = "")))
}

pFigS10 <- ggplot(data=alpha.fit) + theme_bw() +
  geom_point(aes(x=lit_alpha_scale, y= alpha, fill=order, 
                 shape=alpha_type, color=alpha_type), size=3, stroke=1) +
  scale_shape_manual(values=shapez) +
  scale_color_manual(values=colz2) + 
  scale_fill_manual(values=colz) + 
  theme(strip.background = element_rect(fill="white"), 
        axis.title = element_text(size=14), 
        strip.text= element_text(size=14), 
        panel.grid = element_blank(),
        axis.text = element_text(size=12)) +
  facet_grid(tolerance~label,labeller = label_parseall) + geom_abline() +
  ylab(bquote("predicted spillover virulence,"~alpha[S])) + 
  xlab(bquote("observed spillover virulence,"~alpha[S])) +
  guides(fill = guide_legend(override.aes=list(shape=21)), 
         color=guide_legend(title = bquote(alpha[S]~"estimation method"), keywidth=0.1,keyheight=0.5,default.unit="inch"), 
         shape=guide_legend(title = bquote(alpha[S]~"estimation method"), keywidth=0.1,keyheight=0.5,default.unit="inch")) #+
  

pFigS10

ggsave(file = paste0(homewd,"/supp-figs/FigS10.png"),
       plot = pFigS10,
       units="mm",  
       width=120, 
       height=70, 
       scale=3, 
       dpi=300)
# 
# 
# merge(alpha.fit, alpha.melt.orig, by=c("order", "tolerance", "lit_alpha_scale", "alpha_type"))
# out.fit
# 
# #and melt to modulate alpha
# 
# alpha.melt <- melt(out.fit, id.vars = c("order", "tolerance", "fitted_par", "lit_fit", "fit_par_ID", "lit_alpha_scale", "method"))
# 
# names(alpha.melt)[names(alpha.melt)=="variable"] <- "alpha_form"
# names(alpha.melt)[names(alpha.melt)=="value"] <- "alpha"
# alpha.melt$alpha_form <- as.character(alpha.melt$alpha_form)
# #alpha.melt$alpha_form[alpha.melt$alpha_form=="lit_alpha_scale"] <- "alpha-literature"
# alpha.melt$alpha_form[alpha.melt$alpha_form=="alpha_human_fitted_orig_stat"] <- "glm-fit"
# alpha.melt$alpha_form[alpha.melt$alpha_form=="alpha_human_fitted"] <- "profiled-parameter"
# 
# alpha.melt.orig <- subset(alpha.melt, alpha_form=="glm-fit")
# alpha.melt.orig$method = ""
# alpha.melt.orig$fitted_par <- NA
# alpha.melt.orig$lit_fit <- round(alpha.melt.orig$lit_fit, 4)
# alpha.melt.orig$lit_alpha_scale <- round(alpha.melt.orig$lit_alpha_scale, 4)
# alpha.melt.orig$alpha <- round(alpha.melt.orig$alpha, 4)
# 
# alpha.melt.orig <- alpha.melt.orig[!duplicated(alpha.melt.orig),]
# 
# alpha.melt
# alpha.melt$alpha_form <- factor(alpha.melt$alpha_form, levels=c("glm-fit", "profiled-parameter" ))
# 
# #make a constant and a complete plot
# 
# constant <- ggplot(data=subset(alpha.melt, tolerance=="constant")) + theme_bw() +
#   geom_point(aes(x=lit_alpha_scale, y= alpha, color=order, shape=alpha_form), size=3) +
#   facet_grid(method~fit_par_ID) + geom_abline() +ylab("predicted alpha") + xlab("observed alpha")
# 
# constant
# # Prior work overestimated g0 for Eulipotyphla and overestimated Tv_human for bats
# # Fitting suggests that natural populations, excepting bats and carnivores, 
# # should be totally intolerant of the viruses they evolve
# # 
# 
# # And try one other form of the plot...
# 
# 
# profile.plot <- function(par.fit, dat, lower.par, upper.par, tolerance_type, error_rate, stat.par){
#  
#     dat = subset(dat, tolerance==tolerance_type)
#     stat.par = subset(stat.par, tolerance==tolerance_type)
#   
#       
#       
#       par.guess.list = as.list(c(unlist(dat[par.fit])))
#       
#       dat.list <- dlply(dat, .(order))
#       
#       #apply fitting script over all
#       fitted.list <- mapply(profile.single.par, par=par.guess.list, dat=dat.list, 
#                             MoreArgs= list(lower.par=lower.par, upper.par=upper.par, 
#                                            error_rate=error_rate,
#                                            par.fit=par.fit, tolerance_type= tolerance_type), SIMPLIFY = FALSE)
#       
#       
#       fit.df <- data.table::rbindlist(fitted.list)
#       #ggplot(fit.df) + geom_point(aes(x=lit_alpha_scale,  y=alpha_human_fitted, color=order))
#       
#       #and subset to just those of interest
#       
#       #return()
#       
#       #and merge in the old stat.par
#       stat.merge <- dplyr::select(stat.par, order, par.fit, alpha_star_human)
#       names(stat.merge)[(length(names(stat.merge))-1)] <- "lit_fit" #paste0(par.fit,"_lit")
#       names(stat.merge)[length(names(stat.merge))] <- "alpha_human_fitted_orig_stat"
#       fit.df <- merge(fit.df, stat.merge)
#       
#       head(fit.df)
#       names(fit.df)[names(fit.df)=="alpha_star_human"] <- "alpha_human_fitted"
#       names(fit.df)[names(fit.df)==par.fit] <- "par.fit"
#       fit.df <- dplyr::select(fit.df, order, tolerance, par.fit, lit_fit, lit_alpha_scale, alpha_human_fitted, alpha_human_fitted_orig_stat)
#       names(fit.df) <- c("order", "tolerance", "fitted_par", "lit_fit", "lit_alpha_scale", "alpha_human_fitted", "alpha_human_fitted_orig_stat")
#       
#       fit.df$fit_par_ID <- par.fit
#       
#       
#       
#       
#       return(fit.df)  
#       
#       
#     
#     
#   
#   
# }
# 
# dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
# head(dat)
# stat.par <- dplyr::select(dat, order, tolerance, g0, Tw, Tv_human, alpha_star_human)
# 
# dat <- dplyr::select(dat, -(alpha_star_human))
# 
# 
# #Tw 
# out.fit.7 <- profile.plot(par.fit = "Tw",  dat=dat,   1.00001, 
#                           upper.par = 1.99999, tolerance_type = "constant", 
#                           error_rate = .5, stat.par = stat.par)
# 
# #and g0
# out.fit.8 <- profile.plot(par.fit = "g0",  dat=dat,   lower.par = .00001, 
#                       upper.par = 0.99999, tolerance_type = "constant", 
#                       error_rate = .5, stat.par = stat.par)
# 
# 
# #Tv_human 
# out.fit.9 <- profile.plot(par.fit = "Tv_human",  dat=dat, lower.par =  1.00001, 
#                           upper.par = 1.99999, tolerance_type = "constant", 
#                           error_rate = .5, stat.par = stat.par)
# #and complete
# 
# #Tw 
# out.fit.10 <- profile.plot(par.fit = "Tw",  dat=dat,  lower.par = 0.00001, 
#                           upper.par = 0.99999, tolerance_type = "complete", 
#                           error_rate = .5, stat.par = stat.par)
# 
# #and g0
# out.fit.11 <- profile.plot(par.fit = "g0",  dat=dat,   lower.par = .00001, 
#                           upper.par = 0.99999, tolerance_type = "complete", 
#                           error_rate = .5, stat.par = stat.par)
# 
# 
# #Tv_human 
# out.fit.12 <- profile.plot(par.fit = "Tv_human",  dat=dat,   0.00001, 
#                           upper.par = 0.99999, tolerance_type = "complete", 
#                           error_rate = .5, stat.par = stat.par)
# 
# 
# 
# out.fit <- rbind(out.fit.7, out.fit.8, out.fit.9, out.fit.10, out.fit.11, out.fit.12)
# head(out.fit)
# 
# 
# out.fit.sub = subset(out.fit, alpha_human_fitted<2)
# out.fit.sub$tolerance <- factor(out.fit.sub$tolerance, levels=c("constant", "complete"))
# p5 <- ggplot(out.fit.sub) + geom_point(aes(x=fitted_par, y=alpha_human_fitted, color=order)) +
#   geom_point(aes(x=lit_fit, y=  alpha_human_fitted_orig_stat, fill=order), size=3, shape=24) +
#   geom_point(aes(x=lit_fit, y=  lit_alpha_scale, fill=order), size=3, shape=22) +
#   facet_wrap(tolerance~fit_par_ID, scales = "free") 
# p5
# 
# 
# #and 