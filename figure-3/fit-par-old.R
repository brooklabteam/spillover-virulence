
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

wrap.fit <- function(par.fit, dat, tolerance_type, error_rate, stat.par){
  
  if(par.fit!="all"){
    dat = subset(dat, tolerance==tolerance_type)
    stat.par = subset(stat.par, tolerance==tolerance_type)
    
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

#first, 
dat$Tw =1.5
# Now, fit Tw, holding all else constant
out.fit.1 <- wrap.fit(par.fit = "Tw",  dat=dat,  tolerance_type = "constant", error_rate = .1, stat.par = stat.par)

#next, fit g0 holding all else constant
dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
head(dat)
dat <- dplyr::select(dat, -(alpha_star_human))
dat$g0 = .5
out.fit.2 <- wrap.fit(par.fit = "g0",  dat=dat,  tolerance_type = "constant", error_rate = .1, stat.par = stat.par)

#and fit the phylogenetic distance
dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
head(dat)
dat <- dplyr::select(dat, -(alpha_star_human))
dat$Tv_human[dat$tolerance=="constant"] = 1.5
out.fit.3 <- wrap.fit(par.fit = "Tv_human",  dat=dat,  tolerance_type = "constant", error_rate = .1, stat.par = stat.par)

#and with complete
dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
head(dat)
dat <- dplyr::select(dat, -(alpha_star_human))
dat$Tw =0.5
out.fit.4 <- wrap.fit(par.fit = "Tw",  dat=dat,  tolerance_type = "complete", error_rate = .1, stat.par = stat.par)

dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
head(dat)
dat <- dplyr::select(dat, -(alpha_star_human))
dat$g0 = .5
out.fit.5 <- wrap.fit(par.fit = "g0",  dat=dat,  tolerance_type = "complete", error_rate = .1, stat.par = stat.par)

#and fit the phylogenetic distance
dat <- read.csv(file=paste0(homewd,"figure-3/fit-par-dat.csv"), header = T, stringsAsFactors = F)
head(dat)
dat <- dplyr::select(dat, -(alpha_star_human))
dat$Tv_human[dat$tolerance=="conmplete"] = 0.5
out.fit.6 <- wrap.fit(par.fit = "Tv_human",  dat=dat,  tolerance_type = "complete", error_rate = .1, stat.par = stat.par)


out.all <- rbind(out.fit.1, out.fit.2, out.fit.3, out.fit.4, out.fit.5, out.fit.6)

head(out.all)

alpha.merge <- dplyr::select(stat.par, order, tolerance, alpha_star_human)
names(alpha.merge)[names(alpha.merge)=="alpha_star_human"] <- "alpha_human_fitted_stat_orig"

#alpha.merge$alpha_source <- "stat"

#out.all$alpha_source <- "max_likelihood_fitting"
out.all <- merge(out.all, alpha.merge, by=c("order", "tolerance"))


# First plot a comparioson with the alpha fits
p1 <- ggplot(data=out.all) + theme_bw() +
      geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted, color=order), size=3, shape=16) +
      geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted_stat_orig, color=order), size=3,shape=17) +
      facet_grid(tolerance~fitted_par) + geom_abline()

p1


#compare fits
m1 <- lm(lit_alpha_scale~alpha_human_fitted, data= subset(out.all, tolerance=="constant" & fitted_par=="g0"))
m2 <- lm(lit_alpha_scale~alpha_human_fitted_stat_orig, data= subset(out.all, tolerance=="constant" & fitted_par=="g0")) 

AIC(m1, m2)# original fit is better
summary(m1) #not sig
summary(m2) #still sig and pretty good

library(reshape2)
fit.par.df <- dplyr::select(out.all, order, tolerance, Tw, g0, Tv_human, fitted_par)
head(fit.par.df)
fit.par.df <- melt(fit.par.df, id.vars = c("order", "tolerance", "fitted_par"))

#remove anywhere where fitted_par and variable don't align
fit.par.df <- fit.par.df[fit.par.df$fitted_par==fit.par.df$variable,]
fit.par.df <- dplyr::select(fit.par.df, -(variable))
head(fit.par.df)
names(fit.par.df)[3:4] <- c("fitted_par", "fitted_val")

orig.par.df <- dplyr::select(out.all, order, tolerance, fitted_par, lit_fit)
orig.par.df <- melt(orig.par.df, id.vars = c("order", "tolerance", "fitted_par"))

orig.par.df <- dplyr::select(orig.par.df, -(variable))
names(orig.par.df)[names(orig.par.df)=="value"] <- "orig_stat_val"

fit.par.df <- merge(fit.par.df, orig.par.df, by=c("order", "tolerance", "fitted_par"))
head(fit.par.df)

#then plot a comparison with the fits of each parameter
p2 <- ggplot(data=fit.par.df) + theme_bw() +
  geom_point(aes(x=orig_stat_val, y= fitted_val, color=order), size=3, shape=16) +
  facet_grid(tolerance~fitted_par, scales = "free_y")

p2 #fitting is crap for g0 and for Tv_human constant


#can you fit them all???

out.fit.all.1 <- wrap.fit(par.fit = "all",  dat=dat,  tolerance_type = "constant", error_rate = .1, stat.par = stat.par)
out.fit.all.2 <- wrap.fit(par.fit = "all",  dat=dat,  tolerance_type = "complete", error_rate = .1, stat.par = stat.par)
out.fit.all <- rbind(out.fit.all.1, out.fit.all.2)
head(out.fit.all)



fitted.par <- dplyr::select(out.fit.all, order, tolerance, Tw, g0, Tv_human)
fitted.par <- melt(fitted.par, id.vars = c("order", "tolerance"))
names(fitted.par)[3] <- "fitted_par"


orig.par <- dplyr::select(out.fit.all, order, tolerance, Tw_lit, g0_lit, Tv_human_lit)
orig.par <- melt(orig.par, id.vars = c("order", "tolerance"))
head(orig.par)
names(orig.par)[3:4]<- c("fitted_par", "lit_value")
orig.par$fitted_par <- sapply(strsplit(as.character(orig.par$fitted_par), split="_"), '[',1)

all.fit.par <- merge(fitted.par, orig.par, by=c("order", "tolerance", "fitted_par"))
head(all.fit.par)


#and merge in the alphas
alpha.fits <- dplyr::select(out.fit.all, order,tolerance, lit_alpha_scale, alpha_human_fitted)

#fix here with for-loop
# all.fit.par$lit_alpha_scale <- NA
# 
# 
# all.fit.par <- merge(all.fit.par, alpha.fits, by=c("order, tolerance"))

p3 <- ggplot(data=all.fit.par) + theme_bw() +
  geom_point(aes(x=lit_value, y= value, color=order), size=3, shape=16) +
  facet_grid(tolerance~fitted_par, scales = "free_y")

p3


p4 <- ggplot(data=all.fit.par) + theme_bw() +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted, color=order), size=3, shape=16) +
  geom_point(aes(x=lit_alpha_scale, y= alpha_human_fitted_stat_orig, color=order), size=3,shape=17) +
  facet_grid(tolerance~.)

p4



