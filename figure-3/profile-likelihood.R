
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
library(plotly)

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
profile.llik <- function(par, par.fit, order.fit, dat, tolerance_type, error_rate){
  
  
  
  #first, replace the correct element of your parameter set with the ones being optimized
  dat[dat$order==order.fit,][par.fit] <- exp(par)
  
  
  
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

profile.llik.wrap <- function(par.fit, order.fit, par.range, dat, tolerance_type, error_rate, stat.par){
  
    dat = subset(dat, tolerance==tolerance_type)
    stat.par = subset(stat.par, tolerance==tolerance_type)
    
    llik=list()
    
    for (i in 1:length(par.range)){

    par = par.range[i]
    llik[[i]] <- profile.llik(par=log(par), order.fit = order.fit, dat=dat,par.fit=par.fit, tolerance_type=tolerance_type, error_rate = error_rate)

    }
    
    par.out <- cbind.data.frame(parameter=par.range, neg_llik=unlist(llik), par_fitted = par.fit, order = order.fit, tolerance=tolerance_type)
    
    best.par = par.out$parameter[par.out$neg_llik==min(par.out$neg_llik)]
    
    dat[dat$order==order.fit,][par.fit] <- best.par
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
    
    
    return(par.out)    

}
alpha.range.pred <- function(order.fit, par.fit, par.range, dat, tolerance_type){
  
  alpha.list=list()
  
  print(order.fit)
  
  for (i in 1:length(par.range)){
    print(i)
    
    # pick par from range
    print(par.range[i])
    par.temp <- par.range[i]
    
    #first, replace the correct element of your parameter set with the ones being optimized
    dat.rep = subset(dat, order==order.fit)
    dat.rep[par.fit] <- par.temp
    
    dat.other = subset(dat, order!=order.fit)
    
    dat <- rbind(dat.other, dat.rep)
    
    print(dat$Tw)
    
    
    #make alpha prediction
    tmp <-  make.alpha.prediction(par.dat = dat, tolerance_type = tolerance_type)
    print(tmp)
    
    alpha.tmp <- tmp[tmp$order==order.fit,]$alpha_star_human
    print(alpha.tmp)
    alpha.list[[i]] <- alpha.tmp
    #print(alpha.list[[i]])
    
  }
  
  #now make output of changing alpha given the actual
  alpha.out <- cbind.data.frame(par=par.range, alpha_predict=c(unlist(alpha.list)))
  alpha.out$par_fitted = par.fit
  alpha.out$order = order.fit
  alpha.out$lit_alpha_scale <-  dat[dat$order==order.fit,]$lit_alpha_scale
  
  return(alpha.out)
  
}
alpha.range.pred.Tw <- function(order.fit, par.range, dat, tolerance_type){
  
  alpha.list=list()
  
  #print(order.fit)
  
  for (i in 1:length(par.range)){
   
    # pick par from range
    par.temp <- par.range[i]
    
    #first, replace the correct element of your parameter set with the ones being optimized
    dat$Tw[dat$order == order.fit] <- par.temp
    #print(dat$Tw)
    
    
    #make alpha prediction
    tmp <-  make.alpha.prediction(par.dat = dat, tolerance_type = tolerance_type)
    #print(tmp)
    
    alpha.tmp <- tmp[tmp$order==order.fit,]$alpha_star_human
    #print(alpha.tmp)
    alpha.list[[i]] <- alpha.tmp
    #print(alpha.list[[i]])
    
  }
  
  #now make output of changing alpha given the actual
  alpha.out <- cbind.data.frame(par=par.range, alpha_predict=c(unlist(alpha.list)))
  alpha.out$par_fitted = "Tw"
  alpha.out$order = order.fit
  alpha.out$lit_alpha_scale <-  dat[dat$order==order.fit,]$lit_alpha_scale
  
  return(alpha.out)
  
}
alpha.range.pred.g0 <- function(order.fit, par.range, dat, tolerance_type){
  
  alpha.list=list()
  
  #print(order.fit)
  
  for (i in 1:length(par.range)){
    
    # pick par from range
    par.temp <- par.range[i]
    
    #first, replace the correct element of your parameter set with the ones being optimized
    dat$g0[dat$order == order.fit] <- par.temp
    #print(dat$Tw)
    
    
    #make alpha prediction
    tmp <-  make.alpha.prediction(par.dat = dat, tolerance_type = tolerance_type)
    #print(tmp)
    
    alpha.tmp <- tmp[tmp$order==order.fit,]$alpha_star_human
    #print(alpha.tmp)
    alpha.list[[i]] <- alpha.tmp
    #print(alpha.list[[i]])
    
  }
  
  #now make output of changing alpha given the actual
  alpha.out <- cbind.data.frame(par=par.range, alpha_predict=c(unlist(alpha.list)))
  alpha.out$par_fitted = "g0"
  alpha.out$order = order.fit
  alpha.out$lit_alpha_scale <-  dat[dat$order==order.fit,]$lit_alpha_scale
  
  return(alpha.out)
  
}
alpha.range.pred.Tvhum <- function(order.fit, par.range, dat, tolerance_type){
  
  alpha.list=list()
  
  #print(order.fit)
  
  for (i in 1:length(par.range)){
    
    # pick par from range
    par.temp <- par.range[i]
    
    #first, replace the correct element of your parameter set with the ones being optimized
    dat$Tv_human[dat$order == order.fit] <- par.temp
    #print(dat$Tw)
    
    
    #make alpha prediction
    tmp <-  make.alpha.prediction(par.dat = dat, tolerance_type = tolerance_type)
    #print(tmp)
    
    alpha.tmp <- tmp[tmp$order==order.fit,]$alpha_star_human
    #print(alpha.tmp)
    alpha.list[[i]] <- alpha.tmp
    #print(alpha.list[[i]])
    
  }
  
  #now make output of changing alpha given the actual
  alpha.out <- cbind.data.frame(par=par.range, alpha_predict=c(unlist(alpha.list)))
  alpha.out$par_fitted = "Tv_human"
  alpha.out$order = order.fit
  alpha.out$lit_alpha_scale <-  dat[dat$order==order.fit,]$lit_alpha_scale
  
  return(alpha.out)
  
}
profile.alpha.fit.all <- function(par.range.Tw, par.range.g0, par.range.Tv_human, dat, tolerance_type,  stat.par){
  
  
  dat = subset(dat, tolerance==tolerance_type)
  stat.par = subset(stat.par, tolerance==tolerance_type)
  
  
  
  #wrap over this for all the orders
  order.list = as.list(unique(dat$order))
  
  #first, do Tw
  fit.out.Tw <- lapply(X= order.list, FUN =alpha.range.pred.Tw, dat=dat, par.range=par.range.Tw,  tolerance_type=tolerance_type)
  Tw.df <- data.table::rbindlist(fit.out.Tw)
  
  #then, g0
  fit.out.g0 <- lapply(X= order.list, FUN =alpha.range.pred.g0, dat=dat, par.range=par.range.g0,  tolerance_type=tolerance_type)
  g0.df <- data.table::rbindlist(fit.out.g0)
  
  
  #then Tvhuman
  fit.out.Tv_human <- lapply(X= order.list, FUN =alpha.range.pred.Tvhum, dat=dat, par.range=par.range.Tv_human,  tolerance_type=tolerance_type)
  Tvhuman.df <- data.table::rbindlist(fit.out.Tv_human)
  
  #join together
  fit.df <- rbind(Tw.df, g0.df, Tvhuman.df)
  head(fit.df)
  
  p.out <- ggplot(data=fit.df) + geom_point(aes(x=par, y=alpha_predict, color=order)) + facet_grid(~par_fitted)
  
  plot_ly(y=Tw.df$lit_alpha_scale, z=Tw.df$par, x=Tw.df$alpha_predict, type="scatter3d", mode="markers", color=Tw.df$order)
  plot_ly(x=Tw.df$lit_alpha_scale, y=par, z=Tw.df$alpha_predict, type="scatter3d", mode="markers", color=Tw.df$order)
  
  
  
  scatterplot3d(x=alpha.out$lit_alpha_scale,y=alpha.out$par, z=alpha.out$alpha_predict)
  
}

#first, 
par.range.Tw = seq(1.01, 1.99, length=100)
par.range.g0 = seq(.01, .99, length=100)
par.range.Tv_human = seq(1.01, 1.99, length=100)
tolerance_type = "constant"

# Now, fit Tw, holding all else constant
out.fit.1 <- profile.alpha.fit( par.range = par.range,  dat=dat,  , error_rate = .1, stat.par = stat.par)

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
      facet_grid(tolerance~fitted_par)

p1

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


orig.par <- melt(orig.par, id.vars = c("order", "tolerance"))
orig.par <- dplyr::select(out.fit.all, order, tolerance, Tw_lit, g0_lit, Tv_human_lit)
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



