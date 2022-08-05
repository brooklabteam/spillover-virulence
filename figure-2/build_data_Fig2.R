rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(deSolve)
library(reshape2)

# Put this to your own wd
homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-2"
setwd(paste0(homewd, subwd))

# Run the model to get the output - tolerance held steady and variation in r* and alpha*
get.rstar.general <- function(par.dat, tol.type){
  
  
  if(tol.type=="complete-tolerance"){
    
    
    rstar = (par.dat$c*par.dat$g0/par.dat$m) + ((par.dat$c^2)*par.dat$g*par.dat$g0*par.dat$mu)/(sqrt((par.dat$m^2)*(par.dat$c^2)*par.dat$mu*par.dat$g*par.dat$g0*(par.dat$g*par.dat$w+par.dat$v-par.dat$g*par.dat$Tw-par.dat$Tv)))
    
     
    return(rstar)
  }else if(tol.type=="constant-tolerance"){
    
    
    #rstar = ((par.dat$c*par.dat$g0*par.dat$m*par.dat$v*par.dat$Tw)+sqrt(par.dat$c*par.dat$v*par.dat$g*par.dat$Tv*par.dat$Tw*(par.dat$m^2)*(par.dat$m*par.dat$mu*par.dat$Tw+par.dat$g0*par.dat$w)))/(par.dat$v*par.dat$Tw*(par.dat$m^2))
    rstar = (par.dat$c*par.dat$g0/par.dat$m) + (sqrt((par.dat$m^2)*(par.dat$c^2)*par.dat$g*par.dat$g0*par.dat$mu*par.dat$Tw*par.dat$Tv*(par.dat$v*par.dat$Tw+par.dat$g*par.dat$w*par.dat$Tv)))/((par.dat$v*(par.dat$m^2)*par.dat$Tw) +(par.dat$g*par.dat$w*(par.dat$m^2)*par.dat$Tv))
    
    
    return(rstar)
  }
  
  
}
get.betastar.general <- function(par.dat, tol.type){
  
  rstar = get.rstar.general(par.dat = par.dat, tol.type = tol.type)
  
  Lstar = rstar/par.dat$c
  Vstar = ((par.dat$m*rstar)-(par.dat$c*par.dat$g0))/(par.dat$g*par.dat$c)
  
  Lstar[Lstar<0] <- 0
  Vstar[Vstar<0] <- 0
  
  betastar = (par.dat$zeta*Vstar)
  
  return(betastar)
  
}
get.alphastar.general <- function(par.dat, tol.type){
  
  
  if(tol.type=="complete-tolerance"){
    
    
    alphastar = (sqrt((par.dat$c^2)*par.dat$g*par.dat$g0*(par.dat$m^2)*par.dat$mu*(par.dat$v-par.dat$Tv+par.dat$g*(par.dat$w-par.dat$Tw))))/(par.dat$c*par.dat$g*par.dat$m)
    
    return(alphastar)
    
  }else if(tol.type=="constant-tolerance"){
    
    
    #alphastar = ((par.dat$g*par.dat$m*par.dat$w*par.dat$Tv + sqrt(par.dat$c*par.dat$g*par.dat$v*(par.dat$m^2)*par.dat$Tv*par.dat$Tw*(par.dat$m*par.dat$mu*par.dat$Tw+par.dat$g0*par.dat$w)))*(par.dat$c*par.dat$v*par.dat$m*par.dat$Tw*par.dat$g0 + sqrt(par.dat$c*par.dat$g*par.dat$v*(par.dat$m^2)*par.dat$Tv*par.dat$Tw*(par.dat$m*par.dat$mu*par.dat$Tw+par.dat$g0*par.dat$w))))/(par.dat$c*par.dat$g*par.dat$v*(par.dat$m^3)*(par.dat$Tw^2)*par.dat$Tv)
    alphastar = (par.dat$c*par.dat$g0*par.dat$m*par.dat$mu*(par.dat$Tw*par.dat$v+par.dat$g*par.dat$Tv*par.dat$w))/sqrt((par.dat$c^2)*par.dat$g*par.dat$g0*(par.dat$m^2)*par.dat$mu*par.dat$Tv*par.dat$Tw*(par.dat$Tw*par.dat$v+par.dat$g*par.dat$Tv*par.dat$w))
    
    return(alphastar)
  }
  
}
get.prevalence <- function(par.dat, tol.type){
  
  rstar = get.rstar.general(par.dat = par.dat, tol.type = tol.type)
  
  Lstar = rstar/par.dat$c
  Vstar = ((par.dat$m*rstar)-(par.dat$c*par.dat$g0))/(par.dat$g*par.dat$c)
  
  Lstar[Lstar<0] <- 0
  Vstar[Vstar<0] <- 0
  
  betastar = (par.dat$zeta*Vstar)
  
  
  
  
  alphastar = get.alphastar.general(par.dat = par.dat, tol.type = tol.type)
  
  # Istar = Sstar = 0
  
  Istar=(sqrt(betastar*(betastar*((par.dat$b-par.dat$mu-alphastar)^2)+4*par.dat$q*alphastar*(par.dat$mu+ alphastar))) + betastar*(par.dat$b-par.dat$mu-alphastar)-(2*par.dat$q*(par.dat$mu+alphastar)))/(2*betastar*par.dat$q)
  Sstar= (par.dat$mu + alphastar)/betastar
  
  Istar[Istar<0] <- 0
  Sstar[Sstar<0] <- 0
  
  if(!is.na(Istar) & Sstar<Inf){
    
    
    prevalence = NA
    if(Istar >0 | Sstar > 0){
      prevalence = Istar/(Sstar+Istar)  
    }
    
    
  }else{
    prevalence = 0
  }
  
  
  return(prevalence)
  
  
}
get.absInfmort <- function(par.dat, tol.type){
  
  rstar = get.rstar.general(par.dat = par.dat, tol.type = tol.type)
  
  Lstar = rstar/par.dat$c
  Vstar = ((par.dat$m*rstar)-(par.dat$c*par.dat$g0))/(par.dat$g*par.dat$c)
  
  
  Lstar[Lstar<0] <- 0
  Vstar[Vstar<0] <- 0
  
  betastar = (par.dat$zeta*Vstar)
  
  
  alphastar = get.alphastar.general(par.dat = par.dat, tol.type = tol.type)
  
  
  Istar=(sqrt(betastar*(betastar*((par.dat$b-par.dat$mu-alphastar)^2)+4*par.dat$q*alphastar*(par.dat$mu+ alphastar))) + betastar*(par.dat$b-par.dat$mu-alphastar)-(2*par.dat$q*(par.dat$mu+alphastar)))/(2*betastar*par.dat$q)
  Istar[Istar<0] <- 0
  
  if(!is.na(Istar)){
    
    
    absMort = alphastar*Istar
    
  }else{
    absMort = 0
  }
  
  
  
  return(absMort)
  
}
get.relInfmort <- function(par.dat, tol.type){
  
  rstar = get.rstar.general(par.dat = par.dat, tol.type = tol.type)
  
  Lstar = rstar/par.dat$c
  Vstar = ((par.dat$m*rstar)-(par.dat$c*par.dat$g0))/(par.dat$g*par.dat$c)
  
  Lstar[Lstar<0] <- 0
  Vstar[Vstar<0] <- 0
  
  betastar = (par.dat$zeta*Vstar)
  
  
  alphastar = get.alphastar.general(par.dat = par.dat, tol.type = tol.type)
  
  
  Istar=(sqrt(betastar*(betastar*((par.dat$b-par.dat$mu-alphastar)^2)+4*par.dat$q*alphastar*(par.dat$mu+ alphastar))) + betastar*(par.dat$b-par.dat$mu-alphastar)-(2*par.dat$q*(par.dat$mu+alphastar)))/(2*betastar*par.dat$q)
  Istar[Istar<0] <- 0
  
  if(!is.na(Istar)){
    
    
    absMort = alphastar*Istar
    
  }else{
    absMort = 0
  }
  
  
  Sstar= (par.dat$mu + alphastar)/betastar

  Sstar[Sstar<0] <- 0
  
  if (Sstar<Inf){
    Nstar = Istar+Sstar
    relInfMort = absMort/Nstar  
  }else{
    relInfMort = 0
  }
  
  
  
  
  return(relInfMort)
  
}
make.curve <- function(par.dat, tolerance.vect, var.vect, var, tol.shape, tol.type, outcome){
  
  # set up an empty matrix to hold your outputs
  M = matrix(0,nrow=length(var.vect), ncol=length(tolerance.vect))
  
  
  if(tol.shape=="both"){
    
    if(outcome=="rstar"){
      
      # need to loop through the mutants for every value of resident
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.rstar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
      
      
    }else if(outcome=="alphastar"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.alphastar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
      
    }else if(outcome=="betastar"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.betastar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }else if(outcome=="prevalence"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.prevalence(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }else if(outcome=="absInfMort"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.absInfmort(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }
    
  }else if(tol.shape=="Tv"){
    
    if(tol.type=="constant-tolerance"){
      par.dat["Tw"] <- 1
    }else if (tol.type=="complete-tolerance"){
      par.dat["Tw"] <-0
    }
    
    if(outcome=="rstar"){
      
      # need to loop through the mutants for every value of resident
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          # par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.rstar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
      
      
    }else if(outcome=="alphastar"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          # par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.alphastar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
      
    }else if(outcome=="betastar"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.betastar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }else if(outcome=="prevalence"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          # par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.prevalence(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }else if(outcome=="absInfMort"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          # par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.absInfmort(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }
    
  }else if(tol.shape=="Tw"){
    
    
    if(tol.type=="constant-tolerance"){
      par.dat["Tv"] <- 1
    }else if (tol.type=="complete-tolerance"){
      par.dat["Tv"] <-0
    }
    
    if(outcome=="rstar"){
      
      # need to loop through the mutants for every value of resident
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          # par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.rstar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
      
      
    }else if(outcome=="alphastar"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          # par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.alphastar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
      
    }else if(outcome=="betastar"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] =(get.betastar.general(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }else if(outcome=="prevalence"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          # par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.prevalence(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }else if(outcome=="absInfMort"){
      
      for (j in 1:length(tolerance.vect)){ # this is the column (loops through the resident)
        for(i in 1:length(var.vect)){ # this is the row (loops through the mutant)
          par.dat[var] <- var.vect[i]
          par.dat["Tw"] <- tolerance.vect[j]
          # par.dat["Tv"] <- tolerance.vect[j]
          M[i,j] = (get.absInfmort(par.dat = par.dat, tol.type = tol.type))
          # this fills in a column of the matrix
        }}
    }
    
    
  }

  
  colnames(M) = tolerance.vect
  rownames(M) = var.vect
  M.dat <- melt(M)
  names(M.dat) <- c("variable_par", "tolerance_par", "value")
  M.dat$outcome <- outcome
  M.dat$tolerance_type <- tol.type
  M.dat$tolerance_shape <- tol.shape
  M.dat$variable <- var
  return(M.dat)
  
}
make.all.curves <- function(par.dat, var, tol.type, tol.shape){
  
  if(tol.type=="constant-tolerance"){
    # just has to be above 1
    # tolerance.vect = seq(1,2, length.out = 100)
    tolerance.vect = c(1.1,1.5,1.9)
  }else if(tol.type=="complete-tolerance"){
    # has to be less than v and w - theses must be set equal
    # tolerance.vect = seq(.0001,.001, length.out = 100)
    tolerance.vect = c(.2,.5,.9)
  }
  if(var=="mu"){
    var.vect = seq((1/(80*365)), (1/(10*365)), length.out = 100)# different range for host mortality rates
    #var.vect = seq(0.001, 1, length.out = 1000)
    }else{
    var.vect = seq(0.001, 5, length.out = 1000)
  }
  
  
  
  
  map.r <-  make.curve(par.dat=par.dat,
                        tolerance.vect=tolerance.vect,
                        var.vect=var.vect,
                        var=var,
                        tol.shape=tol.shape,
                        tol.type=tol.type,
                        outcome="rstar")
  
  map.alpha <-   make.curve(par.dat=par.dat,
                           tolerance.vect=tolerance.vect,
                           var.vect=var.vect,
                           var=var,
                           tol.shape=tol.shape,
                           tol.type=tol.type,
                            outcome="alphastar")
  
  map.beta <-   make.curve(par.dat=par.dat,
                            tolerance.vect=tolerance.vect,
                            var.vect=var.vect,
                            var=var,
                            tol.shape=tol.shape,
                            tol.type=tol.type,
                            outcome="betastar")
  
  
  map.prev <-  make.curve(par.dat=par.dat,
                                    tolerance.vect=tolerance.vect,
                                    var.vect=var.vect,
                                    var=var,
                                    tol.shape=tol.shape,
                                    tol.type=tol.type,
                           outcome="prevalence")
  
   map.absInfMort <-   make.curve(par.dat=par.dat,
                                  tolerance.vect=tolerance.vect,
                                  var.vect=var.vect,
                                  var=var,
                                  tol.shape=tol.shape,
                                  tol.type=tol.type,
                                  outcome="absInfMort")
  
  #  map.relInfMort <-   make.curve(par.dat=par.dat,
  #                                 tolerance.vect=tolerance.vect,
  #                                 var.vect=var.vect,
  #                                 var=var,
  #                                 tol.shape=tol.shape,
  #                                 tol.type=tol.type,
  #                                 outcome="relInfMort")
  #  
  
  map.combined  <- rbind(map.r, 
                         map.alpha, 
                         map.beta,
                         map.prev, 
                         map.absInfMort)
  
  
  
  return(map.combined)
  
}

combine.all.curves <- function(par.dat){
  
  # now make all maps for each par and combine (c,g,d,gamma,rho,mu)
  c.maps.1 <- make.all.curves(par.dat = par.dat, var="c", tol.type = "constant-tolerance", tol.shape = "Tw")
  c.maps.2 <- make.all.curves(par.dat = par.dat, var="c", tol.type = "complete-tolerance", tol.shape = "Tw")
  
  g.maps.1 <- make.all.curves(par.dat = par.dat, var="g", tol.type = "constant-tolerance", tol.shape = "Tw")
  g.maps.2 <- make.all.curves(par.dat = par.dat, var="g", tol.type = "complete-tolerance", tol.shape = "Tw")
  
  g0.maps.1 <- make.all.curves(par.dat = par.dat, var="g0", tol.type = "constant-tolerance", tol.shape = "Tw")
  g0.maps.2 <- make.all.curves(par.dat = par.dat, var="g0", tol.type = "complete-tolerance", tol.shape = "Tw")
  
  m.maps.1 <- make.all.curves(par.dat = par.dat, var="m", tol.type = "constant-tolerance", tol.shape = "Tw")
  m.maps.2 <- make.all.curves(par.dat = par.dat, var="m", tol.type = "complete-tolerance", tol.shape = "Tw")
  
  mu.maps.1 <- make.all.curves(par.dat = par.dat, var="mu", tol.type = "constant-tolerance", tol.shape = "Tw")
  mu.maps.2 <- make.all.curves(par.dat = par.dat, var="mu", tol.type = "complete-tolerance", tol.shape = "Tw")
  
  
  c.maps.3 <- make.all.curves(par.dat = par.dat, var="c", tol.type = "constant-tolerance", tol.shape = "Tv")
  c.maps.4 <- make.all.curves(par.dat = par.dat, var="c", tol.type = "complete-tolerance", tol.shape = "Tv")
  
  g.maps.3 <- make.all.curves(par.dat = par.dat, var="g", tol.type = "constant-tolerance", tol.shape = "Tv")
  g.maps.4 <- make.all.curves(par.dat = par.dat, var="g", tol.type = "complete-tolerance", tol.shape = "Tv")
  
  g0.maps.3 <- make.all.curves(par.dat = par.dat, var="g0", tol.type = "constant-tolerance", tol.shape = "Tv")
  g0.maps.4 <- make.all.curves(par.dat = par.dat, var="g0", tol.type = "complete-tolerance", tol.shape = "Tv")
  
  m.maps.3 <- make.all.curves(par.dat = par.dat, var="m", tol.type = "constant-tolerance", tol.shape = "Tv")
  m.maps.4 <- make.all.curves(par.dat = par.dat, var="m", tol.type = "complete-tolerance", tol.shape = "Tv")
  
  
  mu.maps.3 <- make.all.curves(par.dat = par.dat, var="mu", tol.type = "constant-tolerance", tol.shape = "Tv")
  mu.maps.4 <- make.all.curves(par.dat = par.dat, var="mu", tol.type = "complete-tolerance", tol.shape = "Tv")
  
  
  
 
  
  out.dat <- rbind(c.maps.1, c.maps.2, c.maps.3, c.maps.4, 
                   g.maps.1, g.maps.2, g.maps.3, g.maps.4, 
                   g0.maps.1, g0.maps.2, g0.maps.3, g0.maps.4, 
                   m.maps.1, m.maps.2, m.maps.3, m.maps.4, 
                   mu.maps.1, mu.maps.2, mu.maps.3, mu.maps.4)
  
  
  
  return(out.dat)
  
}

par.dat <- list(b= .2,# births, per capita per day
                q = .0002,# crowding term, 
                c=.5, # consumption of virus by lymphocytes, per contact, per day
                m=1/(21), # natural mortality of lymphocytes (per day)
                g=.9, # growth of lymphocytes in response to virus, per day
                g0=.3, # rate of constitutive lymphocyte growth, per day
                zeta=.2, # max transmission rate
                v=1, # intrinsic virus virulence
                w=1, # natural damage from immunopathology
                Tv=.005,# tolerance of virus virulence; these get overwritten so values don't matter
                Tw=.005,# tolerance of immunopathology; these get overwritten so values don't matter
                mu = 1/(20*365))# host background death rate, in days

# The range of vectors to plot
out.curves <- combine.all.curves(par.dat = par.dat)
head(out.curves)


out.curves$tolerance <- paste(out.curves$tolerance_shape, out.curves$tolerance_par, sep="=")



out.curves$tolerance <- factor(out.curves$tolerance, levels=c("Tw=1.9", 
                                                              "Tw=1.5", 
                                                              "Tw=1.1",
                                                              "Tv=1.9",
                                                              "Tv=1.5",
                                                              "Tv=1.1",
                                                              "Tw=0.9",
                                                              "Tw=0.5",
                                                              "Tw=0.2",
                                                              "Tv=0.9", 
                                                              "Tv=0.5", 
                                                              "Tv=0.2"))


out.curves$variable <- factor(out.curves$variable, levels=c("mu", "g0", "g", "c", "m"))

names(out.curves)[5:6] <- c("tolerance_shape", "tolerance_type")

# Then, build
save(out.curves, file= paste0(homewd, subwd, "/out.curves.final.2022.Rdata"))



# And make the heatmaps
make.map <- function(par.dat, Tw.vect, Tv.vect,  tol.shape,  outcome){
  
  # set up an empty matrix to hold your outputs
  M = matrix(0,nrow=length(Tv.vect), ncol=length(Tw.vect))
  
  
  if(outcome=="rstar"){
    
    # need to loop through the mutants for every value of resident
    for (j in 1:length(Tw.vect)){ # this is the column (loops through the resident)
      for(i in 1:length(Tv.vect)){ # this is the row (loops through the mutant)
        par.dat["Tw"] <- Tw.vect[j]
        par.dat["Tv"] <- Tv.vect[i]
        #print(i)
        #print(j)
        M[i,j] = (get.rstar.general(par.dat = par.dat, tol.type = tol.shape))
        #print(M[i,j])
        # this fills in a column of the matrix
      }}
    
    
  }else if(outcome=="alphastar"){
    
    for (j in 1:length(Tw.vect)){ # this is the column (loops through the resident)
      for(i in 1:length(Tv.vect)){ # this is the row (loops through the mutant)
        par.dat["Tw"] <- Tw.vect[j]
        par.dat["Tv"] <- Tv.vect[i]
        M[i,j] =(get.alphastar.general(par.dat = par.dat, tol.type = tol.shape))
        # this fills in a column of the matrix
      }}
    
  }else if(outcome=="betastar"){
    
    for (j in 1:length(Tw.vect)){ # this is the column (loops through the resident)
      for(i in 1:length(Tv.vect)){ # this is the row (loops through the mutant)
        par.dat["Tw"] <- Tw.vect[j]
        par.dat["Tv"] <- Tv.vect[i]
        M[i,j] =(get.betastar.general(par.dat = par.dat, tol.type = tol.shape))
        # this fills in a column of the matrix
      }}
    
    
  }else if(outcome=="prevalence"){
    
    for (j in 1:length(Tw.vect)){ # this is the column (loops through the resident)
      for(i in 1:length(Tv.vect)){ # this is the row (loops through the mutant)
        par.dat["Tw"] <- Tw.vect[j]
        par.dat["Tv"] <- Tv.vect[i]
        M[i,j] = (get.prevalence(par.dat = par.dat, tol.type = tol.shape))
        # this fills in a column of the matrix
      }}
  }else if(outcome=="absInfMort"){
    
    for (j in 1:length(Tw.vect)){ # this is the column (loops through the resident)
      for(i in 1:length(Tv.vect)){ # this is the row (loops through the mutant)
        par.dat["Tw"] <- Tw.vect[j]
        par.dat["Tv"] <- Tv.vect[i]
        M[i,j] = (get.absInfmort(par.dat = par.dat, tol.type = tol.shape))
        # this fills in a column of the matrix
      }}
  }
  
  
  
  colnames(M) = Tw.vect
  rownames(M) = Tv.vect
  M.dat <- melt(M)
  names(M.dat) <- c("Tv", "Tw", "value")
  M.dat$outcome <- outcome
  M.dat$tolerance_shape <- tol.shape
  return(M.dat)
  
}
make.all.maps <- function(par.dat, tol.type){
  
  if(tol.type=="constant-tolerance"){
    # just has to be above 1
    Tw.vect = seq(1,2,length.out = 100)
    Tv.vect= Tw.vect
    
  }else if(tol.type=="complete-tolerance"){
    # has to be less than v and w - theses must be set equal
    # tolerance.vect = seq(.0001,.001, length.out = 100)
    Tw.vect = seq(0,.999,length.out = 100)
    Tv.vect= Tw.vect
  }
  
  map.r <-  make.map(par.dat=par.dat,
                     Tw.vect = Tw.vect,
                     Tv.vect = Tv.vect,
                     tol.shape=tol.type,#type
                     outcome="rstar")
  
  map.alpha <-   make.map(par.dat=par.dat,
                          Tw.vect = Tw.vect,
                          Tv.vect = Tv.vect,
                          tol.shape=tol.type,
                          outcome="alphastar")
  
  map.beta <-   make.map(par.dat=par.dat,
                          Tw.vect = Tw.vect,
                          Tv.vect = Tv.vect,
                          tol.shape=tol.type,
                          outcome="betastar")
  
  
  map.prev <-  make.map(par.dat=par.dat,
                        Tw.vect = Tw.vect,
                        Tv.vect = Tv.vect,
                        tol.shape=tol.type,
                        outcome="prevalence")
  
  map.absInfMort <-   make.map(par.dat=par.dat,
                               Tw.vect = Tw.vect,
                               Tv.vect = Tv.vect,
                               tol.shape=tol.type,
                               outcome="absInfMort")
  
  
  map.combined  <- rbind(map.r, 
                         map.alpha,
                         map.beta,
                         map.prev, 
                         map.absInfMort)
  
  
  
  return(map.combined)
  
}
   
# Constant
map.constant <- make.all.maps(par.dat, tol.type = "constant-tolerance")
# Complete
map.complete <- make.all.maps(par.dat, tol.type = "complete-tolerance")
tol.heatmap <- rbind(map.constant, map.complete)
save(tol.heatmap, file = paste0(homewd, subwd, "/tol.heatmap.final.2022.Rdata"))


