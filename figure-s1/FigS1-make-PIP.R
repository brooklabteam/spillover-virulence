library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)


rm(list=ls())

################################################################################
################################################################################


#setwd
homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-s1"
setwd(paste0(homewd, subwd))

# Functions to give you +/- outcome of eigenvalue at different parameters
dot.eigen.out.constant.tol <- function(r_v, r_vm, zeta, v, mu,   g, m, g0, c, w, Tv, Tw){
  
  beta_v = zeta*(m*r_v - c*g0)/(c*g*r_v)
  alpha_v = (r_v*(m*r_v-c*g0))/(c*g*Tv) + (w*(m*r_v-c*g0))/(c*Tw)
  
  beta_vm = zeta*(m*r_vm - c*g0)/(c*g*r_vm)
  alpha_vm = (r_v*(m*r_vm-c*g0))/(c*g*Tv) + (w*(m*r_vm-c*g0))/(c*Tw)
  
  out.eigen = ((mu+alpha_v)/beta_v)*(beta_vm)-mu - alpha_vm
  
  return(out.eigen)
}
dot.eigen.out.complete.tol <- function(r_v, r_vm, zeta, v, mu, rho, g, g0, m, k, c, w, Tv, Tw){
  
  beta_v = zeta*(m*r_v - c*g0)/(c*g*r_v)
  alpha_v = ((r_v*m-c*g0)*(v-Tv))/(c*g) + ((m*r_v-c*g0)*(w-Tw))/c
  
  beta_vm = zeta*(m*r_vm - c*g0)/(c*g*r_vm)
  alpha_vm = ((r_vm*m-c*g0)*(v-Tv))/(c*g) + ((m*r_vm-c*g0)*(w-Tw))/c
  
  out.eigen = ((mu+alpha_v)/beta_v)*(beta_vm)-mu - alpha_vm
  
  return(out.eigen)
}


# This script makes the PIP plot for a given combination of r_v and r_vm vectors

r_v= seq(3.18,3.5,length=100)
r_vm = rev(r_v)

# Functions to make PIPs
make.four.PIP <- function(r_v, r_vm, zeta, v,   mu, g, m, g0, c, w, Tv.lo.com, Tw.lo.com, Tv.hi.com, Tw.hi.com, Tv.lo.cons, Tw.lo.cons, Tv.hi.cons, Tw.hi.cons, filename){
  
  #set up an empty matrix to hold your outputs
  M.hi.cons = M.lo.cons= matrix(0,nrow=length(r_v), ncol=length(r_vm))
  M.hi.com = M.lo.com= matrix(0,nrow=length(r_v), ncol=length(r_vm))
  
  #if (tol.type=="constant-tolerance"){
    
    #need to loop through the mutants for every value of resident
    for (j in 1:length(r_v)){ #this is the column (loops through the resident)
      for(i in 1:length(r_vm)){ #this is the row (loops through the mutabt)
        M.hi.cons[i,j] =(dot.eigen.out.constant.tol(r_v= r_v[j],
                                               r_vm= r_vm[i],
                                               zeta=zeta,
                                               v=v,
                                               g=g,
                                               g0=g0,
                                               mu=mu,
                                               m=m,
                                               c=c,
                                               w=w,
                                               Tv=Tv.hi.cons,
                                               Tw=Tw.hi.cons))#induced vs. constitutive
        #this fills in a column of the matrix
      }}
    
    #need to loop through the mutants for every value of resident
    for (j in 1:length(r_v)){ #this is the column (loops through the resident)
      for(i in 1:length(r_vm)){ #this is the row (loops through the mutabt)
        M.lo.cons[i,j] =(dot.eigen.out.constant.tol(r_v= r_v[j],
                                               r_vm= r_vm[i],
                                               zeta=zeta,
                                               v=v,
                                               g=g,
                                               g0=g0,
                                               mu=mu,
                                               m=m,
                                               c=c,
                                               w=w,
                                               Tv=Tv.lo.cons,
                                               Tw=Tw.lo.cons))#induced vs. constitutive
        #this fills in a column of the matrix
      }}
    
  #}
  
  
 # if (tol.type=="complete-tolerance"){
    
    #need to loop through the mutants for every value of resident
    for (j in 1:length(r_v)){ #this is the column (loops through the resident)
      for(i in 1:length(r_vm)){ #this is the row (loops through the mutabt)
        M.lo.com[i,j] =(dot.eigen.out.complete.tol(r_v= r_v[j],
                                               r_vm= r_vm[i],
                                               zeta=zeta,
                                               v=v,
                                               g=g,
                                               g0=g0,
                                               m=m,
                                               mu=mu,
                                               c=c,
                                               w=w,
                                               Tv=Tv.lo.com,
                                               Tw=Tw.lo.com))#induced vs. constitutive
        #this fills in a column of the matrix
      }}
    
    for (j in 1:length(r_v)){ #this is the column (loops through the resident)
      for(i in 1:length(r_vm)){ #this is the row (loops through the mutabt)
        M.hi.com[i,j] =(dot.eigen.out.complete.tol(r_v= r_v[j],
                                               r_vm= r_vm[i],
                                               zeta=zeta,
                                               v=v,
                                               g=g,
                                               g0=g0,
                                               m=m,
                                               mu=mu,
                                               c=c,
                                               w=w,
                                               Tv=Tv.hi.com,
                                               Tw=Tw.hi.com))#induced vs. constitutive
        #this fills in a column of the matrix
      }}
  #}
  
  M.hi.com = M.hi.com>0 #makes all places where invasion is possible into "TRUE" (value >0) and all others "FALSE"
  M.hi.cons = M.hi.cons >0
  M.lo.com = M.lo.com >0
  M.lo.cons = M.lo.cons >0
  rownames(M.hi.cons) = rownames(M.lo.cons)= r_vm
  rownames(M.hi.com) = rownames(M.lo.com)= r_vm
  colnames(M.hi.cons) = colnames(M.lo.cons) = r_v
  colnames(M.hi.com) = colnames(M.lo.com) = r_v
  M.dat.hi.com <- melt(M.hi.com)
  M.dat.hi.cons <- melt(M.hi.cons)
  M.dat.lo.com <- melt(M.lo.com)
  M.dat.lo.cons <- melt(M.lo.cons)
  names(M.dat.hi.cons) <- names(M.dat.lo.cons) <- c("r_vm", "r_v", "eigen_val")
  names(M.dat.hi.com) <- names(M.dat.lo.com) <- c("r_vm", "r_v", "eigen_val")
  M.dat.hi.com$tol_level <- "high-tolerance"
  M.dat.hi.cons$tol_level <- "high-tolerance"
  M.dat.lo.com$tol_level <- "low-tolerance"
  M.dat.lo.cons$tol_level <- "low-tolerance"
  M.dat.com <- rbind(M.dat.hi.com, M.dat.lo.com)
  M.dat.cons <- rbind(M.dat.hi.cons, M.dat.lo.cons)
  M.dat.cons$tol_form <- "constant-tolerance"
  M.dat.com$tol_form <- "complete-tolerance"
  M.dat <- rbind(M.dat.cons, M.dat.com)
  M.dat$tol_level <- factor(M.dat$tol_level, levels=c("low-tolerance", "high-tolerance"))
  M.dat$tol_form <- factor(M.dat$tol_form, levels=c("complete-tolerance", "constant-tolerance"))
  colz = c('TRUE' = "black", 'FALSE' ="white")

  
  
  p1 <- ggplot(data = M.dat) + geom_point(aes(x=r_v, y=r_vm, color = eigen_val), show.legend = FALSE) + scale_color_manual(values=colz) + 
    facet_grid(tol_form~tol_level) +
    theme_bw() + theme(panel.grid = element_blank(), axis.ticks = element_blank(), 
                       axis.text = element_blank(), strip.text = element_text(size=12),
                       strip.background = element_rect(fill="white"),
                       axis.title = element_text(size=20)) + coord_cartesian(expand=FALSE) +
    ylab(expression("r"["2"])) +xlab(expression("r"["1"]))
  
  
  ggsave(file = filename,
         plot = p1,
         units="mm",  
         width=60, 
         height=50, 
         scale=3, 
         dpi=200)
}


# And make onw with default parameter values:
make.four.PIP(r_v= r_v,
                r_vm= r_vm,
                zeta=.2,
                v=1,
                m=1/(21),
                g=.9,
                g0=.3,
                c=.5,
                mu = 1/(20*365),
                w=1,
                Tv.lo.cons=10,
                Tw.lo.cons=10,
                Tv.lo.com=.5,
                Tw.lo.com=.5,
                Tv.hi.cons=100,
                Tw.hi.cons=100,
                Tv.hi.com=.97,
                Tw.hi.com=.97,
                filename = paste0(homewd,"/supp-figs/FigS1.png"))

