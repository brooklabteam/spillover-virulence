rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(deSolve)
library(scales)

# load the data for the plot
homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-2"
setwd(paste0(homewd, subwd))

load("out.curves.final.2022.Rdata") 
load("tol.heatmap.final.2022.Rdata")

# plotting functions
par.plots <- function(dat.var, tol.shape, sub){
  
  # and add labels for strip labels
  dat.var$label= NA
  
  dat.var$label[dat.var$variable=="c"] <- "c~', virus consumption'"
  dat.var$label[dat.var$variable=="g"] <- "g~', leukocyte activation'"
  dat.var$label[dat.var$variable=="m"] <- "m~', leukocyte mortality'"
  dat.var$label[dat.var$variable=="mu"] <- "mu~', host background mortality'"
  dat.var$label[dat.var$variable=="g0"] <- "g[0]~', constitutive immunity'"
  
  
  dat.var$label = factor(dat.var$label, levels = c("mu~', host background mortality'",
                                                   "g[0]~', constitutive immunity'",
                                                   "g~', leukocyte activation'",
                                                   "c~', virus consumption'",
                                                   "m~', leukocyte mortality'"))
  
  dat.var = subset(dat.var, tolerance_shape == tol.shape)
  if(tol.shape=="constant-tolerance"){
    dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=1.9", "Tv=1.9", "Tw=1.5", "Tv=1.5","Tw=1.1", "Tv=1.1"))
  }else{
    dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=0.9", "Tv=0.9", "Tw=0.5", "Tv=0.5","Tw=0.2", "Tv=0.2"))
  }
  
  #dat.var = subset(dat.var, variable_par <=1)
  dat.var$value[dat.var$value<0]<-0
  
  #dat.var = subset(dat.var, tolerance_shape==tol.shape)
  
  # just rstar for now
  r.dat = subset(dat.var, outcome=="rstar" & tolerance_shape==tol.shape)
  a.dat = subset(dat.var, outcome=="alphastar" & tolerance_shape==tol.shape)
  #p.dat = subset(dat.var, outcome=="prevalence" & tolerance_shape==tol.shape)
  #inf.dat = subset(dat.var, outcome=="absInfMort" & tolerance_shape==tol.shape)
  beta.dat = subset(dat.var, outcome=="betastar" & tolerance_shape==tol.shape)
  
  
  if (tol.shape=="constant-tolerance"){
    
    colz = c("Tw=1.9" = "dodgerblue4", 
             "Tv=1.9" = "firebrick4",
             "Tw=1.5"= "dodgerblue3", 
             "Tv=1.5" = "firebrick3",  
             "Tw=1.1" = "dodgerblue1",
             "Tv=1.1" ="firebrick1")
    
  }else if (tol.shape=="complete-tolerance"){
    
    
    colz = c("Tw=0.9"  = "dodgerblue4",
             "Tv=0.9"=  "firebrick4",
             "Tw=0.5"= "dodgerblue3",
             "Tv=0.5"=  "firebrick3",
             "Tw=0.2"= "dodgerblue1", 
             "Tv=0.2"=  "firebrick1")
    
  }
  
  
  
  
  pA <- ggplot(data=r.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance)) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote("r"^"*"~", virus growth")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size=15),
          axis.title.y = element_text(size=18), strip.background = element_rect(fill="white"), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(.5,.1,0,1.3), "cm"),
          legend.position = c(.1,.8), legend.title = element_blank(), legend.text = element_text(size=10)) + 
          scale_color_manual(values=colz) + guides(color=guide_legend(nrow=2)) + coord_cartesian(ylim = c(0,5))  
 #  print(pA)
  
  
  pB <- ggplot(data=beta.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(beta^"*"~", transmission")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(0,.1,0,.4), "cm")) + coord_cartesian(ylim = c(0,.005)) +
    scale_color_manual(values=colz)
  # print(pB)
 
  
  # and virulence
  
  
  pC <- ggplot(data=a.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(alpha^"*"~", virulence")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text = element_text(size=14),
          plot.margin = unit(c(0,.1,1,.7), "cm")) + coord_cartesian(ylim = c(0,.05)) +
    scale_color_manual(values=colz) + scale_x_continuous(labels = label_comma())
  # print(pC)
  
  #and transmission
  
  
  
  
  out.all<- cowplot:: plot_grid(pA,pB,pC, nrow=3,ncol=1, rel_heights = c(1,.9,1.1))
  # print(out.all)
  
  # return, and plot side by side with the heatmap
  return(out.all)
  
  #  
  #  ggsave(file =filename,
  #         units="mm",  
  #         width=150, 
  #         height=100, 
  #         scale=3, 
  #         dpi=200)
  #  
}
heat.plots <- function(dat.tol, tol.shape){
  
  #dat.tol = subset(dat.tol, Tv<=2 & Tw<=2)
  r.dat = subset(dat.tol, outcome=="rstar" & tolerance_shape==tol.shape)
  a.dat = subset(dat.tol, outcome=="alphastar" & tolerance_shape==tol.shape)
  #p.dat = subset(dat.tol, outcome=="prevalence" & tolerance_shape==tol.shape)
  #inf.dat = subset(dat.tol, outcome=="absInfMort" & tolerance_shape==tol.shape)
  beta.dat = subset(dat.tol, outcome=="betastar" & tolerance_shape==tol.shape)
  
  
  pA <- ggplot(r.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
        scale_fill_gradient(low="yellow", high="red", name=bquote('r'^'*')) +
        theme_bw() + 
        # xlab(bquote('T'[w]~'tolerance of immunopathology'))+
        ylab(bquote('T'[v]~'tolerance of viral pathology'))+
        theme(panel.grid = element_blank(),
              axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1.4,.5,0,1), "cm"))
 #  print(pA)
  
  
  
  pB <- ggplot(beta.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(beta^'*')) +
    theme_bw() + 
    # xlab(bquote('T'[w]~'tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.3,0,1), "cm"))
  
  pC <- ggplot(a.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(alpha^'*')) +
    theme_bw() + 
     xlab(bquote('T'[w]~'tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size =16, color="navy"),
          axis.text = element_text(size=14, color="black"),
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.4,.1,1), "cm"))
  # # print(pB)
  # 
  # 
  
  
  
  p.heat <- cowplot::plot_grid(pA,pB,pC, ncol = 1, nrow = 3, rel_heights = c(1,.9,1.1))
  
  
  p.heat2 <- p.heat + 
             annotate("text", color="firebrick",
                      label=bquote("T"[v]~", tolerance of virus pathology"), 
                      x=.05,y=.55, angle=90, size=6)
             
 
  
  return(p.heat2)
}
plot.join <- function(dat.var, dat.tol, filename, tol.shape, sub){
  
  plot1 <- par.plots(dat.var = dat.var,
                     tol.shape = tol.shape,
                     sub = sub)
  
  # and plot 2
  plot2 <- heat.plots(dat.tol=dat.tol,
                      tol.shape = tol.shape)
  
  
  both.plot <- cowplot::plot_grid(plot1, plot2, ncol=2, nrow = 1, rel_widths = c(6.2,2))
  
  ggsave(file =filename,
         plot = both.plot,
         units="mm",  
         width=170, 
         height=90, 
         scale=3, 
         dpi=200)
          
}

# constant tolerance plot for the main text
plot.join(dat.var = subset(out.curves, variable=="mu" | variable=="g0" | variable =="c"| variable=="g"| variable=="m"),
                       tol.shape= "constant-tolerance",
                       dat.tol=tol.heatmap,
                       sub=TRUE,
                       filename=paste0(homewd,"main-figs/Fig2.png"))

#and pdf for submission
# constant tolerance plot for the main text
plot.join(dat.var = subset(out.curves, variable=="mu" | variable=="g0" | variable =="c"| variable=="g"| variable=="m"),
          tol.shape= "constant-tolerance",
          dat.tol=tol.heatmap,
          sub=TRUE,
          filename=paste0(homewd,"main-figs/Fig2.pdf"))




# and the complete tolerance version for the supplement
heat.plots.complete <- function(dat.tol, tol.shape){
  
  
  r.dat = subset(dat.tol, outcome=="rstar" & tolerance_shape==tol.shape)
  a.dat = subset(dat.tol, outcome=="alphastar" & tolerance_shape==tol.shape)
  #p.dat = subset(dat.tol, outcome=="prevalence" & tolerance_shape==tol.shape)
  #inf.dat = subset(dat.tol, outcome=="absInfMort" & tolerance_shape==tol.shape)
  beta.dat = subset(dat.tol, outcome=="betastar" & tolerance_shape==tol.shape)
  

  
  
  pA <- ggplot(r.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote('r'^'*')) +
    theme_bw() + 
    xlab(bquote('T'[w]~'tolerance of immunopathology'))+
    ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1.4,.9,0,1), "cm"))
  #  print(pA)
  
  
  pB <- ggplot(beta.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(beta^'*')) +
    theme_bw() + 
    # xlab(bquote('T'[w]~'tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.7,0,1), "cm"))
  
  
  
  pC <- ggplot(a.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red",  name=bquote(alpha^'*')) +
    theme_bw() + 
    xlab(bquote('T'[w]~'tolerance of immunopathology'))+
    ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size =16, color="navy"),
          axis.text = element_text(size=14, color="black"),
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.4,.1,1), "cm"))
  # # print(pB)
  # 
  # 
  
  
  
  p.heat <- cowplot::plot_grid(pA,pB,pC, ncol = 1, nrow = 3, rel_heights = c(1,.9,1.1))
  
  
  p.heat2 <- p.heat + 
    annotate("text", color="firebrick",
             label=bquote("T"[v]~", tolerance of virus pathology"), 
             x=.05,y=.55, angle=90, size=6)
  
  
  
  return(p.heat2)
}
par.plots.complete <- function(dat.var, tol.shape, sub, plot.prev){
  
  # and add labels for strip labels
  dat.var$label= NA
  
  dat.var$label[dat.var$variable=="c"] <- "c~', virus consumption'"
  dat.var$label[dat.var$variable=="g"] <- "g~', leukocyte activation'"
  dat.var$label[dat.var$variable=="m"] <- "m~', leukocyte mortality'"
  dat.var$label[dat.var$variable=="mu"] <- "mu~', host background mortality'"
  dat.var$label[dat.var$variable=="g0"] <- "g[0]~', constitutive immunity'"
  
  
  dat.var$label = factor(dat.var$label, levels = c("mu~', host background mortality'",
                                                   "g[0]~', constitutive immunity'",
                                                   "g~', leukocyte activation'",
                                                   "c~', virus consumption'",
                                                   "m~', leukocyte mortality'"))
  
  dat.var = subset(dat.var, tolerance_shape == tol.shape)
  if(tol.shape=="constant-tolerance"){
    dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=1.9", "Tv=1.9", "Tw=1.5", "Tv=1.5","Tw=1.1", "Tv=1.1"))
  }else{
    dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=0.9", "Tv=0.9", "Tw=0.5", "Tv=0.5","Tw=0.2", "Tv=0.2"))
  }
  
  #dat.var = subset(dat.var, variable_par <=1)
  dat.var$value[dat.var$value<0]<-0
  
  #dat.var = subset(dat.var, tolerance_shape==tol.shape)
  
  # just rstar for now
  r.dat = subset(dat.var, outcome=="rstar" & tolerance_shape==tol.shape)
  a.dat = subset(dat.var, outcome=="alphastar" & tolerance_shape==tol.shape)
  #p.dat = subset(dat.var, outcome=="prevalence" & tolerance_shape==tol.shape)
  #inf.dat = subset(dat.var, outcome=="absInfMort" & tolerance_shape==tol.shape)
  beta.dat = subset(dat.var, outcome=="betastar" & tolerance_shape==tol.shape)
  
  
  if (tol.shape=="constant-tolerance"){
    
    colz = c("Tw=1.9" = "dodgerblue4", 
             "Tv=1.9" = "firebrick4",
             "Tw=1.5"= "dodgerblue3", 
             "Tv=1.5" = "firebrick3",  
             "Tw=1.1" = "dodgerblue1",
             "Tv=1.1" ="firebrick1")
    
  }else if (tol.shape=="complete-tolerance"){
    
    
    colz = c("Tw=0.9"  = "dodgerblue4",
             "Tv=0.9"=  "firebrick4",
             "Tw=0.5"= "dodgerblue3",
             "Tv=0.5"=  "firebrick3",
             "Tw=0.2"= "dodgerblue1", 
             "Tv=0.2"=  "firebrick1")
    
  }
  
  
  
  pA <- ggplot(data=r.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance)) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote("r"^"*"~", virus growth")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size=15),
          axis.title.y = element_text(size=18), strip.background = element_rect(fill="white"), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(.5,.1,0,1.3), "cm"),
          legend.position = c(.1,.8), legend.title = element_blank(), legend.text = element_text(size=10)) + 
    scale_color_manual(values=colz) + guides(color=guide_legend(nrow=2)) + coord_cartesian(ylim = c(0,5))  
  #  print(pA)
  
  
  pB <- ggplot(data=beta.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(beta^"*"~", transmission")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(0,.1,0,.4), "cm")) + coord_cartesian(ylim = c(0,.005)) +
    scale_color_manual(values=colz)
  # print(pB)
  
  
  
  pC <- ggplot(data=a.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(alpha^"*"~", virulence")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text = element_text(size=14),
          plot.margin = unit(c(0,.1,1,.7), "cm")) + coord_cartesian(ylim = c(0,.05)) +
    scale_color_manual(values=colz) + scale_x_continuous(labels = label_comma())
  # print(pC)
  
  #and transmission
  
  
  
  
  out.all<- cowplot:: plot_grid(pA,pB,pC, nrow=3,ncol=1, rel_heights = c(1,.9,1.1))
  # print(out.all)
  # return, and plot side by side with the heatmap
  return(out.all)
  
  #  
  #  ggsave(file =filename,
  #         units="mm",  
  #         width=150, 
  #         height=100, 
  #         scale=3, 
  #         dpi=200)
  #  
}
plot.join.complete <- function(dat.var, dat.tol, filename, tol.shape, sub){
  plot1 <- par.plots.complete(dat.var = dat.var,
                     tol.shape = tol.shape,
                     sub = sub)
  
  # and plot 2
  plot2 <- heat.plots.complete(dat.tol=dat.tol,
                      tol.shape = tol.shape)
  
  
  both.plot <- cowplot::plot_grid(plot1, plot2, ncol=2, nrow = 1, rel_widths = c(6,2))
  
  ggsave(file =filename,
         plot = both.plot,
         units="mm",  
         width=170, 
         height=90, 
         scale=3, 
         dpi=200)
}

plot.join.complete(dat.var = subset(out.curves, variable=="mu" | variable=="g0" | variable =="c"| variable=="g"| variable=="m"),
          tol.shape= "complete-tolerance",
          dat.tol=tol.heatmap,
          sub=TRUE,
          filename=paste0(homewd,"supp-figs/FigS2.png"))

