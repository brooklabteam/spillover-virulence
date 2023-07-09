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

load("out.curves.final.2023.Rdata") 
load("tol.heatmap.final.2023.Rdata")



# plotting functions
par.plots <- function(dat.var, tol.shape, sub){
  
  # and add labels for strip labels
  dat.var$label= NA
  
  dat.var$label[dat.var$variable=="c"] <- "atop(c[R]~', virus consumption', 'rate in reservoir host')"
  dat.var$label[dat.var$variable=="g"] <- "atop(g[R]~', reservoir', 'leukocyte activation rate')"
  dat.var$label[dat.var$variable=="m"] <- "atop(m[R]~', reservoir', 'leukocyte mortality rate')"
  dat.var$label[dat.var$variable=="mu"] <- "atop(mu[R]~', reservoir host', 'background mortality rate')"
  dat.var$label[dat.var$variable=="g0"] <- "atop(g[0~R]~', magnitude reservoir', 'constitutive immunity')"
  #dat.var$label[dat.var$variable=="Tv_spill"] <- "T[vS]~', spillover host tolerance of direct virus pathology'"
  #dat.var$label[dat.var$variable=="Tw_spill"] <- "T[wS]~', spillover host tolerance of immunopathology'"
  
  
  dat.var$label = factor(dat.var$label, levels = c("atop(mu[R]~', reservoir host', 'background mortality rate')",
                                                   "atop(g[0~R]~', magnitude reservoir', 'constitutive immunity')",
                                                   "atop(g[R]~', reservoir', 'leukocyte activation rate')",
                                                   "atop(c[R]~', virus consumption', 'rate in reservoir host')",
                                                   "atop(m[R]~', reservoir', 'leukocyte mortality rate')"))#,
                                                   #"T[vS]~', spillover host tolerance of direct virus pathology'",
                                                   #"T[wS]~', spillover host tolerance of immunopathology'"))
  
  dat.var = subset(dat.var, tolerance_shape == tol.shape)
  if(tol.shape=="constant-tolerance"){
    dat.var$tolerance <- as.character(dat.var$tolerance)
    dat.var$tolerance[dat.var$tolerance=="Tw=1.9"] <- bquote('T[wR]==1.9')
    dat.var$tolerance[dat.var$tolerance=="Tv=1.9"] <- bquote('T[vR]==1.9')
    dat.var$tolerance[dat.var$tolerance=="Tw=1.5"] <- bquote('T[wR]==1.5')
    dat.var$tolerance[dat.var$tolerance=="Tv=1.5"] <- bquote('T[vR]==1.5')
    dat.var$tolerance[dat.var$tolerance=="Tw=1.1"] <- bquote('T[wR]==1.1')
    dat.var$tolerance[dat.var$tolerance=="Tv=1.1"] <- bquote('T[vR]==1.1')
    
    dat.var$tolerance = factor(dat.var$tolerance, levels=c(bquote('T[wR]==1.9'), bquote('T[vR]==1.9'), bquote('T[wR]==1.5'), bquote('T[vR]==1.5'),bquote('T[wR]==1.1'), bquote('T[vR]==1.1')))
    
    
    
    #dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=1.9", "Tv=1.9", "Tw=1.5", "Tv=1.5","Tw=1.1", "Tv=1.1"))
  }else{
    
    dat.var$tolerance <- as.character(dat.var$tolerance)
    dat.var$tolerance[dat.var$tolerance=="Tw=0.9"] <- bquote('T[wR]==0.9')
    dat.var$tolerance[dat.var$tolerance=="Tv=0.9"] <- bquote('T[vR]==0.9')
    dat.var$tolerance[dat.var$tolerance=="Tw=0.5"] <- bquote('T[wR]==0.5')
    dat.var$tolerance[dat.var$tolerance=="Tv=0.5"] <- bquote('T[vR]==0.5')
    dat.var$tolerance[dat.var$tolerance=="Tw=0.2"] <- bquote('T[wR]==0.2')
    dat.var$tolerance[dat.var$tolerance=="Tv=0.2"] <- bquote('T[vR]==0.2')
    
    dat.var$tolerance = factor(dat.var$tolerance, levels=c(bquote('T[wR]==0.9'), bquote('T[vR]==0.9'), bquote('T[wR]==0.5'), bquote('T[vR]==0.5'),bquote('T[wR]==0.2'), bquote('T[vR]==0.2')))
    
    
    #dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=0.9", "Tv=0.9", "Tw=0.5", "Tv=0.5","Tw=0.2", "Tv=0.2"))
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
  spill.dat = subset(dat.var, outcome=="alphaSpillover" & tolerance_shape==tol.shape)
  
  
  if (tol.shape=="constant-tolerance"){
    
    
    colz = c('T[wR]==1.9' = "dodgerblue4", 
             'T[vR]==1.9' = "firebrick4",
             'T[wR]==1.5'= "dodgerblue3", 
             'T[vR]==1.5' = "firebrick3",  
             'T[wR]==1.1' = "dodgerblue1",
             'T[vR]==1.1' ="firebrick1")
    
    
    #colz = c(bquote('T[w]==1.9') = "dodgerblue4", 
     #        bquote('T[v]==1.9') = "firebrick4",
      #       bquote('T[w]==1.5')= "dodgerblue3", 
       #      bquote('T[v]==1.5') = "firebrick3",  
        #     bquote('T[w]==1.1') = "dodgerblue1",
         #    bquote('T[v]==1.1') ="firebrick1")
    
    
    
    # 
    # colz = c("Tw=1.9" = "dodgerblue4", 
    #          "Tv=1.9" = "firebrick4",
    #          "Tw=1.5"= "dodgerblue3", 
    #          "Tv=1.5" = "firebrick3",  
    #          "Tw=1.1" = "dodgerblue1",
    #          "Tv=1.1" ="firebrick1")
    
  }else if (tol.shape=="complete-tolerance"){
    
    
    #colz = c("Tw=0.9"  = "dodgerblue4",
     #        "Tv=0.9"=  "firebrick4",
      #       "Tw=0.5"= "dodgerblue3",
       #      "Tv=0.5"=  "firebrick3",
        #     "Tw=0.2"= "dodgerblue1", 
         #    "Tv=0.2"=  "firebrick1")
    
    colz = c('T[wR]==0.9' = "dodgerblue4", 
             'T[vR]==0.9' = "firebrick4",
             'T[wR]==0.5'= "dodgerblue3", 
             'T[vR]==0.5' = "firebrick3",  
             'T[wR]==0.2' = "dodgerblue1",
             'T[vR]==0.2' ="firebrick1")
    
  }
  
  
  
  
  pA <- ggplot(data=r.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance)) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(r[R]^"*"~", virus growth")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size=15),
          axis.title.y = element_text(size=18), strip.background = element_rect(fill="white"), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(.5,.2,0,1.3), "cm"),
          legend.position = c(.1,.2), legend.title = element_blank(), legend.text = element_text(size=8)) + 
          scale_color_manual(values=colz, labels = function(x) parse(text = x)) + guides(color=guide_legend(nrow=2)) + coord_cartesian(ylim = c(0,5))  
 #  print(pA)
  
  
  pB <- ggplot(data=beta.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(beta[r[R]]^"*"~", transmission")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(0,.2,0,.4), "cm")) + coord_cartesian(ylim = c(0,.005)) +
    scale_color_manual(values=colz)
  # print(pB)
 
  
  # and virulence
  
  
  pC <- ggplot(data=a.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(alpha[r[R]]^"*"~", virulence")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.2,.0,.7), "cm")) + coord_cartesian(ylim = c(0,.05)) +
    scale_color_manual(values=colz) + scale_x_continuous(labels = label_comma())
  # print(pC)
  
  #and spillover
   
  pD <- ggplot(data=spill.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(alpha[S]~", spillover virulence")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text = element_text(size=14),
          plot.margin = unit(c(0,.2,1,1), "cm")) + coord_cartesian(ylim = c(0,100)) +
    scale_color_manual(values=colz) + scale_x_continuous(labels = label_comma())
  # print(pD)
  
  
  
  
  out.all<- cowplot:: plot_grid(pA,pB,pC, pD, nrow=4,ncol=1, rel_heights = c(1,.9,.9,1.1))
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
  spill.dat = subset(dat.tol, outcome=="alphaSpillover" & tolerance_shape==tol.shape)
  
  
  pA <- ggplot(r.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
        scale_fill_gradient(low="yellow", high="red", name=bquote(r[R]^'*')) +
        theme_bw() + 
        # xlab(bquote('T'[w]~'tolerance of immunopathology'))+
        ylab(bquote('T'[vR]~'tolerance of viral pathology'))+
        theme(panel.grid = element_blank(),
              axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1.4,.5,0,1), "cm"))
 #  print(pA)
  
  
  
  pB <- ggplot(beta.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(beta[r[R]]^'*')) +
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
    scale_fill_gradient(low="yellow", high="red", name=bquote(alpha[r[R]]^'*')) +
    theme_bw() + 
     xlab(bquote('T'[wR]~', tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.3,0,1), "cm"))
  # # print(pB)
  # 
  # 
  pD <- ggplot(beta.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(alpha[s])) +
    theme_bw() + 
    xlab(bquote('T'[wR]~', tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size =16, color="navy"),
          axis.text = element_text(size=14, color="black"),
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.4,.1,1), "cm"))
  # # print(pD)
  # 
  
  
  
  p.heat <- cowplot::plot_grid(pA,pB,pC, pD, ncol = 1, nrow = 4, rel_heights = c(1,.9,.9,1.1))
  
  
  p.heat2 <- p.heat + 
             annotate("text", color="firebrick",
                      label=bquote("T"[vR]~", tolerance of virus pathology"), 
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
         height=110, 
         scale=3, 
         dpi=300)
          
}

unique(out.curves$variable)
# constant tolerance plot for the main text
plot.join(dat.var = subset(out.curves, variable !="Tv_spill" & variable!="Tw_spill"),
                       tol.shape= "constant-tolerance",
                       dat.tol=tol.heatmap,
                       sub=TRUE,
                       filename=paste0(homewd,"main-figs/Fig2.png"))

#and pdf for submission
# constant tolerance plot for the main text
plot.join(dat.var = subset(out.curves, variable !="Tv_spill" & variable!="Tw_spill"),
          tol.shape= "constant-tolerance",
          dat.tol=tol.heatmap,
          sub=TRUE,
          filename=paste0(homewd,"main-figs/Fig2.pdf"))

#and eps for submission
plot.join(dat.var = subset(out.curves, variable !="Tv_spill" & variable!="Tw_spill"),
          tol.shape= "constant-tolerance",
          dat.tol=tol.heatmap,
          sub=TRUE,
          filename=paste0(homewd,"main-figs/Fig2.eps"))


# and the complete tolerance version for the supplement
heat.plots.complete <- function(dat.tol, tol.shape){
  
  
  r.dat = subset(dat.tol, outcome=="rstar" & tolerance_shape==tol.shape)
  a.dat = subset(dat.tol, outcome=="alphastar" & tolerance_shape==tol.shape)
  #p.dat = subset(dat.tol, outcome=="prevalence" & tolerance_shape==tol.shape)
  #inf.dat = subset(dat.tol, outcome=="absInfMort" & tolerance_shape==tol.shape)
  beta.dat = subset(dat.tol, outcome=="betastar" & tolerance_shape==tol.shape)
  spill.dat = subset(dat.tol, outcome=="alphaSpillover" & tolerance_shape==tol.shape)

  
  
  pA <- ggplot(r.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(r[R]^'*'), trans="log10") +
    theme_bw() + 
    xlab(bquote('T'[wR]~', tolerance of immunopathology'))+
    ylab(bquote('T'[vR]~', tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(1.4,1.1,0,1), "cm"))
  #  print(pA)
  
  
  pB <- ggplot(beta.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(beta[r[R]]^'*'), trans="log10") +
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
    scale_fill_gradient(low="yellow", high="red", name=bquote(alpha[r[R]]^'*')) +
    theme_bw() + 
    xlab(bquote('T'[wR]~', tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size=14, color="black"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.5,0,1), "cm"))
  # # print(pB)
  # 
  # 
  pD <- ggplot(beta.dat) + geom_tile(aes(x=Tw, y=Tv, fill=value)) + 
    scale_fill_gradient(low="yellow", high="red", name=bquote(alpha[s]), trans="log10") +
    theme_bw() + 
    xlab(bquote('T'[wR]~', tolerance of immunopathology'))+
    # ylab(bquote('T'[v]~'tolerance of viral pathology'))+
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size =16, color="navy"),
          axis.text = element_text(size=14, color="black"),
          #axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,.6,.1,1), "cm"))
  # # print(pD)
  # 
  
  
  
  p.heat <- cowplot::plot_grid(pA,pB,pC, pD, ncol = 1, nrow = 4, rel_heights = c(1,.9,.9,1.1))
  
  
  p.heat2 <- p.heat + 
    annotate("text", color="firebrick",
             label=bquote("T"[vR]~", tolerance of virus pathology"), 
             x=.05,y=.55, angle=90, size=6)
  
  
  return(p.heat2)
}
par.plots.complete <- function(dat.var, tol.shape, sub, plot.prev){
  
  
  dat.var$label[dat.var$variable=="c"] <- "atop(c[R]~', virus consumption', 'rate in reservoir host')"
  dat.var$label[dat.var$variable=="g"] <- "atop(g[R]~', reservoir', 'leukocyte activation rate')"
  dat.var$label[dat.var$variable=="m"] <- "atop(m[R]~', reservoir', 'leukocyte mortality rate')"
  dat.var$label[dat.var$variable=="mu"] <- "atop(mu[R]~', reservoir host', 'background mortality rate')"
  dat.var$label[dat.var$variable=="g0"] <- "atop(g[0~R]~', magnitude reservoir', 'constitutive immunity')"
  #dat.var$label[dat.var$variable=="Tv_spill"] <- "T[vS]~', spillover host tolerance of direct virus pathology'"
  #dat.var$label[dat.var$variable=="Tw_spill"] <- "T[wS]~', spillover host tolerance of immunopathology'"
  
  
  dat.var$label = factor(dat.var$label, levels = c("atop(mu[R]~', reservoir host', 'background mortality rate')",
                                                   "atop(g[0~R]~', magnitude reservoir', 'constitutive immunity')",
                                                   "atop(g[R]~', reservoir', 'leukocyte activation rate')",
                                                   "atop(c[R]~', virus consumption', 'rate in reservoir host')",
                                                   "atop(m[R]~', reservoir', 'leukocyte mortality rate')"))#,
  #"T[vS]~', spillover host tolerance of direct virus pathology'",
  #"T[wS]~', spillover host tolerance of immunopathology'"))
  
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
    ylab(bquote(r[R]^"*"~", virus growth")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_text(size=15),
          axis.title.y = element_text(size=18), strip.background = element_rect(fill="white"), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(.5,.1,0,1.3), "cm"),
          legend.position = c(.1,.8), legend.title = element_blank(), legend.text = element_text(size=10)) + 
    scale_color_manual(values=colz) + guides(color=guide_legend(nrow=2)) + coord_cartesian(ylim = c(0,5))  
  #  print(pA)
  
  
  pB <- ggplot(data=beta.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(beta[r[R]]^"*"~", transmission")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.text = element_blank(),
          axis.title.y = element_text(size=18), strip.background = element_blank(), axis.text.y = element_text(size=14),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), plot.margin = unit(c(0,.1,0,.4), "cm")) + coord_cartesian(ylim = c(0,.005)) +
    scale_color_manual(values=colz)
  # print(pB)
  
  
  
  pC <- ggplot(data=a.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(~label, scales="free", labeller = label_parsed) +
    ylab(bquote(alpha[r[R]]^"*"~", virulence")) + theme_bw() +
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
  plot1 <- par.plots(dat.var = dat.var,
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
         height=110, 
         scale=3, 
         dpi=300)
}

plot.join.complete(dat.var = subset(out.curves, variable !="Tv_spill" & variable!="Tw_spill"),
          tol.shape= "complete-tolerance",
          dat.tol=tol.heatmap,
          sub=TRUE,
          filename=paste0(homewd,"supp-figs/FigS2.png"))


plot.join.complete(dat.var = subset(out.curves, variable !="Tv_spill" & variable!="Tw_spill"),
                   tol.shape= "complete-tolerance",
                   dat.tol=tol.heatmap,
                   sub=TRUE,
                   filename=paste0(homewd,"supp-figs/FigS2.eps"))



Fig.S3 <- function(dat.var,  filename){
  
  # and add labels for strip labels
  dat.var$label= NA
  
  # dat.var$label[dat.var$variable=="c"] <- "c~', virus consumption'"
  # dat.var$label[dat.var$variable=="g"] <- "g~', leukocyte activation'"
  # dat.var$label[dat.var$variable=="m"] <- "m~', leukocyte mortality'"
  # dat.var$label[dat.var$variable=="mu"] <- "mu~', host background mortality'"
  # dat.var$label[dat.var$variable=="g0"] <- "g[0]~', constitutive immunity'"
  dat.var$label[dat.var$variable=="Tv_spill"] <- "T[vS]~', spillover host tolerance of direct virus pathology'"
  dat.var$label[dat.var$variable=="Tw_spill"] <- "T[wS]~', spillover host tolerance of immunopathology'"
  
  
  dat.var$label = factor(dat.var$label, levels = c(#"mu~', host background mortality'",
                                                   #"g[0]~', constitutive immunity'",
                                                   #"g~', leukocyte activation'",
                                                   #"c~', virus consumption'",
                                                   #"m~', leukocyte mortality'",
                                                   "T[vS]~', spillover host tolerance of direct virus pathology'",
                                                   "T[wS]~', spillover host tolerance of immunopathology'"))
  
  #dat.var = subset(dat.var, tolerance_shape == tol.shape)
  
    dat.var$tolerance = factor(dat.var$tolerance, levels=c("Tw=1.9", "Tv=1.9", "Tw=1.5", "Tv=1.5","Tw=1.1", "Tv=1.1",
                                                           "Tw=0.9", "Tv=0.9", "Tw=0.5", "Tv=0.5","Tw=0.2", "Tv=0.2"))
  
    
  
  #dat.var = subset(dat.var, variable_par <=1)
  dat.var$value[dat.var$value<0]<-0
  
  #dat.var = subset(dat.var, tolerance_shape==tol.shape)
  
  # just rstar for now
  # r.dat = subset(dat.var, outcome=="rstar" & tolerance_shape==tol.shape)
  # a.dat = subset(dat.var, outcome=="alphastar" & tolerance_shape==tol.shape)
  # #p.dat = subset(dat.var, outcome=="prevalence" & tolerance_shape==tol.shape)
  # #inf.dat = subset(dat.var, outcome=="absInfMort" & tolerance_shape==tol.shape)
  # beta.dat = subset(dat.var, outcome=="betastar" & tolerance_shape==tol.shape)
   spill.dat = subset(dat.var, outcome=="alphaSpillover" )
  
  
  
    
    colz = c("Tw=1.9" = "dodgerblue4", 
             "Tv=1.9" = "firebrick4",
             "Tw=1.5"= "dodgerblue3", 
             "Tv=1.5" = "firebrick3",  
             "Tw=1.1" = "dodgerblue1",
             "Tv=1.1" ="firebrick1", 
             "Tw=0.9"  = "dodgerblue4",
             "Tv=0.9"=  "firebrick4",
             "Tw=0.5"= "dodgerblue3",
             "Tv=0.5"=  "firebrick3",
             "Tw=0.2"= "dodgerblue1", 
             "Tv=0.2"=  "firebrick1")
    
  
  
  
  
  
  
  pSpilloverTolerance <- ggplot(data=spill.dat) +
    geom_line(aes(x=variable_par, y=value, color=tolerance), show.legend = F) +
    facet_grid(tolerance_shape~label, scales="free", labeller = label_parsed) +
    ylab(bquote(alpha[S]~", spillover virulence")) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
          strip.text = element_text(size=14),
          axis.title.y = element_text(size=18), 
          strip.background = element_blank(), axis.text = element_text(size=14),
          plot.margin = unit(c(0,.2,1,1), "cm")) + coord_cartesian(ylim = c(0,100), xlim=c(0,2)) +
    scale_color_manual(values=colz) + scale_x_continuous(labels = label_comma())
  print(pSpilloverTolerance)
  
  ggsave(file =filename,
         plot = pSpilloverTolerance,
         units="mm",  
         width=90, 
         height=70, 
         scale=3, 
         dpi=300)
  
  
  
}

out.spill = subset(out.curves, variable =="Tv_spill" | variable=="Tw_spill")
out.spill$tolerance_shape <- factor(out.spill$tolerance_shape, levels = c("constant-tolerance", "complete-tolerance"))


Fig.S3(dat.var = out.spill,
       filename= paste0(homewd,"supp-figs/FigS3.png"))


Fig.S3(dat.var = out.spill,
       filename= paste0(homewd,"supp-figs/FigS3.eps"))

       