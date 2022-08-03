# This is the script to produce Fig S3. 
# Much of it replicates the brief_script.R


rm(list=ls())

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(RColorBrewer)

homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))


######################################################
######################################################

## load data

predict.dat <- read.csv(paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)

head(predict.dat)

## Now plot mu, Tw (constant/complete), g0, and Tv_human (Tvs) (constant/complete)
## across all the mammalian orders. Do not include order axis on the first three
## so that the plots can be combined


colz = scales::hue_pal()(length(unique(predict.dat$order))-1)
colz=c(colz, "red")

names(colz) <- c(unique(predict.dat$order)[unique(predict.dat$order)!="Chiroptera"], "Chiroptera")

y.int = 0.08894968 #from the brief_script.R file

## First, mu
p1 <- ggplot(data=predict.dat) + 
      geom_point(aes(x=order, y=mu, fill=order, size=N_mu), pch=21) + 
      theme_bw()  + 
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(.5,2.7,0,.5), "lines"), 
            legend.position = c(.8,.91), legend.direction = "horizontal", 
            legend.title = element_blank(),
            panel.grid = element_blank()) +
      scale_fill_manual(values=colz, guide="none") +
      scale_color_manual(values=colz) +
      geom_errorbar(aes(x=order, ymin=mu_lci, ymax=mu_uci, color=order), width=0, show.legend = F) +
      geom_hline(aes(yintercept=y.int), linetype=2) +
      ylab(bquote("predicted annual mortality,"~mu~"("~yrs^-1~")")) 

p1 #will combine and edit below for Fig S3

## Now, Tw (on one plot)


# Constant tolerance (Tw)
p2 <- ggplot(data=predict.dat) + 
  geom_point(aes(x=order, y=Tw_constant, fill=order, size=N_Tw), pch=21) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(.1,.5,0,.7), "lines"), 
        legend.position = c(.82,.91),
        panel.grid = element_blank(),
        legend.direction = "horizontal", legend.title = element_blank()) +
  scale_fill_manual(values=colz, guide="none") +
  scale_color_manual(values=colz) +
  geom_errorbar(aes(x=order, ymin=Tw_constant_lci, ymax=Tw_constant_uci, color=order), width=0, show.legend = F) +
  geom_hline(aes(yintercept=1.5), linetype=2) + ylab(bquote("constant immunopathology tolerance,"~T[w])) +
  scale_y_continuous(sec.axis = sec_axis(~ . -1, name = (bquote("complete immunopathology tolerance,,"~T[w]))))

p2

# Slim your colors down to only those in this dataset
#colz2 = colz[unique(predict.dat$order[!is.na(predict.dat$g0)])]

# g0
p3 <- ggplot(data=predict.dat) + 
  geom_point(aes(x=order, y=g0, fill=order, size=N_g0), pch=21) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(.1,2.8,0,.4), "lines"),
        panel.grid = element_blank(),
        legend.position = c(.75,.92),
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  geom_errorbar(aes(x=order, ymin=g0_lci, ymax=g0_uci, color=order), width=0, show.legend = F) +
  scale_color_manual(values=colz) +
  scale_fill_manual(values=colz, guide="none") +
  geom_hline(aes(yintercept=0.5), linetype=2) + ylab(bquote("magnitude constitutive immunity, "~g[0]))


#  Tvs
p4 <- ggplot(data=predict.dat) + theme_bw() +
        geom_point(aes(x=order, y=Tv_human_constant, fill=order), pch=21, show.legend = F, size=3) + 
        scale_fill_manual(values=colz) + ylab(bquote("constant viral spillover tolerance,"~T[vS])) +
        theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
        plot.margin = unit(c(0,.4,.5,.5), "lines"),
        panel.grid = element_blank()) +
        geom_hline(aes(yintercept=1.5), linetype=2) +
        scale_y_continuous(sec.axis = sec_axis(~ . -1, name = (bquote("complete viral spillover tolerance,"~T[vS]))))



# And put together

FigS3 <- cowplot::plot_grid(p1,p2,p3,p4, nrow=4, ncol=1, rel_heights = c(1.05,1,1,1.3), labels = c("A", "B", "C", "D"), label_size = 22, label_x = -.01)


#and save
ggsave(file = paste0(homewd, "supp-figs/FigS3.png"),
       plot = FigS3,
       units="mm",  
       width=60, 
       height=140, 
       scale=3, 
       dpi=300)


