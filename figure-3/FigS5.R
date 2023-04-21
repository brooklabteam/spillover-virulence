
# This script produces Figure 3 of the main text.

rm(list=ls())

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(ggimage)
library(sjPlot)
library(lme4)

homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))




predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)

#and sum to add with order.date

# Plot rstar complete tolerance

#Complete tolerance
predict.dat <- arrange(predict.dat, desc(rstar_complete))
predict.dat$order = factor(predict.dat$order, levels=rev(unique(predict.dat$order)))

# And add images by order

order.dat <- read.csv(file=paste0(homewd,"/phylo-tree/Timetree_ReservoirMapping.csv"), header = T, stringsAsFactors = F)
head(order.dat)

order.dat$order

# Rename
order.dat$species
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species[order.dat$species=="Homo_sapiens"] <-  "Macaca_mulatta"
order.dat$species[order.dat$species=="Rattus_rattus"] <-  "Mus_musculus_domesticus"
order.dat$species[order.dat$species=="Antilocapra_americana"] <-  "Sus_scrofa"
order.dat$species[order.dat$species=="Phascolarctos_cinereus"] <-  "Macropus_rufus"
order.dat$species[order.dat$species=="Caenolestes_sangay"]<- "Caenolestes_convelatus"
order.dat$species[order.dat$species=="Sarcophilus_harrisii"] <- "Dasyurus_viverrinus"
order.dat$species[order.dat$species=="Oryctolagus_cuniculus"] <- "Ochotona_princeps"
order.dat$species[order.dat$species=="Tupaia_glis"] <- "Dermoptera"

# Load images - will take a moment
pic.df <- ggimage::phylopic_uid(order.dat$species) 

pic.df$order <- order.dat$order

unique(predict.dat$order)
unique(pic.df$order)

# And select those orders from which we can make a prediction:

pic.df = subset(pic.df, order=="Afrosoricida" | order == "Carnivora" | order=="Cetartiodactyla" | order=="Chiroptera" |order== "Cingulata" | 
                  order == "Dasyuromorphia" | order == "Didelphimorphia" | order=="Diprotodontia"  | order=="Eulipotyphla" | order=="Hyracoidea" |
                  order=="Monotremata" | order == "Peramelemorphia" | order == "Perissodactyla" | order=="Pilosa" | order=="Primates"  |
                  order == "Proboscidea" | order=="Rodentia" | order == "Scandentia" | order == "Tubulidentata" )

plot.df = subset(predict.dat, !is.na(g0))
setdiff(pic.df$order, plot.df$order)
setdiff(plot.df$order, pic.df$order)


#load color bar
load(paste0(homewd, subwd, "/color-bar.Rdata"))

pA <- ggplot(data=plot.df) + 
  geom_errorbar(aes(x=order, ymin=rstar_complete_lci, ymax=rstar_complete_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(x=order, y=rstar_complete, fill=order, size=N_cumulative), pch=23) +
  scale_color_manual(values=colz, guide="none")+  coord_flip(ylim=c(-0.05,1.1)) +
  scale_fill_manual(values=colz, guide="none")+theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 13, angle = 90),
        axis.title.x = element_text(size=18),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.title = element_blank(), legend.direction = "horizontal",
        legend.position = c(.82,.05),
        legend.background = element_rect(fill='transparent'),
        plot.margin = unit(c(1.8,2,.5,.5), "lines")) +
  ylab(bquote(atop("",r~"*, optimal virus growth rate in reservoir (complete tolerance)"))) +
  geom_phylopic(data=pic.df, aes(x=order, y = -.045, image=uid, color=order), size=.05)
#ylab(bquote(atop("optimal virus growth rate", "in reservoir,"~r~"* (constant tolerance)")))


#and R star

predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)

#

load(paste0(homewd,"source/gam.dat.Guth.et.al.2021.Rdata"))

#rank by descending alpha, then confidence
gam.dat <- arrange(gam.dat, desc(alpha), desc(Nobs))


# Merge our nested model predictions with those from Guth et al. 2022

# Make columns match
head(gam.dat)
names(gam.dat)[1] <- "order"
names(gam.dat)[3] <- "N"
names(predict.dat)[names(predict.dat)=="N_cumulative"] <- "N"
head(predict.dat)

# And rescale alpha in both vectors

# Constant
predict.dat$alpha_star_human_constant[!is.na(predict.dat$alpha_star_human_constant)] <- scales::rescale(x =predict.dat$alpha_star_human_constant[!is.na(predict.dat$alpha_star_human_constant)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_constant_lci[!is.na(predict.dat$alpha_star_human_constant_lci)] <- scales::rescale(x =predict.dat$alpha_star_human_constant_lci[!is.na(predict.dat$alpha_star_human_constant_lci)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_constant_uci[!is.na(predict.dat$alpha_star_human_constant_uci)] <- scales::rescale(x =predict.dat$alpha_star_human_constant_uci[!is.na(predict.dat$alpha_star_human_constant_uci)], from=c(min(predict.dat$alpha_star_human_constant_lci, na.rm = T), max(predict.dat$alpha_star_human_constant_uci, na.rm = T)), to =c(0,1)) 


# Complete
predict.dat$alpha_star_human_complete[!is.na(predict.dat$alpha_star_human_complete)] <- scales::rescale(x =predict.dat$alpha_star_human_complete[!is.na(predict.dat$alpha_star_human_complete)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_complete_lci[!is.na(predict.dat$alpha_star_human_complete_lci)] <- scales::rescale(x =predict.dat$alpha_star_human_complete_lci[!is.na(predict.dat$alpha_star_human_complete_lci)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 
predict.dat$alpha_star_human_complete_uci[!is.na(predict.dat$alpha_star_human_complete_uci)] <- scales::rescale(x =predict.dat$alpha_star_human_complete_uci[!is.na(predict.dat$alpha_star_human_complete_uci)], from=c(min(predict.dat$alpha_star_human_complete_lci, na.rm = T), max(predict.dat$alpha_star_human_complete_uci, na.rm = T)), to =c(0,1)) 


# And do the same for gam.dat
gam.dat$alpha <- scales::rescale(x =gam.dat$alpha, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
gam.dat$alpha_lci <- scales::rescale(x =gam.dat$alpha_lci, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
gam.dat$alpha_uci <- scales::rescale(x =gam.dat$alpha_uci, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 


# And merge the data
gam.plot.dat <-  dplyr::select(gam.dat, order, N, alpha, alpha_lci, alpha_uci)
share.dat.constant <- dplyr::select(predict.dat,order, N, alpha_star_human_constant, alpha_star_human_constant_lci, alpha_star_human_constant_uci)
share.dat.complete <- dplyr::select(predict.dat,order, N, alpha_star_human_complete, alpha_star_human_complete_lci, alpha_star_human_complete_uci)
names(share.dat.complete) <- names(share.dat.constant) <- names(gam.plot.dat)

share.dat.constant <- arrange(share.dat.constant, desc(alpha))
share.dat.complete <- arrange(share.dat.complete, desc(alpha))


#gam.plot.dat[,3:5] <- -1*gam.plot.dat[,3:5] 
share.dat.constant[,3:5] <- -1*share.dat.constant[,3:5] 
share.dat.complete[,3:5] <- -1*share.dat.complete[,3:5] 
share.dat.complete$tolerance = "complete"
share.dat.constant$tolerance = "constant"
gam.plot.dat$tolerance = "natural"
share.dat.complete$source <- share.dat.constant$source <-  "predicted from\nnested model"
gam.plot.dat$source <- "predicted from zoonoses"
plot.dat <- rbind(gam.plot.dat, share.dat.complete, share.dat.constant)
plot.dat$source <- factor(plot.dat$source, levels=c("predicted from zoonoses", "predicted from\nnested model"))

# Reorder, ranked by the complete data
share.dat.complete <- arrange(share.dat.complete, desc(alpha))
plot.dat$order <- factor(plot.dat$order, levels=unique(share.dat.complete$order))

head(plot.dat)
plot.dat


# And take only the complete data
# first edit
plot.dat <- plot.dat[complete.cases(plot.dat),]

#shapez <- c("predicted from zoonoses"=25, "predicted from\nnested model"=24)

# And add images by order

order.dat <- read.csv(file=paste0(homewd,"/phylo-tree/Timetree_ReservoirMapping.csv"), header = T, stringsAsFactors = F)
head(order.dat)

order.dat$order

# Rename
order.dat$species
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species[order.dat$species=="Homo_sapiens"] <-  "Macaca_mulatta"
order.dat$species[order.dat$species=="Rattus_rattus"] <-  "Mus_musculus_domesticus"
order.dat$species[order.dat$species=="Antilocapra_americana"] <-  "Sus_scrofa"
order.dat$species[order.dat$species=="Phascolarctos_cinereus"] <-  "Macropus_rufus"
order.dat$species[order.dat$species=="Caenolestes_sangay"]<- "Caenolestes_convelatus"
order.dat$species[order.dat$species=="Sarcophilus_harrisii"] <- "Dasyurus_viverrinus"
order.dat$species[order.dat$species=="Oryctolagus_cuniculus"] <- "Ochotona_princeps"
order.dat$species[order.dat$species=="Tupaia_glis"] <- "Dermoptera"

# Load images - will take a moment
pic.df <- ggimage::phylopic_uid(order.dat$species) 

pic.df$order <- order.dat$order

unique(plot.dat$order)
unique(pic.df$order)

# And select those orders from which we can make a prediction:

pic.df = subset(pic.df, order=="Afrosoricida" | order == "Carnivora" | order=="Cetartiodactyla" | order=="Chiroptera" |order== "Cingulata" | 
                  order == "Dasyuromorphia" | order == "Didelphimorphia" | order=="Diprotodontia"  | order=="Eulipotyphla" | order=="Hyracoidea" |
                  order=="Monotremata" | order == "Peramelemorphia" | order == "Perissodactyla" | order=="Pilosa" | order=="Primates"  |
                  order == "Proboscidea" | order=="Rodentia" | order == "Scandentia" | order == "Tubulidentata" )

setdiff(pic.df$order, plot.dat$order)
setdiff(plot.dat$order, pic.df$order)

pBa <- ggplot(data=subset(plot.dat, tolerance=="complete"))  +  geom_hline(aes(yintercept=0), size=.2) +
  geom_errorbar(aes(x=order, ymin=alpha_lci, ymax=alpha_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(order, alpha, fill=order, size=N), shape=24) + 
  theme_bw() +
  scale_color_manual(values=colz, guide="none") +
  scale_fill_manual(values=colz, guide="none") +
  #scale_shape_manual(values=shapez, guide="none") +
  scale_x_discrete(position = "top") +
  #facet_grid(source~., scales = "free_y") +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size=18), 
        legend.direction = "horizontal", legend.position = c(.15,.05), legend.title = element_blank(),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        plot.margin = unit(c(1,.5,.5,3), "cm")) + 
  ylab(bquote(atop("","relative spillover virulence,"~alpha[S]~" (complete tolerance)"))) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_flip(ylim=c(-1.1,1.1), clip = "off") +
  geom_phylopic(data=pic.df, aes(x=order, y = 1.095, image=uid, color=order), size=.05)

pB <- pBa + geom_text(x=20.5, y=0, label="From nested model                    From zoonotic literature", angle=0, nudge_y = 2, size=5.5, fontface="bold") + 
  coord_flip(ylim=c(-1.1,1.1), clip = "off") 


### And combine

pFigS5 <- cowplot::plot_grid(pA,pB, ncol = 2, nrow = 1, labels = c("A", "B"), label_size = 22, rel_widths = c(.8,1), label_x = c(0,.07))



ggsave(file = paste0(homewd,"/supp-figs/FigS5.png"),
       plot = pFigS5,
       bg="white",
       units="mm",  
       width=160, 
       height=50, 
       scale=3, 
       dpi=300)

