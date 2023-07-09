
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

homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))


######################################################
######################################################
## Get your color bar to match

# Now make the comparisons
homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"

load(paste0(homewd, subwd, "/color-bar.Rdata"))

predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)

######################################################
######################################################



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

# And merge
setdiff(unique(predict.dat$order), unique(gam.dat$order)) # here, we lack data but can make predictions

# remove these from the comparison
gam.choice <- dplyr::select(gam.dat, order)

predict.dat = merge(predict.dat, gam.choice, by="order")


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


gam.plot.dat[,3:5] <- -1*gam.plot.dat[,3:5] 
share.dat.complete$tolerance = "complete"
share.dat.constant$tolerance = "constant"
gam.plot.dat$tolerance = "natural"
share.dat.complete$source <- share.dat.constant$source <-  "predicted from\nnested model"
gam.plot.dat$source <- "predicted from zoonoses"
plot.dat <- rbind(gam.plot.dat, share.dat.complete, share.dat.constant)
plot.dat$source <- factor(plot.dat$source, levels=c("predicted from zoonoses", "predicted from\nnested model"))

# Reorder, ranked by the constant data
plot.dat$order <- factor(plot.dat$order, levels=unique(arrange(share.dat.constant, desc(alpha))$order))

head(plot.dat)

# And take only the complete data
plot.dat <- subset(plot.dat, !is.na(alpha))#plot.dat[complete.cases(plot.dat),]

#and now compare predictions
unique(plot.dat$tolerance)
unique(plot.dat$source[plot.dat$tolerance=="natural"])
nat.dat = subset(plot.dat, tolerance=="natural")
nat.dat <- dplyr::select(nat.dat, order, N, alpha, alpha_lci, alpha_uci)
names(nat.dat) <- c("order", "lit_N", "lit_alpha", "lit_alpha_lci", "lit_alpha_uci")

dat.compare=subset(plot.dat, tolerance!="natural")
unique(dat.compare$source)
head(dat.compare)

dat.compare <- merge(dat.compare, nat.dat, by=c("order"), all.x = T)
dat.compare$lit_alpha <- dat.compare$lit_alpha*-1 

p1 <- ggplot(data = dat.compare) + theme_bw() +
      geom_point(aes(x=lit_alpha, y=alpha, color=order), size=3) +
      facet_grid(~tolerance)

m1 <- lm(alpha~lit_alpha, dat=subset(dat.compare, tolerance=="complete"))
summary(m1) #Multiple R-squared:  0.6208. Adjusted R-squared:  0.5577.  p-value: 0.02021

m2 <- lm(alpha~lit_alpha, dat=subset(dat.compare, tolerance=="constant"))
summary(m2) #Multiple R-squared:  0.5696,	Adjusted R-squared:  0.4979. p-value: 0.03044

#and add to the plot
dat.compare$predicted_alpha <- NA

dat.compare$predicted_alpha[dat.compare$tolerance=="complete" & !is.na(dat.compare$lit_alpha)] <- predict(m1)
dat.compare$predicted_alpha[dat.compare$tolerance=="constant" & !is.na(dat.compare$lit_alpha)] <- predict(m2)

droplevels(dat.compare$order)

dat.compare$residuals = dat.compare$alpha-dat.compare$predicted_alpha
dat.compare$residuals_round <- round(dat.compare$residuals, 2)
dat.compare$residuals_round[dat.compare$residuals_round ==0.00 & !is.na(dat.compare$residuals)] <- round(dat.compare$residuals[dat.compare$residuals_round ==0.00 & !is.na(dat.compare$residuals)], 3)
dat.compare$label_point_y= dat.compare$alpha - (dat.compare$residuals/2)
dat.compare$label_point_x = dat.compare$lit_alpha
dat.compare$label_point_x[dat.compare$residuals<0 & !is.na(dat.compare$residuals)] <- dat.compare$label_point_x[dat.compare$residuals<0 & !is.na(dat.compare$residuals)] +.025
dat.compare$label_point_x[dat.compare$residuals>0 & !is.na(dat.compare$residuals)] <- dat.compare$label_point_x[dat.compare$residuals>0 & !is.na(dat.compare$residuals)] -.025
dat.compare$label_point_x[dat.compare$order=="Chiroptera"] <- dat.compare$lit_alpha[dat.compare$order=="Chiroptera"] -.05
dat.compare$label_point_x[dat.compare$order=="Rodentia"] <- dat.compare$lit_alpha[dat.compare$order=="Rodentia"] -.01
dat.compare$label_point_y[dat.compare$order=="Rodentia"] <- dat.compare$predicted_alpha[dat.compare$order=="Rodentia"] +.02


trendline.dat <- cbind.data.frame(r_sq=c(summary(m1)$r.squared, summary(m2)$r.squared), tolerance=c("complete", "constant"))
trendline.dat$label= paste0("R^'2'~'='~", round(trendline.dat$r_sq, 2))

trendline.dat$tolerance <- factor(trendline.dat$tolerance,levels = c("constant", "complete"))    
dat.compare$tolerance <- factor(dat.compare$tolerance, levels=c("constant", "complete"))

FigS6 <- ggplot(data = subset(dat.compare, !is.na(lit_alpha))) + theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.title = element_blank(), legend.position = c(.87,.15), 
        legend.background = element_rect(color="black"), strip.text = element_text(size = 18),
        axis.title = element_text(size=18), axis.text = element_text(size=14)) +
  geom_linerange(aes(x=lit_alpha, ymin = predicted_alpha, ymax=alpha), linetype=2) +
  geom_line(aes(x=lit_alpha, y=predicted_alpha)) +
  ylab("relative virulence predicted from nested model") +
  xlab("relative virulence predicted from zoonoses in the literature") + 
  geom_label(aes(x=label_point_x, y=label_point_y, label = residuals_round), label.size = NA,size=4) +
  geom_label(data=trendline.dat, aes(x=.09, y=.61, label=label), parse = T, fontface="italic") +
  geom_point(aes(x=lit_alpha, y=alpha, color=order), size=3) +
  facet_grid(tolerance~.) + scale_color_manual(values=colz)

FigS6

ggsave(file =paste0(homewd,"supp-figs/FigS6.png"),
       plot = FigS6,
       units="mm",  
       width=65, 
       height=80, 
       scale=3, 
       dpi=300)

ggsave(file =paste0(homewd,"supp-figs/FigS6.eps"),
       plot = FigS6,
       units="mm",  
       width=65, 
       height=80, 
       scale=3, 
       dpi=300)

#and save this comparison data for fitting - but accompany with the parameters of interes
names(predict.dat)
predict.merge <- dplyr::select(predict.dat, order, mu, mu_lci, mu_uci, Tw_constant, Tw_constant_lci, Tw_constant_uci, Tw_complete, Tw_complete_lci, Tw_complete_uci, g0, g0_lci, g0_uci, Tv_human_constant, Tv_human_complete)

dat.merge <- merge(dat.compare, predict.merge, by="order", all.x = T)


dat.merge <- dplyr::select(dat.merge, -(source), -(residuals), -(predicted_alpha), -(residuals_round), -(label_point_x), -(label_point_y))

dat.merge$Tw_complete[dat.merge$tolerance=="constant"] <- dat.merge$Tw_constant[dat.merge$tolerance=="constant"]
dat.merge$Tw_complete_lci[dat.merge$tolerance=="constant"] <- dat.merge$Tw_constant_lci[dat.merge$tolerance=="constant"]
dat.merge$Tw_complete_uci[dat.merge$tolerance=="constant"] <- dat.merge$Tw_constant_uci[dat.merge$tolerance=="constant"]

dat.merge$Tw_constant[dat.merge$tolerance=="complete"] <- dat.merge$Tw_complete[dat.merge$tolerance=="complete"]
dat.merge$Tw_constant_lci[dat.merge$tolerance=="complete"] <- dat.merge$Tw_complete_lci[dat.merge$tolerance=="complete"]
dat.merge$Tw_constant_uci[dat.merge$tolerance=="complete"] <- dat.merge$Tw_complete_uci[dat.merge$tolerance=="complete"]

dat.merge$Tv_human_complete[dat.merge$tolerance=="constant"] <- dat.merge$Tv_human_constant[dat.merge$tolerance=="constant"] 
dat.merge$Tv_human_constant[dat.merge$tolerance=="complete"] <- dat.merge$Tv_human_complete[dat.merge$tolerance=="complete"] 

dat.merge <- dplyr::select(dat.merge, -(Tv_human_complete), -(Tw_complete), -(Tw_complete_lci), -(Tw_complete_uci))
names(dat.merge)[names(dat.merge)=="Tw_constant"] <- "Tw"
names(dat.merge)[names(dat.merge)=="Tw_constant_lci"] <- "Tw_lci"
names(dat.merge)[names(dat.merge)=="Tw_constant_uci"] <- "Tw_uci"

names(dat.merge)[names(dat.merge)=="Tv_human_constant"] <- "Tv_human"

write.csv(dat.merge, file = paste0(homewd, "figure-3/compare_predictions.csv"), row.names = F)



######################################################
######################################################

rm(list=ls())

# Now make the comparisons
homewd =  "/Users/carabrook/Developer/spillover-virulence/"
subwd = "figure-3"
setwd(paste0(homewd, subwd))


#rerun section 1 above
load(paste0(homewd, subwd, "/color-bar.Rdata"))



######################################################
######################################################
## Get your color bar to match

predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)

# and with rabies excluded
load(paste0(homewd,"source/SI.gam.dat.no.rabies.Guth.et.al.2021.Rdata"))

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


gam.plot.dat[,3:5] <- -1*gam.plot.dat[,3:5] 
share.dat.complete$tolerance = "complete"
share.dat.constant$tolerance = "constant"
gam.plot.dat$tolerance = "natural"
share.dat.complete$source <- share.dat.constant$source <-  "predicted from\nnested model"
gam.plot.dat$source <- "predicted from zoonoses"
plot.dat <- rbind(gam.plot.dat, share.dat.complete, share.dat.constant)
plot.dat$source <- factor(plot.dat$source, levels=c("predicted from zoonoses", "predicted from\nnested model"))

# Reorder, ranked by the constant data
plot.dat$order <- factor(plot.dat$order, levels=unique(arrange(share.dat.constant, desc(alpha))$order))

head(plot.dat)

# And take only the complete data
plot.dat <- subset(plot.dat, !is.na(alpha))#plot.dat <- plot.dat[complete.cases(plot.dat),]

#and now compare predictions
unique(plot.dat$tolerance)
unique(plot.dat$source[plot.dat$tolerance=="natural"])
nat.dat = subset(plot.dat, tolerance=="natural")
nat.dat <- dplyr::select(nat.dat, order, N, alpha, alpha_lci, alpha_uci)
names(nat.dat) <- c("order", "lit_N", "lit_alpha", "lit_alpha_lci", "lit_alpha_uci")

dat.compare=subset(plot.dat, tolerance!="natural")
unique(dat.compare$source)
head(dat.compare)

dat.compare <- merge(dat.compare, nat.dat, by=c("order"), all.x = T)
dat.compare$lit_alpha <- dat.compare$lit_alpha*-1 

p1 <- ggplot(data = dat.compare) + theme_bw() +
  geom_point(aes(x=lit_alpha, y=alpha, color=order), size=3) +
  facet_grid(~tolerance)

m1 <- lm(alpha~lit_alpha, dat=subset(dat.compare, tolerance=="complete" & !is.na(lit_N)))
summary(m1) #Multiple R-squared:  0.5918,	Adjusted R-squared:  0.5238. p-value: 0.02563

m2 <- lm(alpha~lit_alpha, dat=subset(dat.compare, tolerance=="constant" & !is.na(lit_N)))
summary(m2) #Multiple R-squared:  0.5564,	Adjusted R-squared:  0.4825 . p-value: 0.03359

#and add to the plot
dat.compare$predicted_alpha <- NA

dat.compare$predicted_alpha[dat.compare$tolerance=="complete" & !is.na(dat.compare$lit_alpha)] <- predict(m1)
dat.compare$predicted_alpha[dat.compare$tolerance=="constant" & !is.na(dat.compare$lit_alpha)] <- predict(m2)

droplevels(dat.compare$order)

dat.compare$residuals = dat.compare$alpha-dat.compare$predicted_alpha
dat.compare$residuals_round <- round(dat.compare$residuals, 2)
dat.compare$residuals_round[dat.compare$residuals_round ==0.00 & !is.na(dat.compare$residuals)] <- round(dat.compare$residuals[dat.compare$residuals_round ==0.00 & !is.na(dat.compare$residuals)], 3)
dat.compare$label_point_y= dat.compare$alpha - (dat.compare$residuals/2)
dat.compare$label_point_x = dat.compare$lit_alpha
dat.compare$label_point_x[dat.compare$residuals<0 & !is.na(dat.compare$residuals)] <- dat.compare$label_point_x[dat.compare$residuals<0 & !is.na(dat.compare$residuals)] +.025
dat.compare$label_point_x[dat.compare$residuals>0 & !is.na(dat.compare$residuals)] <- dat.compare$label_point_x[dat.compare$residuals>0 & !is.na(dat.compare$residuals)] -.025
dat.compare$label_point_x[dat.compare$order=="Chiroptera"] <- dat.compare$lit_alpha[dat.compare$order=="Chiroptera"] -.05
dat.compare$label_point_y[dat.compare$order=="Chiroptera"] <- dat.compare$label_point_y[dat.compare$order=="Chiroptera"] + .01
dat.compare$label_point_x[dat.compare$order=="Rodentia"] <- dat.compare$lit_alpha[dat.compare$order=="Rodentia"] -.01
dat.compare$label_point_y[dat.compare$order=="Rodentia"] <- dat.compare$predicted_alpha[dat.compare$order=="Rodentia"] +.01
dat.compare$label_point_y[dat.compare$order=="Rodentia"& dat.compare$tolerance!="complete"] <- dat.compare$predicted_alpha[dat.compare$order=="Rodentia" & dat.compare$tolerance!="complete"] -.025

dat.compare$label_point_x[dat.compare$order=="Cetartiodactyla"] <- dat.compare$lit_alpha[dat.compare$order=="Cetartiodactyla"] -.03
dat.compare$label_point_y[dat.compare$order=="Cetartiodactyla"] <- dat.compare$label_point_y[dat.compare$order=="Cetartiodactyla"] -.01
dat.compare$label_point_y[dat.compare$order=="Carnivora" & dat.compare$tolerance=="complete"] <- dat.compare$label_point_y[dat.compare$order=="Carnivora" & dat.compare$tolerance=="complete"] +.01
dat.compare$label_point_y[dat.compare$order=="Carnivora" & dat.compare$tolerance!="complete"] <- dat.compare$label_point_y[dat.compare$order=="Carnivora" & dat.compare$tolerance!="complete"] + .01
dat.compare$label_point_x[dat.compare$order=="Carnivora" & dat.compare$tolerance!="complete"] <- dat.compare$lit_alpha[dat.compare$order=="Carnivora" & dat.compare$tolerance!="complete"] +.02
dat.compare$label_point_y[dat.compare$order=="Primates"] <- dat.compare$label_point_y[dat.compare$order=="Primates"] +.01

trendline.dat <- cbind.data.frame(r_sq=c(summary(m1)$r.squared, summary(m2)$r.squared), tolerance=c("complete", "constant"))
trendline.dat$label= paste0("R^'2'~'='~", round(trendline.dat$r_sq, 2))


trendline.dat$tolerance <- factor(trendline.dat$tolerance,levels = c("constant", "complete"))    
dat.compare$tolerance <- factor(dat.compare$tolerance, levels=c("constant", "complete"))


FigS8 <- ggplot(data = subset(dat.compare, !is.na(lit_alpha))) + theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        legend.title = element_blank(), legend.position = c(.87,.15), 
        legend.background = element_rect(color="black"), strip.text = element_text(size = 18),
        axis.title = element_text(size=18), axis.text = element_text(size=14)) +
  geom_line(aes(x=lit_alpha, y=predicted_alpha)) +
  ylab("relative virulence predicted from nested model") +
  xlab("relative virulence predicted from zoonoses in the literature") + 
  geom_label(aes(x=label_point_x, y=label_point_y, label = residuals_round), label.size = NA,size=4) +
  geom_label(data=trendline.dat, aes(x=.05, y=.41, label=label), parse = T, fontface="italic") +
  geom_linerange(aes(x=lit_alpha, ymin = predicted_alpha, ymax=alpha), linetype=2) +
  geom_point(aes(x=lit_alpha, y=alpha, color=order), size=3) +
  facet_grid(tolerance~.) + scale_color_manual(values=colz)


FigS8

ggsave(file =paste0(homewd,"supp-figs/FigS8.png"),
       plot = FigS8,
       units="mm",  
       width=65, 
       height=80, 
       scale=3, 
       dpi=300)

ggsave(file =paste0(homewd,"supp-figs/FigS8.eps"),
       plot = FigS8,
       units="mm",  
       width=65, 
       height=80, 
       scale=3, 
       dpi=300)







