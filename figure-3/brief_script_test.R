# This is the script to accompany the brief


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
## 1. estimate mu

pan.dat <- read.csv(file="PanTHERIA.csv", header = T, stringsAsFactors = F)
head(pan.dat)

# And select only those columns of interest 
pan.dat <- dplyr::select(pan.dat, 1:5, X5.1_AdultBodyMass_g, X17.1_MaxLongevity_m, X10.1_PopulationGrpSize, X21.1_PopulationDensity_n.km2, X22.1_HomeRange_km2, X18.1_BasalMetRate_mLO2hr)#, X5.2_BasalMetRateMass_g)
head(pan.dat)


# Rename
names(pan.dat) <- c("order", "family", "genus", "species", "binomial", "mass_g", "max_lifespan_months", "pop_group_size", "pop_density_N_km2", "homerange_km2", "BMR_mLO2hr")#, "BMR_mg")
pan.dat$max_lifespan_yrs <- pan.dat$max_lifespan_months/12

# Add dash for binomial nomemclature
pan.dat$binomial <- sub(" ", "_", pan.dat$binomial)


# Correct errors
pan.dat$max_lifespan_yrs[pan.dat$binomial=="Galeopterus_variegates"] <- 17.5 #from animaldiversity.org

# Remove humans from dataset
pan.dat = subset(pan.dat, binomial!="Homo_sapiens")

# Rename some incorrect orders with contemporary taxonomy
pan.dat$order[pan.dat$order=="Artiodactyla"] <- "Cetartiodactyla"
pan.dat$order[pan.dat$order=="Cetacea"] <- "Cetartiodactyla"
pan.dat$order[pan.dat$order=="Soricomorpha"] <- "Eulipotyphla"
pan.dat$order[pan.dat$order=="Erinaceomorpha"] <- "Eulipotyphla"



# Convert BMR to W/G
pan.dat$BMR_mLO2hr <- (pan.dat$BMR_mLO2hr/60/60*20.1)/pan.dat$mass_g
names(pan.dat)[names(pan.dat)=="BMR_mLO2hr"] <- "BMR_W_g"

# Correct errors
pan.dat$BMR_W_g[pan.dat$binomial=="Acrobates_pygmaeus"] <- 0.084/pan.dat$mass_g[pan.dat$binomial=="Acrobates_pygmaeus"]  # from AnAge database

#add healy data
healy.dat <- read.csv(file = "Healy.csv", header = T, stringsAsFactors = F)
head(healy.dat)
healy.dat = subset(healy.dat, class=="Mammalia")
healy.dat <- dplyr::select(healy.dat, species, maximum_lifespan_yr, mass_g, BMR)

names(healy.dat) <- c("binomial", "healy_max_lifespan_yrs", "healy_mass_g",  "healy_BMR_W_g")

setdiff(healy.dat$binomial, pan.dat$binomial)
#setdiff(pan.dat$binomial, healy.dat$binomial)

#and merge the two
all.dat <- merge(pan.dat, healy.dat, all = T, by = c("binomial"))
head(all.dat)

#and average among mass and lifespan
pan.dat <- ddply(all.dat, .(order, family, genus, species, binomial), summarise, mass_g = mean(c(mass_g, healy_mass_g), na.rm=T), max_lifespan_yrs=mean(c(max_lifespan_yrs, healy_max_lifespan_yrs), na.rm=T), BMR_W_g=mean(c(BMR_W_g, healy_BMR_W_g), na.rm=T),pop_group_size = unique(pop_group_size),  pop_density_N_km2= unique(pop_density_N_km2), homerange_km2 = unique(homerange_km2))

head(pan.dat)

# And try the raw BMR
pan.dat$BMR_W <- pan.dat$BMR_W_g*pan.dat$mass_g

# Remove any data for which you lack mass:
pan.dat = subset(pan.dat, !is.na(mass_g))

# Remove any data for which you lack lifespan:
pan.dat = subset(pan.dat, !is.na(max_lifespan_yrs))

# How many have lifespan too?
length(pan.dat$max_lifespan_yrs[!is.na(pan.dat$max_lifespan_yrs)]) #1055

# And how many have BMR?
length(pan.dat$BMR_W_g[!is.na(pan.dat$BMR_W_g)]) #629

# Eventually will have wbc data too, so load it here to
# make your color bar

colz = scales::hue_pal()(length(sort(unique(pan.dat$order)))-1)
colz=c(colz, "red")

names(colz) <- c(sort(unique(pan.dat$order))[sort(unique(pan.dat$order))!="Chiroptera"], "Chiroptera")


p1 <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,   fill=order), size =3, pch=21) +  
  scale_fill_manual(values=colz) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  theme(panel.grid = element_blank()) +
  xlab("mass (g)") + ylab("max lifespan (yrs)")
p1

# Save plot


# Now, run a GAM to determine the partial effect of host order
# (phylogeny) on maximum longevity, and predict lifespan 
# rates by order:

#library(mgcv)
pan.dat$log10mass_g <- log10(pan.dat$mass_g)
pan.dat$order <- as.factor(pan.dat$order)
pan.dat$log10_max_lifespan_yrs <- log10(pan.dat$max_lifespan_yrs)

library(lme4)
library(sjPlot)
m1 <- lmer(log10_max_lifespan_yrs~ log10mass_g + (1|order), data = pan.dat) #mass longevity linear interaction with separate intercepts per order
plot_model(m1, type = "re")
summary(m1) 


# And, plot the partial effects of order:
# Source partial effects script from Mollentze and Streicker 2020:
# This script includes two plotting functions, added by me!
# 
# source(paste0(homewd,"source/mollentze-streicker-2020-functions.R"))
# 
# 
# order.dat <- get_partial_effects(m1, var="order")
# mass.dat <- get_partial_effects_continuous(m1, var="log10mass_g")
# 
# p2a <- plot.partial(df=order.dat, var="order", response_var = "max lifespan (yrs)")
# p2b <- plot.partial.cont(df=mass.dat, var="log10mass_g",  response_var = "max lifespan (yrs)", alt_var = "mass (g)", log = T)
# 
# p2 <- cowplot::plot_grid(p2a, p2b, ncol=1, nrow = 2, labels=c("A", "B"), label_x = .1)


newdata.df <- cbind.data.frame(log10mass_g=3.229679, order="Chiroptera")


# Significant positive associations with lifespan: 
# Chiroptera, Cingulata, Monotremata, Primates

# Significant negative associations with lifespan: 
# Afrosoricida, Dasyuromorphia, Didelphimorphia,
# Eulipotyphla, Notoryctemorphia, Peramelemorphia 


# And plot with predictions:

pan.dat.sum <- cbind.data.frame(log10mass_g=seq(log10(min(pan.dat$mass_g)), max(log10(pan.dat$mass_g)), length=10), order = "Chiroptera")
#pan.dat.sum$predict_lifespan <- 10^(predict.gam(m1, newdata = pan.dat.sum, exclude = "s(order)"))
pan.dat.sum$predict_lifespan <- 10^(predict(m1, newdata = pan.dat.sum, re.form=NA)) #predicts without random effects
pan.dat.sum$mass_g <- 10^(pan.dat.sum$log10mass_g)
pan.dat.sum = pan.dat.sum[pan.dat.sum$mass_g==min(pan.dat.sum$mass_g, na.rm=T) | pan.dat.sum$mass_g==max(pan.dat.sum$mass_g, na.rm=T),]

p1b <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,   fill=order), size =3, pch=21) +  
  scale_fill_manual(values=colz) +
  geom_line(data = pan.dat.sum, aes(x=mass_g, y=predict_lifespan), size=1) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  theme(panel.grid = element_blank()) +
  xlab("mass (g)") + ylab("max lifespan (yrs)")
p1b



# Now, to parameterize mu, the annual mortality rate by order for
# our within-host model, we can refit a glm simply using order as a direct
# predictor of lifespan excluding the effects of body mass, then take the inverse (but be sure to
# express in timesteps of days to make comparable to the viral dynamics):


# Remake glm with no involvement of mass
m1b <- lmer(log10_max_lifespan_yrs~ (1|order), data = pan.dat) 
summary(m1b)

plot_model(m1b, type="re")
# 
# 
# m1b <- gam(log10_max_lifespan_yrs~ s(order, bs="re"), data = pan.dat) 
# summary(m1b)


predict.dat <- cbind.data.frame(order = unique(pan.dat$order))#, log10mass_g=1)

#predict.dat <- cbind.data.frame(order = order.dat[[1]]$order, log10mass_g=1)
# You can insert any value you like for "log10mass_g" since we won't 
# be using it anyhow

predict.dat$mu <- 1/((10^predict(m1b, newdata = predict.dat, type="response"))*365)

# There is no option for computing standard errors of predictions for lmer
# because it is difficult to define an efficient method that incorporates uncertainty
# in the variance parameters; R help recommends bootMer for this task.

# We do that here:
predFun <- function(fit) {
  1/((10^predict(fit,predict.dat))*365)
}
bb <- bootMer(m1b,nsim=200,FUN=predFun,seed=101)

#now get CIs from bb
sim.out <- bb$t #gives a matrix of 200 sims of all the data

#get standard error:
predict.dat$mu_se <- apply(sim.out, 2, sd)/sqrt(200)   
predict.dat$mu_lci <- predict.dat$mu - 1.96*predict.dat$mu_se
predict.dat$mu_uci <- predict.dat$mu + 1.96*predict.dat$mu_se



predict.dat$mu_lci[predict.dat$mu_lci<0] <- 0
predict.dat$mu_uci[predict.dat$mu_uci>1] <- 1

# Also compute the number of data entries per order used to determine this:
mu.sum <- ddply(pan.dat, .(order), summarise, N_mu = length(binomial))

predict.dat <- merge(predict.dat, mu.sum, by="order", all.x = T)
predict.dat <- dplyr::select(predict.dat, -(mu_se))


# and plot your predictions by order for mass

# first, get your null
predict.dat$mu[which(predict.dat$mu==median(predict.dat$mu))]
y.int = median(predict.dat$mu)




# Now plot predictions against null:

p3 <- ggplot(data=predict.dat) + 
  theme_bw()  + 
  geom_hline(aes(yintercept=y.int), linetype=2) +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
        legend.position = c(.7,.91), legend.direction = "horizontal", 
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values=colz, guide="none") +
  scale_color_manual(values=colz) +
  geom_errorbar(aes(x=order, ymin=mu_lci, ymax=mu_uci, color=order), width=0, show.legend = F) +
  geom_point(aes(x=order, y=mu, fill=order, size=N_mu), pch=21) + 
  ylab(bquote("predicted annual mortality rate by order,"~mu~"("~days^-1~")")) 

p3

# And save

######################################################
######################################################
## 2. estimate Tw

# We can extract tolerance as the order-level effect on the 
# above interaction (rather than the prediction)

# Look at p1 again

p1 <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,   fill=order), size =3, pch=21) +  
  scale_fill_manual(values=colz) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  xlab("mass (g)") + ylab("max lifespan (yrs)")
p1


# m1 from above should already be appropriate,
# and Example Fig. 2 is also appropriate

# Make tolerance (Tw) a scaled version of the order effect in this relationship
# First, get high and low confidence and estimate from this plot and scale as Tw

tmp.dat = cbind.data.frame(order= plot_model(m1, type="re")$data$term, estimate=plot_model(m1, type="re")$data$estimate, conf.low=plot_model(m1, type="re")$data$conf.low, conf.high=plot_model(m1, type="re")$data$conf.high) 


#tmp.dat <- cbind.data.frame(, estimate=order.dat[[1]]$y, conf.low=order.dat[[1]]$ylower, conf.high=order.dat[[1]]$yupper)

# Our goal is for Tw to span 0 to 1 for complete tolerance assumptions
# and be >1 for constant tolerance assumptions. We only allow for 
# linear transformations of the data, in order to retain differences
# in the magnitude of effect. We add to the predict.dat database from above
# among the orders. 


# Linear transformation: 
# Constant tolerance: Make all effects positive, then add 1
# Complete tolerance: Make all effects positive
tmp.dat$Tw_constant <- tmp.dat$estimate + abs(min(tmp.dat$conf.low)) + 1
tmp.dat$Tw_complete <- tmp.dat$estimate + abs(min(tmp.dat$conf.low))
tmp.dat$Tw_constant_lci <- tmp.dat$conf.low + abs(min(tmp.dat$conf.low)) + 1
tmp.dat$Tw_complete_lci <- tmp.dat$conf.low + abs(min(tmp.dat$conf.low))
tmp.dat$Tw_constant_uci <- tmp.dat$conf.high + abs(min(tmp.dat$conf.low)) + 1
tmp.dat$Tw_complete_uci <- tmp.dat$conf.high + abs(min(tmp.dat$conf.low))

# Now merge with predict.dat
tmp.dat <- dplyr::select(tmp.dat, order, Tw_constant, Tw_complete, Tw_constant_lci, Tw_constant_uci, Tw_complete_lci, Tw_complete_uci)

predict.dat <- merge(predict.dat, tmp.dat, by="order", all.x = T)
head(predict.dat)



# Now calculate "N" (the number of observations upon which 
# each parameter estimate is based) for Tw
Tw.sum <- ddply(pan.dat, .(order), summarise, N_Tw=length(species))
predict.dat <- merge(predict.dat, Tw.sum, by="order", all.x = T)

head(predict.dat)
# Now, plot it

# Constant tolerance (Tw)
p4a <- ggplot(data=predict.dat) + 
  geom_hline(aes(yintercept=1.5), linetype=2) +
  geom_errorbar(aes(x=order, ymin=Tw_constant_lci, ymax=Tw_constant_uci, color=order), width=0, show.legend = F) +
  geom_point(aes(x=order, y=Tw_constant, fill=order, size=N_Tw), pch=21) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(.1,.1,0,.5), "lines"), 
        legend.position = c(.82,.91),
        panel.grid = element_blank(),
        legend.direction = "horizontal", legend.title = element_blank()) +
  scale_fill_manual(values=colz, guide="none") +
  scale_color_manual(values=colz) +
  ylab(bquote("constant immunopathology tolerance,"~T[w]))

# Complete tolerance (Tw)
p4b <- ggplot(data=predict.dat) + 
  geom_hline(aes(yintercept=.5), linetype=2) +
  geom_errorbar(aes(x=order, ymin=Tw_complete_lci, ymax=Tw_complete_uci, color=order), width=0, show.legend = F) +
  geom_point(aes(x=order, y=Tw_complete, fill=order, size=N_Tw), pch=21, show.legend = F) + 
  theme_bw() +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz) +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,.1,.1,.5), "lines"),
        panel.grid = element_blank()) +
  ylab(bquote("complete immunopathology tolerancem,"~T[w]))

# Visualize:
p4 <- cowplot::plot_grid(p4a, p4b, ncol=1, nrow = 2, labels=c("A", "B"), rel_heights = c(1,1.3), label_x = -0.01)



# Significant positive associations with tolerance: 
# Chiroptera, Cingulata, Monotremata, Primates

# Significant negative associations with tolerance: 
# Afrosoricida, Dasyuromorphia, Didelphimorphia,
# Eulipotyphla, Notoryctemorphia, Peramelemorphia 



######################################################
######################################################
## 3. Now estimate g0

# We know that constitutive immunity (baseline neutrophil count)
# scales positively with (a) body mass, # (b) negatively with BMR
# and mass-specific BMR, and positively with (c) higher host
# exposure/ promiscuity/ gregariousness

# We scale constitutive immunity by neutrophils concentrations
# reported in the Species360 data

# First, load the data
wbc.dat <- read.csv(file = "species_360_neut.csv", header=TRUE, stringsAsFactors = FALSE)
head(wbc.dat)

names(wbc.dat)[names(wbc.dat)=="phylo"] <- "binomial"

unique(wbc.dat$Order)

#Rename
names(wbc.dat)[names(wbc.dat)=="Segmented.neutrophils"] <- "neutro_conc"
#names(wbc.dat)[names(wbc.dat)=="AdultBodyMass"] <- "mass_g"
names(wbc.dat)[1:4] <- c("order", "family", "genus", "species")

#and merge with features from pan.dat
pan.bmr <- dplyr::select(pan.dat, binomial, mass_g, log10mass_g, max_lifespan_yrs, log10_max_lifespan_yrs, BMR_W, BMR_W_g)#, pop_group_size, pop_density_N_km2, homerange_km2)

wbc.dat <- merge(wbc.dat, pan.bmr, by="binomial", all.x = T)

head(wbc.dat)
wbc.dat$ln10neutro <- log10(wbc.dat$neutro_conc)
length(wbc.dat$mass_g[is.na(wbc.dat$mass_g)]) #102 with no mass
length(wbc.dat$mass_g[is.na(wbc.dat$max_lifespan_yrs)]) #102 with no lifespan
length(wbc.dat$mass_g[is.na(wbc.dat$BMR_W)]) #323 with no BMR


# Slim your colors down to only those in this dataset
colz2 = sort(colz[unique(wbc.dat$order)])

wbc.dat$order <- as.factor(wbc.dat$order)

# Plot with mass
p5 <- ggplot(wbc.dat) + 
  geom_point(aes(x=mass_g, y=neutro_conc, fill=order), size =3, pch=21) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank()) +
  xlab("mass (g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
p5



# BMR plot 


p6 <- ggplot(wbc.dat) + 
  geom_point(aes(x=BMR_W_g, y=neutro_conc, fill=order), size =3, pch=21) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank()) +
  xlab("mass-specific BMR (W/g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
p6





# plot with predictions
# faux gam first
#wbc.dat$ln10BMR_W_g <- log10(wbc.dat$mass_g)/wbc.dat$BMR_W
wbc.dat$ln10BMR_W_g <- log10(wbc.dat$BMR_W_g)

wbc.dat$order <- as.factor(wbc.dat$order)


m2a <- lmer(ln10neutro~ln10BMR_W_g + (1|order), data=wbc.dat)

plot_model(m2a, type="re")
# m2a <- gam(ln10neutro~s(ln10BMR_W_g, bs="tp") +
#             s(order, bs="re"),
#           
# summary(m2a) #51.8%; n =144

wbc.pred = cbind.data.frame(ln10BMR_W_g=seq(min(wbc.dat$ln10BMR_W_g, na.rm=T),  max(wbc.dat$ln10BMR_W_g, na.rm=T), length=2),  order = "Chiroptera")
#wbc.pred$predict_neut <- 10^(predict.gam(m2a, newdata = wbc.pred, exclude = "s(order)"))
wbc.pred$predict_neut <- 10^(predict(m2a, newdata = wbc.pred, re.form=NA))
wbc.pred$BMR_W_g <- 10^(wbc.pred$ln10BMR_W_g)


p6b <- ggplot(wbc.dat) + 
  geom_point(aes(x=BMR_W_g, y=neutro_conc, fill=order), size =3, pch=21) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  geom_line(data = wbc.pred, aes(x=BMR_W_g , y=predict_neut), size=1) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank()) +
  xlab("mass-specific BMR (W/g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
p6b



wbc.dat$order <- as.factor(wbc.dat$order)
wbc.dat$ln10mass <- log10(wbc.dat$mass_g)
# 
# m2 <- lmer(ln10neutro~ln10mass + (1|order), data=wbc.dat)
# # m2 <- gam(ln10neutro~s(ln10mass, bs="tp") +
# #                         s(order, bs="re"),
# #                         data=wbc.dat)
# summary(m2) #41.8%; n =365
# 
# plot_model(m2, type="re")
# #order.dat <- get_partial_effects(m2, var="order")
# #plot.partial(order.dat, var="order", response_var = "log10 neutrophil conc.")



# now look at mass and BMR independently (model is stronger)
library(lmerTest)
m3 <- lmer(ln10neutro~ln10mass+BMR_W + (1|order), data=wbc.dat)

summary(m3) 
AIC(m2a, m3) #m3 best



#and finally, calculate g0 as the order effects from this model

tmp.dat = cbind.data.frame(order= plot_model(m3, type="re")$data$term, estimate=plot_model(m3, type="re")$data$estimate, conf.low=plot_model(m3, type="re")$data$conf.low, conf.high=plot_model(m3, type="re")$data$conf.high) 



# tmp.dat <- cbind.data.frame(order= order.dat[[1]]$order, estimate=order.dat[[1]]$y, conf.low=order.dat[[1]]$ylower, conf.high=order.dat[[1]]$yupper)

# Our goal is for g0 to span 0 to 1 for model parameterization.
# We only allow for linear transformations of the data, in order 
# to retain differences in the magnitude of effect. We add to 
# the predict.dat database from above among the orders. 


# Linear transformation: 
# g0: Make all effects positive
tmp.dat$g0 <- tmp.dat$estimate + abs(min(tmp.dat$conf.low)) +.0000001
tmp.dat$g0_lci <- tmp.dat$conf.low + abs(min(tmp.dat$conf.low))+.0000001
tmp.dat$g0_uci <- tmp.dat$conf.high + abs(min(tmp.dat$conf.low))+.0000001

#and rescale to fall within a range that causes normalized impact
tmp.dat$g0 <-  scales::rescale(x =tmp.dat$g0, from=c(min(tmp.dat$g0_lci), max(tmp.dat$g0_uci)), to =c(.0000001,.1)) 
tmp.dat$g0_lci <-  scales::rescale(x =tmp.dat$g0_lci, from=c(min(tmp.dat$g0_lci), max(tmp.dat$g0_uci)), to =c(.0000001,.1)) 
tmp.dat$g0_uci <-  scales::rescale(x =tmp.dat$g0_uci, from=c(min(tmp.dat$g0_lci), max(tmp.dat$g0_uci)), to =c(.0000001,.1)) 


# Now merge with predict.dat
tmp.dat <- dplyr::select(tmp.dat, order, g0, g0_lci, g0_uci)


predict.dat <- merge(predict.dat, tmp.dat, by="order", all.x = T)
head(predict.dat)



# Now calculate "N" (the number of observations upon which 
# each parameter estimate is based) for g0
g0.sum <- ddply(wbc.dat, .(order), summarise, N_g0=length(binomial))
predict.dat <- merge(predict.dat, g0.sum, by="order", all.x = T)

head(predict.dat)




# g0
p8 <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  geom_hline(aes(yintercept=0.05), linetype=2) +
  geom_point(aes(x=order, y=g0, fill=order, size=N_g0), pch=21) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),
        plot.margin = unit(c(.1,.1,.1,1.1), "lines"),
        panel.grid = element_blank(),
        legend.position = c(.75,.92),
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  geom_errorbar(aes(x=order, ymin=g0_lci, ymax=g0_uci, color=order), width=0, show.legend = F) +
  scale_color_manual(values=colz2) +
  scale_fill_manual(values=colz2, guide="none") +
  ylab(bquote("magnitude constitutive immunity, "~g[0]))
p8
######################################################
######################################################
## 4. Now, estimate Tvs

# Finally, and add in phylogenetic distance from Primates from MRCA tree,
# one entry unique to each order

# Load the tree data
load(paste0(homewd, "phylo-tree/phylo.dat.final.Rdata"))

# Merge the phylogenetic distance (eta_R) into the existing tree
predict.dat <- merge(predict.dat, phylo.dat, by="order", all.x = T, sort=F)

head(predict.dat)

# Now, estimate Tv_human (Tv_S), using similar scaling as above for Tw:
# Constant tolerance needs to be >1 and 
# Complete tolerance needs to be >0 and <1

# Here, small phylogenetic distance equates to high tolerance, so we
# subtract from 2 and 1 rather than adding:
predict.dat$Tv_human_constant <- 2-(predict.dat$phylo_dist/max(predict.dat$phylo_dist, na.rm=T))
predict.dat$Tv_human_complete <- 1-(predict.dat$phylo_dist/max(predict.dat$phylo_dist, na.rm=T))

# Visualize both constant (p4a) and complete (p4b) estimates for Tvs
p9a <- ggplot(data=predict.dat) + theme_bw() +
        geom_point(aes(x=order, y=Tv_human_constant, fill=order), pch=21, show.legend = F, size=3) + 
        scale_fill_manual(values=colz) + ylab(bquote("constant viral spillover tolerance,"~T[vS])) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              plot.margin = unit(c(.1,.1,0,.5), "lines"),
              panel.grid = element_blank()) +
        geom_hline(aes(yintercept=1.5), linetype=2)
p9b <- ggplot(data=predict.dat) + theme_bw() +
       geom_point(aes(x=order, y=Tv_human_complete, fill=order), pch=21, show.legend = F, size=3) + 
       scale_fill_manual(values=colz) + ylab(bquote("complete viral spillover tolerance,"~T[vS])) +
       theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
             plot.margin = unit(c(0,.1,.1,.5), "lines"),
             panel.grid = element_blank()) +
      geom_hline(aes(yintercept=0.5), linetype=2)

p9 <- cowplot::plot_grid(p9a, p9b, ncol=1, nrow = 2, labels=c("A", "B"), rel_heights = c(1,1.3), label_x = -0.01)




######################################################
######################################################
## 5. Now, estimate r*

predict.dat <- read.csv(file=paste0(homewd, subwd, "/predict_pars.csv"), header = T, stringsAsFactors = F)
head(predict.dat)


# Now take all the outputs and merge them into a prediction for rstar and for alphastar-hum
# Set all other within-host parameters to the same default parameter values
# used to generate Fig. 2 in the main text.
# The values for g0, Tw, Tv, and mu will get overwritten, so it does not 
# really matter what you put here.

vir.par <- list(b= .2, q = .0002, c=.5,  m=1/(21),  g=.9, g0 = .3, zeta=.2,
                v=1, w=1,  Tv=.005, Tw=.005,mu = 1/(20*365))


# Constant tolerance rstar predictions using order-specific g0, Tw, and mu
# and default values for all other parameters

vir.par$Tv = 1.005
predict.dat$rstar_constant <- (vir.par$c*predict.dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu*predict.dat$Tw_constant*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))
predict.dat$rstar_constant_lci <- (vir.par$c*predict.dat$g0_lci/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0_lci*predict.dat$mu_lci*predict.dat$Tw_constant_lci*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant_lci+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant_lci) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))
predict.dat$rstar_constant_uci <- (vir.par$c*predict.dat$g0_uci/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0_uci*predict.dat$mu_uci*predict.dat$Tw_constant_uci*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant_uci+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant_uci) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))


# Complete rstar predictions using  order-specific g0, Tw, and mu
# and default values for all other parameters
vir.par$Tv=.005

predict.dat$rstar_complete <- (vir.par$c*predict.dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu*vir.par$g*predict.dat$g0*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete-vir.par$Tv)))
predict.dat$rstar_complete_lci <- (vir.par$c*predict.dat$g0_lci/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0_lci*predict.dat$mu_lci)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu_lci*vir.par$g*predict.dat$g0_lci*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete_lci-vir.par$Tv)))
#predict.dat$rstar_complete_lci[is.na(predict.dat$rstar_complete_lci)] <- 0
predict.dat$rstar_complete_uci <- (vir.par$c*predict.dat$g0_uci/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0_uci*predict.dat$mu_uci)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu_uci*vir.par$g*predict.dat$g0_uci*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete_uci-vir.par$Tv)))


# First, get mean N across all the factors that went in to each prediction
predict.dat$N_cumulative <-rowMeans(cbind(predict.dat$N_mu, predict.dat$N_Tw, predict.dat$N_g0))
predict.dat <- arrange(predict.dat, desc(rstar_constant), desc(N_cumulative))
predict.dat$order <- factor(predict.dat$order, levels = unique(predict.dat$order))

# Visualize rstar by order

#Constant tolerance
p10a <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  geom_errorbar(aes(x=order, ymin=rstar_constant_lci, ymax=rstar_constant_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(x=order, y=rstar_constant, fill=order, size=N_cumulative), pch=21) +
  scale_color_manual(values=colz)+
  scale_fill_manual(values=colz, guide="none")+theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(), legend.direction = "horizontal",
        legend.position = c(.75,.9),
        plot.margin = unit(c(.1,.1,0,.8), "lines")) +
  ylab(bquote(atop("optimal virus growth rate", "in reservoir,"~r~"* (constant tolerance)")))


p10b <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  geom_errorbar(aes(x=order, ymin=rstar_complete_lci, ymax=rstar_complete_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(x=order, y=rstar_complete, fill=order, size=N_cumulative), pch=21, show.legend = F) +
  scale_color_manual(values=colz)+
  scale_fill_manual(values=colz, guide="none")+theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.8), "lines")) +
  ylab(bquote(atop("optimal virus growth rate", "in reservoir,"~r~"* (complete tolerance)")))


p10 <- cowplot::plot_grid(p10a, p10b, ncol=1, nrow = 2, labels=c("A", "B"), rel_heights = c(1,1.3), label_x = -0.01, align = "v")


######################################################
######################################################
## 5. Now, estimate alpha* in spillover human host

# Virulence in the human host is a function of phylogenetic distance,
# represented in order-specific Tvs (shown here as Tv_human_constant)
# Tolerance of immunopathology is now human, so we hold constant across 
# all reservoir host orders.

# First, calculate Vmax in the human host:
predict.dat$Vs_max_constant <- predict.dat$rstar_constant/(vir.par$g*vir.par$c) - predict.dat$rstar_constant/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_constant*vir.par$g)
predict.dat$Vs_max_constant_lci <- predict.dat$rstar_constant_lci/(vir.par$g*vir.par$c) - predict.dat$rstar_constant_lci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_constant_lci*vir.par$g)
predict.dat$Vs_max_constant_uci <- predict.dat$rstar_constant_uci/(vir.par$g*vir.par$c) - predict.dat$rstar_constant_uci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_constant_uci*vir.par$g)

predict.dat$Vs_max_complete <- predict.dat$rstar_complete/(vir.par$g*vir.par$c) - predict.dat$rstar_complete/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_complete*vir.par$g)
predict.dat$Vs_max_complete_lci <- predict.dat$rstar_complete_lci/(vir.par$g*vir.par$c) - predict.dat$rstar_complete_lci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_complete_lci*vir.par$g)
predict.dat$Vs_max_complete_uci <- predict.dat$rstar_complete_uci/(vir.par$g*vir.par$c) - predict.dat$rstar_complete_uci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_complete_uci*vir.par$g)


# Next, take that viral load and calculate spillover virulence, here for the constant tolerance assumption:
vir.par$Tw = 1
predict.dat$alpha_star_human_constant <- ((predict.dat$rstar_constant*vir.par$v)/predict.dat$Tv_human_constant  + (vir.par$g*vir.par$w*predict.dat$rstar_constant)/vir.par$Tw)*predict.dat$Vs_max_constant
predict.dat$alpha_star_human_constant_lci <- ((predict.dat$rstar_constant_lci*vir.par$v)/predict.dat$Tv_human_constant  + (vir.par$g*vir.par$w*predict.dat$rstar_constant_lci)/vir.par$Tw)*predict.dat$Vs_max_constant_lci
predict.dat$alpha_star_human_constant_uci <- ((predict.dat$rstar_constant_uci*vir.par$v)/predict.dat$Tv_human_constant  + (vir.par$g*vir.par$w*predict.dat$rstar_constant_uci)/vir.par$Tw)*predict.dat$Vs_max_constant_uci


# And here for complete tolerance:
vir.par$Tw = 0
predict.dat$alpha_star_human_complete <- (predict.dat$rstar_complete*(vir.par$v-predict.dat$Tv_human_complete) + predict.dat$rstar_complete*vir.par$g*(vir.par$w-vir.par$Tw))*predict.dat$Vs_max_complete
predict.dat$alpha_star_human_complete_lci <- (predict.dat$rstar_complete_lci*(vir.par$v-predict.dat$Tv_human_complete) + predict.dat$rstar_complete_lci*vir.par$g*(vir.par$w-vir.par$Tw))*predict.dat$Vs_max_complete_lci
predict.dat$alpha_star_human_complete_uci <- (predict.dat$rstar_complete_uci*(vir.par$v-predict.dat$Tv_human_complete) + predict.dat$rstar_complete_uci*vir.par$g*(vir.par$w-vir.par$Tw))*predict.dat$Vs_max_complete_uci

# Now, we rank by virulence, then confidence, and plot...
predict.dat <- arrange(predict.dat, desc(alpha_star_human_constant), desc(N_cumulative))
predict.dat$order <- factor(predict.dat$order, levels = unique(predict.dat$order))

# Visualize alphastar in humans

# Constant
p11a <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  #geom_errorbar(aes(x=order, ymin=alpha_star_human_constant_lci, ymax=alpha_star_human_constant_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(x=order, y=alpha_star_human_constant, fill=order, size=N_cumulative), pch=21) +
  scale_color_manual(values=colz)+
  scale_fill_manual(values=colz, guide="none")+theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(), legend.direction = "horizontal",
        legend.position = c(.75,.9),
        plot.margin = unit(c(.1,.1,0,.8), "lines")) +
  ylab(bquote("spillover virulence,"~alpha[S]~"(constant tolerance)"))
  

predict.dat <- arrange(predict.dat, desc(alpha_star_human_constant), desc(N_cumulative))
predict.dat$order <- factor(predict.dat$order, levels = unique(predict.dat$order))

# Complete
p11b <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  #geom_errorbar(aes(x=order, ymin=alpha_star_human_complete_lci, ymax=alpha_star_human_complete_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(x=order, y=alpha_star_human_complete, fill=order, size=N_cumulative), pch=21, show.legend = F) +
  scale_color_manual(values=colz)+
  scale_fill_manual(values=colz, guide="none")+theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.5), "lines")) +
  ylab(bquote("spillover virulence,"~alpha[S]~"(complete tolerance)"))

# Visualize together
p11 <- cowplot::plot_grid(p11a, p11b, ncol=1, nrow = 2, labels=c("A", "B"), rel_heights = c(1,1.3), label_x = -0.01)






######################################################
######################################################
## 6. Now, compare with the literature


load(paste0(homewd,"source/gam.dat.Guth.et.al.2021.Rdata"))

# Plot alpha_S for real
p12 <- ggplot(data=gam.dat) + 
  geom_errorbar(aes(x=hOrder, ymin=alpha_lci, ymax=alpha_uci, color=hOrder),  width=0, linetype=3, show.legend=F) +
  geom_point(aes(hOrder, alpha, color=hOrder, size=Nobs), show.legend=F) + 
  theme_bw() +
  scale_color_manual(values=colz) +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12, angle = 90), legend.title = element_blank()) + 
  ylab(bquote("predicted"~alpha[S]~"by order ("~days^-1~")")) #+ coord_cartesian(ylim=c(0,100))

p12



#rank by descending alpha, then confidence
gam.dat <- arrange(gam.dat, desc(alpha), desc(Nobs))


# Merge our nested model predictions with those from Guth et al. 2022

# Make columns match
head(gam.dat)
names(gam.dat)[1] <- "order"
names(gam.dat)[3] <- "N"
names(predict.dat)[names(predict.dat)=="N_cumulative"] <- "N"
head(predict.dat)


# Before this, save for fitting!

sub.dat <- dplyr::select(predict.dat, order, mu, Tw_constant, Tw_complete, g0, Tv_human_constant, Tv_human_complete, alpha_star_human_constant, alpha_star_human_complete)

cons.sub = dplyr::select(sub.dat, order, mu, Tw_constant, g0, Tv_human_constant, alpha_star_human_constant)
comp.sub  = dplyr::select(sub.dat, order, mu, Tw_complete, g0, Tv_human_complete, alpha_star_human_complete)

cons.sub$tolerance="constant"
comp.sub$tolerance="complete"

names(cons.sub) <- names(comp.sub) <- c("order", "mu", "Tw", "g0", "Tv_human", "alpha_star_human", "tolerance")
cons.all <- rbind(cons.sub, comp.sub)

#and merge with gam
gam.merge <- dplyr::select(gam.dat, order, alpha)


dat.comp <- merge(cons.all, gam.merge, by="order")
dat.comp$alpha_star_human[dat.comp$tolerance=="constant"]
dat.comp$alpha_star_human[dat.comp$tolerance=="complete"]

#rescale alpha lit by range of our predictions
dat.comp$alpha[dat.comp$tolerance=="constant"] <- scales::rescale(x=dat.comp$alpha[dat.comp$tolerance=="constant"], from=c(min(dat.comp$alpha[dat.comp$tolerance=="constant"]), max(dat.comp$alpha[dat.comp$tolerance=="constant"])), to =c(min(dat.comp$alpha_star_human[dat.comp$tolerance=="constant"]), max(dat.comp$alpha_star_human[dat.comp$tolerance=="constant"])))
dat.comp$alpha[dat.comp$tolerance=="complete"] <- scales::rescale(x=dat.comp$alpha[dat.comp$tolerance=="complete"], from=c(min(dat.comp$alpha[dat.comp$tolerance=="complete"]), max(dat.comp$alpha[dat.comp$tolerance=="complete"])), to =c(min(dat.comp$alpha_star_human[dat.comp$tolerance=="complete"]), max(dat.comp$alpha_star_human[dat.comp$tolerance=="complete"])))

dat.comp$alpha[dat.comp$tolerance=="constant"]
dat.comp$alpha[dat.comp$tolerance=="complete"]


#check it
plot(x=dat.comp$alpha[dat.comp$tolerance=="constant"], y=dat.comp$alpha_star_human[dat.comp$tolerance=="constant"])
plot(x=dat.comp$alpha[dat.comp$tolerance=="complete"], y=dat.comp$alpha_star_human[dat.comp$tolerance=="complete"])



#and the general fitting par
dat=dat.comp

dat$b=.2
dat$q=.0002
dat$c=.5
dat$m=1/21
dat$g1=.9
dat$zeta=.2
dat$v=1
dat$w=1
dat$Tv_bat=NA
dat$Tv_bat[dat$tolerance=="complete"] <-  .005
dat$Tv_bat[dat$tolerance=="constant"] <-  1.005

head(dat)
names(dat)[names(dat)=="alpha"] <- "lit_alpha_scale"

write.csv(dat, file=paste0(homewd,"figure-3/fit-par-dat.csv"), row.names = F)






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
plot.dat <- plot.dat[complete.cases(plot.dat),]

shapez <- c("predicted from zoonoses"=25, "predicted from\nnested model"=24)


# And add images by order
plot.dat <- dplyr::select(plot.dat, -(N), -(source))
dat.comp =subset(plot.dat, tolerance!="natural")
dat.nat = subset(plot.dat, tolerance=="natural")
head(dat.comp)
head(dat.nat)

dat.nat$alpha= -1*dat.nat$alpha
dat.nat$alpha_lci= -1*dat.nat$alpha_lci
dat.nat$alpha_uci= -1*dat.nat$alpha_uci
names(dat.nat) <- c("order", "lit_alpha", "lit_alpha_lci", "lit_alpha_uci", "tolerance")
dat.nat <- dplyr::select(dat.nat, -(tolerance)) 

dat.comp <- merge(dat.comp, dat.nat, by="order")

#add in from above
head(predict.dat)

