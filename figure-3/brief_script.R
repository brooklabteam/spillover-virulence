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

colz = scales::hue_pal()(length(unique(pan.dat$order))-1)
colz=c(colz, "red")

names(colz) <- c(unique(pan.dat$order)[unique(pan.dat$order)!="Chiroptera"], "Chiroptera")


p1 <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,   fill=order), size =3, pch=21) +  
  scale_fill_manual(values=colz) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  theme(panel.grid = element_blank()) +
  xlab("mass (g)") + ylab("max lifespan (yrs)")
p1

# Save plot

ggsave(file = paste0(homewd, subwd, "/brief/Fig1.png"),
       plot = p1,
       units="mm",  
       width=60, 
       height=30, 
       scale=3, 
       dpi=300)


# Now, run a GAM to determine the partial effect of host order
# (phylogeny) on maximum longevity, and predict lifespan 
# rates by order:

library(mgcv)
pan.dat$log10mass_g <- log10(pan.dat$mass_g)
pan.dat$order <- as.factor(pan.dat$order)
pan.dat$log10_max_lifespan_yrs <- log10(pan.dat$max_lifespan_yrs)

m1 <- gam(log10_max_lifespan_yrs~s(log10mass_g, bs = "tp") +
                                  s(order, bs="re"),
                                  data = pan.dat)

summary(m1) #deviance explained = 75.8%


# And, plot the partial effects of order:
# Source partial effects script from Mollentze and Streicker 2020:
# This script includes two plotting functions, added by me!

source(paste0(homewd,"source/mollentze-streicker-2020-functions.R"))


order.dat <- get_partial_effects(m1, var="order")
mass.dat <- get_partial_effects_continuous(m1, var="log10mass_g")

p2a <- plot.partial(df=order.dat, var="order", response_var = "max lifespan (yrs)")
p2b <- plot.partial.cont(df=mass.dat, var="log10mass_g",  response_var = "max lifespan (yrs)", alt_var = "mass (g)", log = T)

p2 <- cowplot::plot_grid(p2a, p2b, ncol=1, nrow = 2, labels=c("A", "B"), label_x = .1)

ggsave(file = paste0(homewd, subwd, "/brief/Fig2.png"),
       plot = p2,
       units="mm",  
       width=25, 
       height=40, 
       scale=5, 
       dpi=300)



# Significant positive associations with lifespan: 
# Chiroptera, Cingulata, Monotremata, Primates

# Significant negative associations with lifespan: 
# Afrosoricida, Dasyuromorphia, Didelphimorphia,
# Eulipotyphla, Notoryctemorphia, Peramelemorphia 


# And plot with predictions:

pan.dat.sum <- cbind.data.frame(log10mass_g=seq(log10(min(pan.dat$mass_g)), max(log10(pan.dat$mass_g)), length=10), order = "Chiroptera")
pan.dat.sum$predict_lifespan <- 10^(predict.gam(m1, newdata = pan.dat.sum, exclude = "s(order)"))
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

ggsave(file = paste0(homewd, subwd, "/brief/Fig1w_line.png"),
       plot = p1b,
       units="mm",  
       width=60, 
       height=30, 
       scale=3, 
       dpi=300)


# Now, to parameterize mu, the annual mortality rate by order for
# our within-host model, we can simply predict lifespan from our fitted
# GAM, excluding the effects of body mass, then take the inverse (but be sure to
# express in timesteps of days to make comparable to the viral dynamics):

predict.dat <- cbind.data.frame(order = order.dat[[1]]$order, log10mass_g=1)
# You can insert any value you like for "log10mass_g" since we won't 
# be using it anyhow

predict.dat$mu <- 1/((10^predict.gam(m1, newdata = predict.dat, exclude = "s(log10mass_g)", type="response"))*365)
predict.dat$mu_lci <- 1/((10^(predict.gam(m1, newdata = predict.dat, exclude = "s(log10mass_g)", type="response") -1.96*predict.gam(m1, newdata = predict.dat, exclude = "s(log10mass_g)", type="response", se.fit = T)$se))*365)
predict.dat$mu_uci <- 1/((10^(predict.gam(m1, newdata = predict.dat, exclude = "s(log10mass_g)", type="response") +1.96*predict.gam(m1, newdata = predict.dat, exclude = "s(log10mass_g)", type="response", se.fit = T)$se))*365)

predict.dat$mu_lci[predict.dat$mu_lci<0] <- 0
predict.dat$mu_uci[predict.dat$mu_uci>1] <- 1

# Also compute the number of data entries per order used to determine this:
mu.sum <- ddply(pan.dat, .(order), summarise, N_mu = length(binomial))

predict.dat <- merge(predict.dat, mu.sum, by="order", all.x = T)
predict.dat <- dplyr::select(predict.dat, -(log10mass_g))


# and plot your predictions by order for mass

# first, get your null
y.int = 1/((10^(predict.gam(m1, 
                    newdata = cbind.data.frame(order="Primates", log10mass_g = unique(order.dat[[1]]$log10mass_g)), 
                    exclude = "s(order)", type = "response")))*365)

# Now plot predictions against null:

p3 <- ggplot(data=predict.dat) + 
      geom_point(aes(x=order, y=mu, fill=order, size=N_mu), pch=21) + 
      theme_bw()  + 
      theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
            legend.position = c(.8,.91), legend.direction = "horizontal", 
            legend.title = element_blank(),
            panel.grid = element_blank()) +
      scale_fill_manual(values=colz, guide="none") +
      scale_color_manual(values=colz) +
      geom_errorbar(aes(x=order, ymin=mu_lci, ymax=mu_uci, color=order), width=0, show.legend = F) +
      geom_hline(aes(yintercept=y.int), linetype=2) +
      ylab(bquote("predicted annual mortality rate by order,"~mu~"("~days^-1~")")) 

p3

# And save

ggsave(file = paste0(homewd, subwd, "/brief/Fig3.png"),
       plot = p3,
       units="mm",  
       width=45, 
       height=40, 
       scale=3, 
       dpi=300)


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

tmp.dat <- cbind.data.frame(order= order.dat[[1]]$order, estimate=order.dat[[1]]$y, conf.low=order.dat[[1]]$ylower, conf.high=order.dat[[1]]$yupper)

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
  geom_errorbar(aes(x=order, ymin=Tw_constant_lci, ymax=Tw_constant_uci, color=order), width=0, show.legend = F) +
  geom_hline(aes(yintercept=1.5), linetype=2) + ylab(bquote("constant immunopathology tolerance,"~T[w]))

# Complete tolerance (Tw)
p4b <- ggplot(data=predict.dat) + 
  geom_point(aes(x=order, y=Tw_complete, fill=order, size=N_Tw), pch=21, show.legend = F) + 
  theme_bw() +
  scale_fill_manual(values=colz) +
  scale_color_manual(values=colz) +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,.1,.1,.5), "lines"),
        panel.grid = element_blank()) +
  geom_errorbar(aes(x=order, ymin=Tw_complete_lci, ymax=Tw_complete_uci, color=order), width=0, show.legend = F) +
  geom_hline(aes(yintercept=.5), linetype=2) + ylab(bquote("complete immunopathology tolerancem,"~T[w]))

# Visualize:
p4 <- cowplot::plot_grid(p4a, p4b, ncol=1, nrow = 2, labels=c("A", "B"), rel_heights = c(1,1.3), label_x = -0.01)


ggsave(file = paste0(homewd, subwd, "/brief/Fig4.png"),
       plot = p4,
       units="mm",  
       width=40, 
       height=70, 
       scale=3, 
       dpi=300)

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

# We scale by (a) after Cynthia Downs' 2020 AmNat paper 

# First, load the data
wbc.dat <- read.csv(file = "Downs_et_al_2020.csv", header=TRUE, stringsAsFactors = FALSE)
head(wbc.dat)

# Drop unnecessary columns
wbc.dat <- dplyr::select(wbc.dat, Order, Family, Genus, Species, Common.Name, 
                      AdultBodyMass, TotalWBC, Lymphocytes, Segmented.neutrophils,
                      ln10Mass, ln10WBC, ln10lympho, ln10neutro) 

head(wbc.dat) 

wbc.dat$binomial <- paste0(wbc.dat$Genus, "_", wbc.dat$Species)

wbc.dat <- dplyr::select(wbc.dat, Order, Family, Genus, Species, binomial, AdultBodyMass, Segmented.neutrophils, ln10Mass, ln10WBC, ln10lympho, ln10neutro)

# Reclassify the Afrotherians and the Xenartha
wbc.dat$Order <- as.character(wbc.dat$Order)
wbc.dat$binomial[wbc.dat$Order=="Afrotheria"]
wbc.dat$binomial[wbc.dat$Order=="Xenarthra"]
wbc.dat$Order[wbc.dat$binomial=="Elephas_maximus"] <- "Proboscidea"
wbc.dat$Order[wbc.dat$binomial=="Loxodonta_africana"]<- "Proboscidea"
wbc.dat$Order[wbc.dat$binomial=="Orycteropus_afer"]<- "Tubulidentata"
wbc.dat$Order[wbc.dat$binomial=="Procavia_capensis"] <- "Hyracoidea"
wbc.dat$Order[wbc.dat$binomial=="Trichechus_manatus"] <- "Sirenia"
wbc.dat$Order[wbc.dat$binomial=="Dasypus_novemcinctus"] <- "Cingulata"

#Rename
names(wbc.dat)[names(wbc.dat)=="Segmented.neutrophils"] <- "neutro_conc"
names(wbc.dat)[names(wbc.dat)=="AdultBodyMass"] <- "mass_g"
names(wbc.dat)[1:4] <- c("order", "family", "genus", "species")


#and link to other pan.dat features: BMR and exposure and pop size
pan.bmr <- dplyr::select(pan.dat, binomial, BMR_W, BMR_W_g, pop_group_size, pop_density_N_km2, homerange_km2)

wbc.dat <- merge(wbc.dat, pan.bmr, by="binomial", all.x = T)

# Slim your colors down to only those in this dataset
colz2 = colz[unique(wbc.dat$order)]

# Plot with mass
p5 <- ggplot(wbc.dat) + 
  geom_point(aes(x=mass_g, y=neutro_conc, fill=order), size =3, pch=21) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank()) +
  xlab("mass (g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
p5


# Save plot
ggsave(file = paste0(homewd, subwd, "/brief/Fig5.png"),
       plot = p5,
       units="mm",  
       width=60, 
       height=30, 
       scale=3, 
       dpi=300)


# BMR plot 


p6 <- ggplot(wbc.dat) + 
  geom_point(aes(x=BMR_W_g, y=neutro_conc, fill=order), size =3, pch=21) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank()) +
  xlab("mass-specific BMR (W/g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
p6

ggsave(file = paste0(homewd, subwd, "/brief/Fig6.png"),
       plot = p6,
       units="mm",  
       width=60, 
       height=30, 
       scale=3, 
       dpi=300)



# plot with predictions
# faux gam first
wbc.dat$ln10BMR_W_g <- log10(wbc.dat$mass_g)/wbc.dat$BMR_W
wbc.dat$ln10BMR_W_g <- log10(wbc.dat$BMR_W_g)

wbc.dat$order <- as.factor(wbc.dat$order)

m2a <- gam(ln10neutro~s(ln10BMR_W_g, bs="tp") +
            s(order, bs="re"),
          data=wbc.dat)
summary(m2a) #54.2%; n =97

wbc.pred = cbind.data.frame(ln10BMR_W_g=seq(min(wbc.dat$ln10BMR_W_g, na.rm=T),  max(wbc.dat$ln10BMR_W_g, na.rm=T), length=2),  order = "Chiroptera")
wbc.pred$predict_neut <- 10^(predict.gam(m2a, newdata = wbc.pred, exclude = "s(order)"))
wbc.pred$BMR_W_g <- 10^(wbc.pred$ln10BMR_W_g)


p6b <- ggplot(wbc.dat) + 
  geom_point(aes(x=BMR_W_g, y=neutro_conc, fill=order), size =3, pch=21) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  geom_line(data = wbc.pred, aes(x=BMR_W_g , y=predict_neut), size=1) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank()) +
  xlab("mass-specific BMR (W/g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
p6b

ggsave(file = paste0(homewd, subwd, "/brief/Fig6_w_line.png"),
       plot = p6b,
       units="mm",  
       width=50, 
       height=30, 
       scale=3, 
       dpi=300)

# Now, model effects of order as a GAM, using just mass
wbc.dat$order <- as.factor(wbc.dat$order)
m2 <- gam(ln10neutro~s(ln10Mass, bs="tp") +
                        s(order, bs="re"),
                        data=wbc.dat)
summary(m2) #47.1%; n =260
order.dat <- get_partial_effects(m2, var="order")
plot.partial(order.dat, var="order", response_var = "log10 neutrophil conc.")

# Significant negative associations with:
# Cetartiodactyal, Diprotodontia, Proboscidea

#Significant positive associations with:
# Monotremata, Primates



# Or try including BMR
m3 <- gam(ln10neutro~s(ln10Mass, bs="tp") +
                      s(BMR_W, bs="tp") +
                      s(order, bs="re"),
                      data=wbc.dat)
summary(m3) #69.3%; n=97. much stronger model fit
order.dat <- get_partial_effects(m3, var="order")
mass.dat <- get_partial_effects_continuous(m3, var="ln10Mass")
BMR.dat <- get_partial_effects_continuous(m3, var="BMR_W")

p7a <- plot.partial(order.dat, var="order", response_var = "log10 neutrophils")
p7b <- plot.partial.cont(mass.dat, var="ln10Mass", log = T, alt_var = "mass (g)", response_var = "log10 neutrophils")
p7c <- plot.partial.cont(BMR.dat, var="BMR_W", log = F, alt_var = "BMR (W)", response_var = "log10 neutrophils")


#and plot together

p7 <- cowplot::plot_grid(p7a, p7b, p7c, ncol=1, nrow = 3, labels=c("A", "B", "C"), label_x = 0.09)


ggsave(file = paste0(homewd, subwd, "/brief/Fig7.png"),
       plot = p7,
       units="mm",  
       width=50, 
       height=100, 
       scale=3, 
       dpi=300)


#and finally, calculate g0 as the order effects from this model


tmp.dat <- cbind.data.frame(order= order.dat[[1]]$order, estimate=order.dat[[1]]$y, conf.low=order.dat[[1]]$ylower, conf.high=order.dat[[1]]$yupper)

# Our goal is for g0 to span 0 to 1 for model parameterization.
# We only allow for linear transformations of the data, in order 
# to retain differences in the magnitude of effect. We add to 
# the predict.dat database from above among the orders. 


# Linear transformation: 
# g0: Make all effects positive
tmp.dat$g0 <- tmp.dat$estimate + abs(min(tmp.dat$conf.low))
tmp.dat$g0_lci <- tmp.dat$conf.low + abs(min(tmp.dat$conf.low))
tmp.dat$g0_uci <- tmp.dat$conf.high + abs(min(tmp.dat$conf.low))

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
  geom_hline(aes(yintercept=0.5), linetype=2) + ylab(bquote("magnitude constitutive immunity, "~g[0]))

ggsave(file = paste0(homewd, subwd, "/brief/Fig8.png"),
       plot = p8,
       units="mm",  
       width=45, 
       height=40, 
       scale=3, 
       dpi=300)



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



ggsave(file = paste0(homewd, subwd, "/brief/Fig9.png"),
       plot = p9,
       units="mm",  
       width=30, 
       height=60, 
       scale=3, 
       dpi=300)


######################################################
######################################################
## 5. Now, estimate r*

# Now take all the outputs and merge them into a prediction for rstar and for alphastar-hum
# Set all other within-host parameters to the same default parameter values
# used to generate Fig. 2 in the main text.
# The values for g0, Tw, Tv, and mu will get overwritten, so it does not 
# really matter what you put here.

vir.par <- list(b= .2, q = .0002, c=.5,  m=1/(21),  g=.9, g0 = .3, zeta=.2,
                v=1, w=1,  Tv=.005, Tw=.005,mu = 1/(20*365))


# Constant tolerance rstar predictions using order-specific g0, Tw, and mu
# and default values for all other parameters
predict.dat$rstar_constant <- (vir.par$c*predict.dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu*predict.dat$Tw_constant*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))
predict.dat$rstar_constant_lci <- (vir.par$c*predict.dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu_uci*predict.dat$Tw_constant_lci*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant_lci+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant_lci) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))
predict.dat$rstar_constant_uci <- (vir.par$c*predict.dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu_lci*predict.dat$Tw_constant_uci*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant_uci+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant_uci) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))


# Complete rstar predictions using  order-specific g0, Tw, and mu
# and default values for all other parameters
predict.dat$rstar_complete <- (vir.par$c*predict.dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu*vir.par$g*predict.dat$g0*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete-vir.par$Tv)))
predict.dat$rstar_complete_lci <- (vir.par$c*predict.dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu_uci)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu_lci*vir.par$g*predict.dat$g0*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete_lci-vir.par$Tv)))
predict.dat$rstar_complete_lci[is.na(predict.dat$rstar_complete_lci)] <- 0
predict.dat$rstar_complete_uci <- (vir.par$c*predict.dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu_lci)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu_uci*vir.par$g*predict.dat$g0*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete_uci-vir.par$Tv)))


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
  ylab(bquote("optimal virus growth rate in reservoir,"~r~"* (constant)"))


p10b <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  geom_errorbar(aes(x=order, ymin=rstar_complete_lci, ymax=rstar_complete_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(x=order, y=rstar_complete, fill=order, size=N_cumulative), pch=21, show.legend = F) +
  scale_color_manual(values=colz)+
  scale_fill_manual(values=colz, guide="none")+theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.5), "lines")) +
  ylab(bquote("optimal virus growth rate in reservoir,"~r~"* (complete)"))


p10 <- cowplot::plot_grid(p10a, p10b, ncol=1, nrow = 2, labels=c("A", "B"), rel_heights = c(1,1.3), label_x = -0.01)

#and save

ggsave(file = paste0(homewd, subwd, "/brief/Fig10.png"),
       plot = p10,
       units="mm",  
       width=50, 
       height=70, 
       scale=3, 
       dpi=300)



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
predict.dat <- arrange(predict.dat, desc(alpha_star_human_complete), desc(N_cumulative))
predict.dat$order <- factor(predict.dat$order, levels = unique(predict.dat$order))

# Visualize alphastar in humans

# Constant
p11a <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  geom_errorbar(aes(x=order, ymin=alpha_star_human_constant_lci, ymax=alpha_star_human_constant_uci, color=order),  width=0, linetype=3, show.legend = F) +
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
  

predict.dat <- arrange(predict.dat, desc(alpha_star_human_complete), desc(N_cumulative))
predict.dat$order <- factor(predict.dat$order, levels = unique(predict.dat$order))

# Complete
p11b <- ggplot(data=subset(predict.dat, !is.na(g0))) + 
  geom_errorbar(aes(x=order, ymin=alpha_star_human_complete_lci, ymax=alpha_star_human_complete_uci, color=order),  width=0, linetype=3, show.legend = F) +
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


#and save
ggsave(file = paste0(homewd, subwd, "/brief/Fig11.png"),
       plot = p11,
       units="mm",  
       width=50, 
       height=70, 
       scale=3, 
       dpi=300)

# And save the predict.dat for Fig S3

write.csv(predict.dat, file = paste0(homewd, subwd, "/predict_pars.csv"), row.names = F)

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

#and save
ggsave(file = paste0(homewd, subwd, "/brief/Fig12.png"),
       plot = p12,
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)


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
plot.dat <- plot.dat[complete.cases(plot.dat),]

shapez <- c("predicted from zoonoses"=25, "predicted from\nnested model"=24)


# And add images by order

order.dat <- read.csv(file=paste0(homewd,"/phylo-tree/Timetree_ReservoirMapping.csv"), header = T, stringsAsFactors = F)
head(order.dat)

# Rename
order.dat$species
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species[order.dat$species=="Homo_sapiens"] <-  "Gorilla_gorilla"
order.dat$species[order.dat$species=="Rattus_rattus"] <-  "Mus_musculus_domesticus"
order.dat$species[order.dat$species=="Antilocapra_americana"] <-  "Sus_scrofa"
order.dat$species[order.dat$species=="Phascolarctos_cinereus"] <-  "Macropus_rufus"
order.dat$species[order.dat$species=="Caenolestes_sangay"]<- "Caenolestes_convelatus"
order.dat$species[order.dat$species=="Sarcophilus_harrisii"] <- "Dasyurus_viverrinus"

# Load images - will take a moment
pic.df <- ggimage::phylopic_uid(order.dat$species) 

pic.df$order <- order.dat$order

unique(plot.dat$order)
unique(pic.df$order)
#and only plot those for which there is a comparison available
#plot.dat = subset(plot.dat, order=="Carnivora"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia" | order=="Primates"  | order=="Rodentia")
#pic.df = subset(pic.df, order=="Carnivora"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia"  | order=="Primates"  | order=="Rodentia")
pic.df = subset(pic.df, order=="Monotremata" | order=="Hyracoidea" | order == "Didelphimorphia" | order=="Proboscidea" | order=="Pilosa" | order== "Cingulata" | order == "Tubulidentata" | order=="Carnivora" | order=="Perissodactyla"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia"  | order=="Primates"  | order=="Rodentia" | order=="Eulipotyphla")


# This is a part of Fig. 3 in the main text
p13a <- ggplot(data=subset(plot.dat, tolerance!="complete"))  +  geom_hline(aes(yintercept=0), size=.2) +
  geom_errorbar(aes(x=order, ymin=alpha_lci, ymax=alpha_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(order, alpha, fill=order, size=N, shape=source)) + 
  theme_bw() +
  scale_color_manual(values=colz, guide="none") +
  scale_fill_manual(values=colz, guide="none") +
  scale_shape_manual(values=shapez, guide="none") +
  #facet_grid(source~., scales = "free_y") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16), 
        legend.direction = "horizontal", legend.position = c(.74,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,1.5,1,.2), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(constant tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.06)



p13 <- p13a + geom_text(x=16, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 


#and save
ggsave(file = paste0(homewd, subwd, "/brief/Fig13.png"),
       plot = p13,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)

# Then, plot complete tolerance assumptions:

# Reorder, ranked by the complete data
plot.dat$order <- factor(plot.dat$order, levels=unique(arrange(share.dat.complete, desc(alpha))$order))


# Plot 
p14a <- ggplot(data=subset(plot.dat, tolerance!="constant"))  +  geom_hline(aes(yintercept=0), size=.2) +
  geom_errorbar(aes(x=order, ymin=alpha_lci, ymax=alpha_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(order, alpha, fill=order, size=N, shape=source)) + 
  theme_bw() +
  scale_color_manual(values=colz, guide="none") +
  scale_fill_manual(values=colz, guide="none") +
  scale_shape_manual(values=shapez, guide="none") +
  #facet_grid(source~., scales = "free_y") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16), 
        legend.direction = "horizontal", legend.position = c(.74,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,1.5,1,.2), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(complete tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.06)



p14 <- p14a + geom_text(x=16, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 

ggsave(file = paste0(homewd, subwd, "/brief/Fig14.png"),
       plot = p14,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)

#this is also Figure S5 of the SI

ggsave(file = paste0(homewd, "/supp-figs/FigS5.png"),
       plot = p14,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)



######################################################
######################################################

# And now constant tolerance assumptions excluding rabies:


load(paste0(homewd, "source/SI.gam.dat.no.rabies.Guth.et.al.2021.Rdata"))

# Rescale alpha:
gam.dat$alpha <- scales::rescale(x =gam.dat$alpha, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
gam.dat$alpha_lci <- scales::rescale(x =gam.dat$alpha_lci, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 
gam.dat$alpha_uci <- scales::rescale(x =gam.dat$alpha_uci, from=c(min(gam.dat$alpha_lci), max(gam.dat$alpha_uci)), to =c(0,1)) 

#rank by descending alpha, then confidence
gam.dat <- arrange(gam.dat, desc(alpha), desc(Nobs))


# Merge our nested model predictions with those from Guth et al. 2022

# Make columns match
head(gam.dat)
names(gam.dat)[1] <- "order"
names(gam.dat)[3] <- "N"
names(predict.dat)[names(predict.dat)=="N_cumulative"] <- "N"
head(predict.dat)

# And merge the data again
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


# And add images by order
order.dat <- read.csv(file=paste0(homewd,"/phylo-tree/Timetree_ReservoirMapping.csv"), header = T, stringsAsFactors = F)
head(order.dat)

# Rename
order.dat$species
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species[order.dat$species=="Homo_sapiens"] <-  "Gorilla_gorilla"
order.dat$species[order.dat$species=="Rattus_rattus"] <-  "Mus_musculus_domesticus"
order.dat$species[order.dat$species=="Antilocapra_americana"] <-  "Sus_scrofa"
order.dat$species[order.dat$species=="Phascolarctos_cinereus"] <-  "Macropus_rufus"
order.dat$species[order.dat$species=="Caenolestes_sangay"]<- "Caenolestes_convelatus"
order.dat$species[order.dat$species=="Sarcophilus_harrisii"] <- "Dasyurus_viverrinus"

# Load images - will take a moment
pic.df <- ggimage::phylopic_uid(order.dat$species) 

pic.df$order <- order.dat$order

unique(plot.dat$order)
unique(pic.df$order)
#and only plot those for which there is a comparison available
#plot.dat = subset(plot.dat, order=="Carnivora"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia" | order=="Primates"  | order=="Rodentia")
#pic.df = subset(pic.df, order=="Carnivora"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia"  | order=="Primates"  | order=="Rodentia")
pic.df = subset(pic.df, order=="Monotremata" | order=="Hyracoidea" | order == "Didelphimorphia" | order=="Proboscidea" | order=="Pilosa" | order== "Cingulata" | order == "Tubulidentata" | order=="Carnivora" | order=="Perissodactyla"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia"  | order=="Primates"  | order=="Rodentia" | order=="Eulipotyphla")



shapez <- c("predicted from zoonoses"=25, "predicted from\nnested model"=24)

plot.dat$order <- factor(plot.dat$order, levels=unique(arrange(share.dat.constant, desc(alpha))$order))

# Plot Fig
p15a <- ggplot(data=subset(plot.dat, tolerance!="complete"))  +  geom_hline(aes(yintercept=0), size=.2) +
  geom_errorbar(aes(x=order, ymin=alpha_lci, ymax=alpha_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(order, alpha, fill=order, size=N, shape=source)) + 
  theme_bw() +
  scale_color_manual(values=colz, guide="none") +
  scale_fill_manual(values=colz, guide="none") +
  scale_shape_manual(values=shapez, guide="none") +
  #facet_grid(source~., scales = "free_y") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16), 
        legend.direction = "horizontal", legend.position = c(.74,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,0,1,1.3), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(constant tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.06)



p15 <- p15a + geom_text(x=16, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 



#and save
ggsave(file = paste0(homewd, subwd, "/brief/Fig15.png"),
       plot = p15,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)


# Then, plot complete tolerance assumptions:

# Reorder, ranked by the complete data
plot.dat$order <- factor(plot.dat$order, levels=unique(arrange(share.dat.complete, desc(alpha))$order))


# Plot 
p16a <- ggplot(data=subset(plot.dat, tolerance!="constant"))  +  geom_hline(aes(yintercept=0), size=.2) +
  geom_errorbar(aes(x=order, ymin=alpha_lci, ymax=alpha_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(order, alpha, fill=order, size=N, shape=source)) + 
  theme_bw() +
  scale_color_manual(values=colz, guide="none") +
  scale_fill_manual(values=colz, guide="none") +
  scale_shape_manual(values=shapez, guide="none") +
  #facet_grid(source~., scales = "free_y") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16), 
        legend.direction = "horizontal", legend.position = c(.74,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,1.5,1,1.3), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(complete tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.06)



p16 <- p16a + geom_text(x=16, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 

ggsave(file = paste0(homewd, subwd, "/brief/Fig16.png"),
       plot = p16,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)

# and pair them together to get Fig S6 of the SI

FigS6 <- cowplot::plot_grid(p15, p16, ncol = 2, nrow = 1, labels = c("A", "B"), label_size = 22)#, label_x = .12, label_y = .99)


ggsave(file = paste0(homewd, "/supp-figs/FigS6.png"),
       plot = FigS6,
       units="mm",  
       width=100, 
       height=60, 
       scale=3, 
       dpi=300)
