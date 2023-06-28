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



#now


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

#save this color bar
save(colz, file = paste0(homewd, "/figure-3/color-bar.Rdata"))

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


# Now, run a glm  to determine the partial effect of host order
# (phylogeny) on maximum longevity, and predict lifespan 
# rates by order:


pan.dat$log10mass_g <- log10(pan.dat$mass_g)
pan.dat$order <- as.factor(pan.dat$order)
pan.dat$log10_max_lifespan_yrs <- log10(pan.dat$max_lifespan_yrs)

library(lme4)
library(sjPlot)
m1 <- lmer(log10_max_lifespan_yrs~ log10mass_g + (1|order), data = pan.dat) #mass longevity linear interaction with separate intercepts per order
plot_model(m1, type = "re")


summary(m1) 

order.dat <- cbind.data.frame(order = plot_model(m1, type = "re")$data$term, estimate= plot_model(m1, type = "re")$data$estimate, lci= plot_model(m1, type = "re")$data$conf.low, uci= plot_model(m1, type = "re")$data$conf.high)
order.dat$sig <- "no"
order.dat$sig[order.dat$lci>0 & order.dat$uci>0] <- "pos"
order.dat$sig[order.dat$lci<0 & order.dat$uci<0] <- "neg"
colz2 = c('no' = "grey", 'pos'="red3", 'neg'='blue4')


mass.dat <- cbind.data.frame(mass= 10^(plot_model(m1, type = "pred")$log10mass$data$x), lifespan=10^(plot_model(m1, type = "pred")$log10mass$data$predicted), lifespan_uci=10^(plot_model(m1, type = "pred")$log10mass$data$conf.high), lifespan_lci=10^(plot_model(m1, type = "pred")$log10mass$data$conf.low))

# And, plot the partial effects of order:

newdata.df <- cbind.data.frame(log10mass_g=3.229679, order="Chiroptera")

p2a <- ggplot(data=order.dat) + geom_point(aes(x=order, y=estimate, color=sig), size=3, show.legend = F) + 
        geom_linerange(aes(x=order, ymin=lci, ymax=uci, color=sig), show.legend = F) + 
        scale_color_manual(values=colz2) +
        ylab("order y-intercept on mass : lifespan relationship") +
        geom_hline(aes(yintercept = 0), linetype=2) +
        coord_flip() + theme_bw() + theme(panel.grid = element_blank(),
                                         axis.text.x = element_text(size = 13),
                                         axis.text.y = element_text(size = 10),
                                         plot.margin = unit(c(.3,.3,.1,.5), "lines"), 
                                         axis.title.y = element_blank(),
                                         axis.title.x = element_text(size=14))

mass.dat$mass_kg <- mass.dat$mass/1000

p2b <-ggplot(data=mass.dat) + geom_line(aes(x=mass_kg, y=lifespan), size=1, color="red3") + 
   geom_ribbon(aes(x=mass_kg, ymin=lifespan_lci, ymax=lifespan_uci), alpha=.3, fill="red3") + 
   xlab("mass (kg)") + scale_y_log10() +
   scale_x_log10(breaks = c(.1, 100, 1*10^5), labels=c(".1", "100","100000")) +
   ylab("max lifespan (years)") + 
   geom_hline(aes(yintercept = 0), linetype=2) +
   theme_bw() + 
   theme(panel.grid = element_blank(), 
         axis.text = element_text(size = 13),
         plot.margin = unit(c(.1,.3,.3,3.7), "lines"), 
         axis.title = element_text(size=14))
 
 
 p2 <- cowplot::plot_grid(p2a, p2b, ncol=1, nrow = 2, labels=c("A", "B"), label_x = .002, label_size = 22)
# 
 ggsave(file = paste0(homewd, subwd, "/brief/Fig2.png"),
        plot = p2,
        units="mm",  
        width=30, 
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

ggsave(file = paste0(homewd, subwd, "/brief/Fig1w_line.png"),
       plot = p1b,
       units="mm",  
       width=60, 
       height=30, 
       scale=3, 
       dpi=300)


# Now, to parameterize mu, the annual mortality rate by order for
# our within-host model, we can refit a glm simply using order as a direct
# predictor of lifespan excluding the effects of body mass, then take the inverse (but be sure to
# express in timesteps of days to make comparable to the viral dynamics):


# Remake lm with no involvement of mass - 
#first express annual mortality rate in 
pan.dat$mortality_rate_per_year <- 1/(pan.dat$max_lifespan_yrs)


m1b <- lm(mortality_rate_per_year~ order, data = pan.dat) 
summary(m1b)

pan.test = pan.dat
pan.test$order <- factor(pan.test$order, levels=c("Diprotodontia", unique(as.character(pan.test$order[as.character(pan.test$order)!="Diprotodontia"]))))

m1b2 <- lm(mortality_rate_per_year~ order, data = pan.test) 
summary(m1b2)

tables2A <- as.data.frame(summary(m1b2)$coefficients)
tables2A$order <- rownames(tables2A) 
tables2A$order[tables2A$order=="(Intercept)"] <- "Intercept (Diprotodontia)"
tables2A$order <- sub(pattern = "order", replacement = "", x=tables2A$order, fixed = T)
names(tables2A)<- c("est", "se", "tval", "pval", "order")
tables2A$lci <- tables2A$est-1.96*tables2A$se
tables2A$uci <- tables2A$est+1.96*tables2A$se
tables2A$est <- round(tables2A$est,2)
tables2A$lci <- round(tables2A$lci,2)
tables2A$uci <- round(tables2A$uci,2)
tables2A$lci_uci <- paste0(tables2A$lci," - ",tables2A$uci)

tables2A <- dplyr::select(tables2A, order, est, lci_uci, tval, pval)
tables2A$tval <- round(tables2A$tval,2)
tables2A$pval <- round(tables2A$pval,2)

tables2A$pval[tables2A$pval==0.00] <- "<0.001***"

write.csv(tables2A, file = paste0(homewd, "figure-3/tables2A.csv"), row.names = F)

plot_model(m1b, type="est") #this is using Afrosoricida as a reference, so almost everything looks long-lived by comparison
plot_model(m1b, type="pred")



# n.sum = ddply(pan.dat[!is.na(pan.dat$mortality_rate_per_year),], .(order), summarise, Nmort = length(mortality_rate_per_year))
# 
# pan.dat <- merge(pan.dat, n.sum, by="order")
# 
# 
# m1b <- gam(log10_max_lifespan_yrs~ s(order, bs="re"), data = pan.dat) 
# summary(m1b)

pan.dat$order <- factor(pan.dat$order, levels = c(sort(unique(as.character(pan.dat$order)))))

predict.dat <- cbind.data.frame(order = unique(pan.dat$order))#, log10mass_g=1)

#now express mortality rate in timestep of day
predict.dat$mu <- predict(m1b, newdata = predict.dat, type="response")/365
predict.dat$mu_se <- predict(m1b, newdata = predict.dat, type="response", se.fit = T)$se.fit


#get standard error:
predict.dat$mu_lci <- (predict(m1b, newdata = predict.dat, type="response") -1.96*predict.dat$mu_se)/365
predict.dat$mu_uci <- (predict(m1b, newdata = predict.dat, type="response") +1.96*predict.dat$mu_se)/365


# 

predict.dat$mu_lci[predict.dat$mu_lci<0] <- 0
predict.dat$mu_uci[predict.dat$mu_uci>1] <- 1

# Also compute the number of data entries per order used to determine this:
mu.sum <- ddply(pan.dat, .(order), summarise, N_mu = length(binomial))

predict.dat <- merge(predict.dat, mu.sum, by="order", all.x = T)
#predict.dat_b <- merge(predict.dat_b, mu.sum, by="order", all.x = T)
predict.dat <- dplyr::select(predict.dat, -(mu_se))


#what is the middle one? ~diprotodontia
arrange(predict.dat, mu)



# and plot your predictions by order for mass

# first, get your null
# y.int = 1/((10^(predict.gam(m1, 
#                     newdata = cbind.data.frame(order="Primates", log10mass_g = unique(order.dat[[1]]$log10mass_g)), 
#                     exclude = "s(order)", type = "response")))*365)

#set this y-int at the median of all the orders

#y.int = (max(predict.dat$mu) - min(predict.dat$mu))/2 + min(predict.dat$mu)
#y.int = predict.dat$mu[predict.dat$order=="Diprotodontia"]

y.int=median(predict.dat$mu)

pan.dat.test = pan.dat
pan.dat.test$order <- factor(as.character(pan.dat.test$order), levels=c("Diprotodontia", unique(as.character(pan.dat.test$order[as.character(pan.dat.test$order)!="Diprotodontia"]))))
m1c <- lm(mortality_rate_per_year~ order, data = pan.dat.test) 
out.sum = as.data.frame(summary(m1c)$coefficients)
names(out.sum) <- c("est", "se", "t_val", "p_val")

subset(out.sum, p_val<.001)


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

#or try it as the residual

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
        legend.position = c(.7,.91),
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
# Afrosoricida, Cetartiodactyla, Dasyuromorphia, Didelphimorphia,
# Eulipotyphla, Notoryctemorphia, Peramelemorphia, Rodentia



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
wbc.dat$ln10BMR_W_g <- log10(wbc.dat$BMR_W_g)

wbc.dat$order <- as.factor(wbc.dat$order)


m2a <- lmer(ln10neutro~ln10BMR_W_g + (1|order), data=subset(wbc.dat, !is.na(ln10BMR_W_g)))
plot_model(m2a, type="re")



wbc.pred = cbind.data.frame(ln10BMR_W_g=seq(min(wbc.dat$ln10BMR_W_g, na.rm=T),  max(wbc.dat$ln10BMR_W_g, na.rm=T), length=2),  order = "Chiroptera")
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

ggsave(file = paste0(homewd, subwd, "/brief/Fig6_w_line.png"),
       plot = p6b,
       units="mm",  
       width=50, 
       height=30, 
       scale=3, 
       dpi=300)


 wbc.dat$order <- as.factor(wbc.dat$order)
 wbc.dat$ln10mass <- log10(wbc.dat$mass_g)



# now look at mass and BMR independently 
m3 <- lmer(ln10neutro~ln10mass+BMR_W + (1|order), data=subset(wbc.dat, !is.na(ln10BMR_W_g)))
summary(m3)

plot_model(m3, type="re")
plot_model(m3, type="pred")$ln10mass
plot_model(m3, type="pred")$BMR_W


summary(m3) 
AIC(m2a, m3)

order.dat <- cbind.data.frame(order = plot_model(m3, type = "re")$data$term, estimate= plot_model(m3, type = "re")$data$estimate, lci= plot_model(m3, type = "re")$data$conf.low, uci= plot_model(m3, type = "re")$data$conf.high)
order.dat$sig <- "no"
order.dat$sig[order.dat$lci>0 & order.dat$uci>0] <- "pos"
order.dat$sig[order.dat$lci<0 & order.dat$uci<0] <- "neg"

colz3 = c('no' = "grey", 'pos'="red3", 'neg'='blue4')



mass.dat <- cbind.data.frame(mass= 10^(plot_model(m3, type = "pred")$ln10mass$data$x), neutro=10^(plot_model(m3, type = "pred")$ln10mass$data$predicted), neutro_uci=10^(plot_model(m3, type = "pred")$ln10mass$data$conf.high), neutro_lci=10^(plot_model(m3, type = "pred")$ln10mass$data$conf.low))
BMR.dat <- cbind.data.frame(BMR_W= (plot_model(m3, type = "pred")$BMR_W$data$x), neutro=10^(plot_model(m3, type = "pred")$BMR_W$data$predicted), neutro_uci=10^(plot_model(m3, type = "pred")$BMR_W$data$conf.high), neutro_lci=10^(plot_model(m3, type = "pred")$BMR_W$data$conf.low))



p7a <- ggplot(data=order.dat) + geom_point(aes(x=order, y=estimate, color=sig), size=3, show.legend = F) + 
  geom_linerange(aes(x=order, ymin=lci, ymax=uci, color=sig), show.legend = F) + 
  scale_color_manual(values=colz3) +
  ylab("order y-intercept on mass/BMR/neutrophil relationship") +
  geom_hline(aes(yintercept = 0), linetype=2) +
  coord_flip() + theme_bw() + theme(panel.grid = element_blank(),
                                    axis.text.x = element_text(size = 13),
                                    axis.text.y = element_text(size = 10),
                                    plot.margin = unit(c(.3,.3,.3,.5), "lines"), 
                                    axis.title.y = element_blank(),
                                    axis.title.x = element_text(size=14))

mass.dat$mass_kg <- mass.dat$mass/1000

p7b <-ggplot(data=mass.dat) + geom_line(aes(x=mass_kg, y=neutro), size=1, color="red3") + 
  geom_ribbon(aes(x=mass_kg, ymin=neutro_lci, ymax=neutro_uci), alpha=.3, fill="red3") + 
  xlab("mass (kg)") + scale_y_log10() +
  scale_x_log10(breaks = c(.1, 10, 1000), labels=c(".1", "10","1000")) +
  ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))+
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 13),
        plot.margin = unit(c(.2,.3,.3,3.7), "lines"), 
        axis.title = element_text(size=14))

p7c <-ggplot(data=BMR.dat) + geom_line(aes(x=BMR_W, y=neutro), size=1, color="blue4") + 
  geom_ribbon(aes(x=BMR_W, ymin=neutro_lci, ymax=neutro_uci), alpha=.3, fill="blue4") + 
  xlab("BMR (W)") + scale_y_log10() +
  #scale_x_log10(breaks = c(.1, 10, 1000), labels=c(".1", "10","1000")) +
  ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))+
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 13),
        plot.margin = unit(c(.2,.3,.3,3.7), "lines"), 
        axis.title = element_text(size=14))

# #and plot together

 p7 <- cowplot::plot_grid(p7a, p7b, p7c, ncol=1, nrow = 3, labels=c("A", "B", "C"), label_x = 0.01)
# 

 ggsave(file = paste0(homewd, subwd, "/brief/Fig7.png"),
        plot = p7,
        units="mm",  
        width=50, 
        height=100, 
        scale=3, 
        dpi=300)

# Significant negative associations with:
#  Cetartiodactyal, Dasyuromorphia, Diprotodontia, Scandentia


#Significant positive associations with:
# Monotremata, Primates



#and finally, calculate g0 as the order effects from this model

tmp.dat = cbind.data.frame(order= plot_model(m3, type="re")$data$term, estimate=plot_model(m3, type="re")$data$estimate, conf.low=plot_model(m3, type="re")$data$conf.low, conf.high=plot_model(m3, type="re")$data$conf.high) 

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
        legend.position = c(.65,.92),
        legend.direction = "horizontal",
        legend.title = element_blank()) +
  geom_errorbar(aes(x=order, ymin=g0_lci, ymax=g0_uci, color=order), width=0, show.legend = F) +
  scale_color_manual(values=colz2) +
  scale_fill_manual(values=colz2, guide="none") +
  ylab(bquote("magnitude constitutive immunity, "~g[0]))

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

vir.par$Tv = 1.005
predict.dat$rstar_constant <- (vir.par$c*predict.dat$g0/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu*predict.dat$Tw_constant*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))
predict.dat$rstar_constant_lci <- (vir.par$c*predict.dat$g0_lci/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0_lci*predict.dat$mu_lci*predict.dat$Tw_constant_lci*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant_lci+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant_lci) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))
predict.dat$rstar_constant_uci <- (vir.par$c*predict.dat$g0_uci/vir.par$m) + (sqrt((vir.par$m^2)*(vir.par$c^2)*vir.par$g*predict.dat$g0_uci*predict.dat$mu_uci*predict.dat$Tw_constant_uci*vir.par$Tv*(vir.par$v*predict.dat$Tw_constant_uci+vir.par$g*vir.par$w*vir.par$Tv)))/((vir.par$v*(vir.par$m^2)*predict.dat$Tw_constant_uci) +(vir.par$g*vir.par$w*(vir.par$m^2)*vir.par$Tv))


# Complete rstar predictions using  order-specific g0, Tw, and mu
# and default values for all other parameters
vir.par$Tv=.005
predict.dat$rstar_complete <- (vir.par$c*predict.dat$g0/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0*predict.dat$mu)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu*vir.par$g*predict.dat$g0*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete-vir.par$Tv)))
predict.dat$rstar_complete_lci <- (vir.par$c*predict.dat$g0_lci/vir.par$m) + ((vir.par$c^2)*vir.par$g*predict.dat$g0_lci*predict.dat$mu_lci)/(sqrt((vir.par$m^2)*(vir.par$c^2)*predict.dat$mu_lci*vir.par$g*predict.dat$g0_lci*(vir.par$g*vir.par$w+vir.par$v-vir.par$g*predict.dat$Tw_complete_lci-vir.par$Tv)))
#predict.dat$rstar_complete_lci[is.na(predict.dat$rstar_complete_lci)] <- .00000001
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

#and save

ggsave(file = paste0(homewd, subwd, "/brief/Fig10.png"),
       plot = p10,
       units="mm",  
       width=50, 
       height=70, 
       scale=3, 
       dpi=300)

# This is also Fig. S5 of the main text:

ggsave(file = paste0(homewd, "/supp-figs/FigS5.png"),
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
# predict.dat$Vs_max_constant <- predict.dat$rstar_constant/(vir.par$g*vir.par$c) - predict.dat$rstar_constant/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_constant*vir.par$g)
# predict.dat$Vs_max_constant_lci <- predict.dat$rstar_constant_lci/(vir.par$g*vir.par$c) - predict.dat$rstar_constant_lci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_constant_lci*vir.par$g)
# predict.dat$Vs_max_constant_uci <- predict.dat$rstar_constant_uci/(vir.par$g*vir.par$c) - predict.dat$rstar_constant_uci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_constant_uci*vir.par$g)
# 
# predict.dat$Vs_max_complete <- predict.dat$rstar_complete/(vir.par$g*vir.par$c) - predict.dat$rstar_complete/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_complete*vir.par$g)
# predict.dat$Vs_max_complete_lci <- predict.dat$rstar_complete_lci/(vir.par$g*vir.par$c) - predict.dat$rstar_complete_lci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_complete_lci*vir.par$g)
# predict.dat$Vs_max_complete_uci <- predict.dat$rstar_complete_uci/(vir.par$g*vir.par$c) - predict.dat$rstar_complete_uci/(2*vir.par$g*vir.par$c) + 1 - 1/(vir.par$g) + vir.par$c/(2*predict.dat$rstar_complete_uci*vir.par$g)

#rep as VS_avg
predict.dat$Vs_avg_constant <- (predict.dat$rstar_constant*(sqrt((vir.par$c^2)+2*predict.dat$rstar_constant*vir.par$c*(vir.par$g-1)+(predict.dat$rstar_constant^2))) +(predict.dat$rstar_constant^2) + 4*(predict.dat$rstar_constant)*(vir.par$c)*(vir.par$g-1) + 2*(vir.par$c^2))/(6*vir.par$c*vir.par$g*(predict.dat$rstar_constant))
predict.dat$Vs_avg_constant_lci <- (predict.dat$rstar_constant_lci*(sqrt((vir.par$c^2)+2*predict.dat$rstar_constant_lci*vir.par$c*(vir.par$g-1)+(predict.dat$rstar_constant_lci^2))) +(predict.dat$rstar_constant_lci^2) + 4*(predict.dat$rstar_constant_lci)*(vir.par$c)*(vir.par$g-1) + 2*(vir.par$c^2))/(6*vir.par$c*vir.par$g*(predict.dat$rstar_constant_lci))
predict.dat$Vs_avg_constant_uci <- (predict.dat$rstar_constant_uci*(sqrt((vir.par$c^2)+2*predict.dat$rstar_constant_uci*vir.par$c*(vir.par$g-1)+(predict.dat$rstar_constant_uci^2))) +(predict.dat$rstar_constant_uci^2) + 4*(predict.dat$rstar_constant_uci)*(vir.par$c)*(vir.par$g-1) + 2*(vir.par$c^2))/(6*vir.par$c*vir.par$g*(predict.dat$rstar_constant_uci))

predict.dat$Vs_avg_complete <- (predict.dat$rstar_complete*(sqrt((vir.par$c^2)+2*predict.dat$rstar_complete*vir.par$c*(vir.par$g-1)+(predict.dat$rstar_complete^2))) +(predict.dat$rstar_complete^2) + 4*(predict.dat$rstar_complete)*(vir.par$c)*(vir.par$g-1) + 2*(vir.par$c^2))/(6*vir.par$c*vir.par$g*(predict.dat$rstar_complete))
predict.dat$Vs_avg_complete_lci <- (predict.dat$rstar_complete_lci*(sqrt((vir.par$c^2)+2*predict.dat$rstar_complete_lci*vir.par$c*(vir.par$g-1)+(predict.dat$rstar_complete_lci^2))) +(predict.dat$rstar_complete_lci^2) + 4*(predict.dat$rstar_complete_lci)*(vir.par$c)*(vir.par$g-1) + 2*(vir.par$c^2))/(6*vir.par$c*vir.par$g*(predict.dat$rstar_complete_lci))
predict.dat$Vs_avg_complete_uci <- (predict.dat$rstar_complete_uci*(sqrt((vir.par$c^2)+2*predict.dat$rstar_complete_uci*vir.par$c*(vir.par$g-1)+(predict.dat$rstar_complete_uci^2))) +(predict.dat$rstar_complete_uci^2) + 4*(predict.dat$rstar_complete_uci)*(vir.par$c)*(vir.par$g-1) + 2*(vir.par$c^2))/(6*vir.par$c*vir.par$g*(predict.dat$rstar_complete_uci))



# Next, take that viral load and calculate spillover virulence, here for the constant tolerance assumption:
vir.par$Tw = 1
predict.dat$alpha_star_human_constant <- ((predict.dat$rstar_constant*vir.par$v)/predict.dat$Tv_human_constant  + (vir.par$g*vir.par$w*predict.dat$rstar_constant)/vir.par$Tw)*predict.dat$Vs_avg_constant
predict.dat$alpha_star_human_constant_lci <- ((predict.dat$rstar_constant_lci*vir.par$v)/predict.dat$Tv_human_constant  + (vir.par$g*vir.par$w*predict.dat$rstar_constant_lci)/vir.par$Tw)*predict.dat$Vs_avg_constant_lci
predict.dat$alpha_star_human_constant_uci <- ((predict.dat$rstar_constant_uci*vir.par$v)/predict.dat$Tv_human_constant  + (vir.par$g*vir.par$w*predict.dat$rstar_constant_uci)/vir.par$Tw)*predict.dat$Vs_avg_constant_uci


# And here for complete tolerance:
vir.par$Tw = 0
predict.dat$alpha_star_human_complete <- (predict.dat$rstar_complete*(vir.par$v-predict.dat$Tv_human_complete) + predict.dat$rstar_complete*vir.par$g*(vir.par$w-vir.par$Tw))*predict.dat$Vs_avg_complete
predict.dat$alpha_star_human_complete_lci <- (predict.dat$rstar_complete_lci*(vir.par$v-predict.dat$Tv_human_complete) + predict.dat$rstar_complete_lci*vir.par$g*(vir.par$w-vir.par$Tw))*predict.dat$Vs_avg_complete_lci
predict.dat$alpha_star_human_complete_uci <- (predict.dat$rstar_complete_uci*(vir.par$v-predict.dat$Tv_human_complete) + predict.dat$rstar_complete_uci*vir.par$g*(vir.par$w-vir.par$Tw))*predict.dat$Vs_avg_complete_uci

# Now, we rank by virulence, then confidence, and plot...
predict.dat <- arrange(predict.dat, desc(alpha_star_human_constant), desc(N_cumulative))
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
  

predict.dat <- arrange(predict.dat, desc(alpha_star_human_constant), desc(N_cumulative))
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
plot.dat$alpha_lci[is.na(plot.dat$alpha_lci) & !is.na(plot.dat$alpha)] <- 0
plot.dat$alpha_uci[is.na(plot.dat$alpha_uci) & !is.na(plot.dat$alpha)] <- 1

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
order.dat$species[order.dat$species=="Oryctolagus_cuniculus"] <- "Ochotona_princeps"
order.dat$species[order.dat$species=="Tupaia_glis"] <- "Dermoptera"



# Load images - will take a moment
pic.df <- ggimage::phylopic_uid(order.dat$species) 

pic.df$order <- order.dat$order

unique(plot.dat$order)
unique(pic.df$order)

setdiff(unique(pic.df$order), unique(plot.dat$order)) #these are those that don't overlap
sort(intersect(unique(pic.df$order), unique(plot.dat$order))) #these are those that do overlap
#and only plot those for which there are predictions 
pic.df = subset(pic.df, order=="Afrosoricida" | order == "Carnivora" | order=="Cetartiodactyla" | order=="Chiroptera" |order== "Cingulata" | 
                order == "Dasyuromorphia" | order == "Didelphimorphia" | order=="Diprotodontia"  | order=="Eulipotyphla" | order=="Hyracoidea" |
                order=="Monotremata" | order == "Peramelemorphia" | order == "Perissodactyla" | order=="Pilosa" | order=="Primates"  |
                order == "Proboscidea" | order=="Rodentia" | order == "Scandentia" | order == "Tubulidentata" )
                  
library(ggimage)

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
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.05)



p13 <- p13a + geom_text(x=21, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
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
share.dat.complete$alpha_lci[share.dat.complete$alpha_lci=="NaN"] <- 0
share.dat.complete <- arrange(share.dat.complete, desc(alpha))
plot.dat$order <- factor(plot.dat$order, levels=unique(as.character(share.dat.complete$order)))

plot.comp = subset(plot.dat, tolerance!="constant")
plot.comp <- arrange(plot.comp, order)
pic.df <- arrange(pic.df, order= unique(as.character(plot.dat$order)))

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
        legend.direction = "horizontal", legend.position = c(.7,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,1.5,1,.2), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(complete tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.05)



p14 <- p14a + geom_text(x=21, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 

ggsave(file = paste0(homewd, subwd, "/brief/Fig14.png"),
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
plot.dat$alpha_lci[!is.na(plot.dat$alpha) & is.na(plot.dat$alpha_lci)] <- 0
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
order.dat$species[order.dat$species=="Oryctolagus_cuniculus"] <- "Ochotona_princeps"
order.dat$species[order.dat$species=="Tupaia_glis"] <- "Dermoptera"

# Load images - will take a moment
pic.df <- ggimage::phylopic_uid(order.dat$species) 

pic.df$order <- order.dat$order

unique(plot.dat$order)
unique(pic.df$order)
#and only plot those for which there is a comparison available
#plot.dat = subset(plot.dat, order=="Carnivora"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia" | order=="Primates"  | order=="Rodentia")
#pic.df = subset(pic.df, order=="Carnivora"  | order=="Cetartiodactyla" | order=="Chiroptera" |  order=="Diprotodontia"  | order=="Primates"  | order=="Rodentia")
pic.df = subset(pic.df, order=="Afrosoricida" | order == "Carnivora" | order=="Cetartiodactyla" | order=="Chiroptera" |order== "Cingulata" | 
                  order == "Dasyuromorphia" | order == "Didelphimorphia" | order=="Diprotodontia"  | order=="Eulipotyphla" | order=="Hyracoidea" |
                  order=="Monotremata" | order == "Peramelemorphia" | order == "Perissodactyla" | order=="Pilosa" | order=="Primates"  |
                  order == "Proboscidea" | order=="Rodentia" | order == "Scandentia" | order == "Tubulidentata" )



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
        legend.direction = "horizontal", legend.position = c(.7,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,1.5,1,1), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(constant tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.05)



p15 <- p15a + geom_text(x=21, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
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
        legend.direction = "horizontal", legend.position = c(.65,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.2,1.5,1,1.3), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S]~"(complete tolerance)")) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.05)



p16 <- p16a + geom_text(x=21, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 

ggsave(file = paste0(homewd, subwd, "/brief/Fig16.png"),
       plot = p16,
       units="mm",  
       width=50, 
       height=60, 
       scale=3, 
       dpi=300)


