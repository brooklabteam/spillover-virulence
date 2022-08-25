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

# How many have lifespan too?
length(pan.dat$max_lifespan_yrs[!is.na(pan.dat$max_lifespan_yrs)]) #1055

# And how many have BMR?
length(pan.dat$BMR_W_g[!is.na(pan.dat$BMR_W_g)]) #629


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
plot.dat <- order.dat[[1]]
plot.dat$variable <- "max_lifespan"


# We know that constitutive immunity (baseline neutrophil count)
# scales positively with (a) body mass, # (b) negatively with BMR
# and mass-specific BMR, and positively with (c) higher host
# exposure/ promiscuity/ gregariousness

# We scale by (a) neutrophil data from ZIMs

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

wbc.dat$order <- as.factor(wbc.dat$order)

wbc.dat$ln10BMR_W_g <- log10(wbc.dat$BMR_W_g)
wbc.dat$ln10mass <- log10(wbc.dat$mass_g)

# Or try including BMR
m3 <- gam(ln10neutro~s(ln10mass, bs="tp") +
            s(BMR_W, bs="tp") +
            s(order, bs="re"),
          data=wbc.dat)
summary(m3) #63.5%; n=144. much stronger model fit

order.dat <- get_partial_effects(m3, var="order")
plot.dat2 <- order.dat[[1]]
head(plot.dat2)
plot.dat2$variable <- "neutrophil_conc"

plot.dat <- dplyr::select(plot.dat, y, se, ylower, yupper, IsSignificant, order, variable)
plot.dat2 <- dplyr::select(plot.dat2, y, se, ylower, yupper, IsSignificant, order, variable)

plot.dat$tag <- "A"
plot.dat2$tag <- "B"

plot.dat <- rbind(plot.dat, plot.dat2)

plot.dat$variable[plot.dat$variable=="max_lifespan"] <- "maximum lifespan"
plot.dat$variable[plot.dat$variable=="neutrophil_conc"] <- "neutrophil concentration"

fillz = c("No"="gray70", "Yes" = "skyblue3")

FigS3 <- ggplot(data=plot.dat, aes(order, y)) + 
  geom_crossbar(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                alpha=.4) +
  facet_grid(variable~., scales = "free_y", switch = "y") +
  #geom_point(aes(x=var, y=y, color=var), size=5) +
  #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
  scale_fill_manual(values = fillz, name = "Significant") +
  geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(size=12, angle = 90),
        plot.margin = unit(c(.1,.1,.5,.1), "cm"),
        strip.placement = "outside",
        legend.title = element_text(size=10),
        legend.direction = "horizontal",
        legend.text = element_text(size=8),
        strip.text = element_text(size = 14),
        legend.position = c(.84,.96),
        strip.background = element_blank()) +
  ylab(paste0("partial effect of order on response variable:")) + 
  geom_label(aes(x=1.5,y=.49, label=tag), label.size = 0, size=10)



#and save
ggsave(file = paste0(homewd, "supp-figs/FigS3.png"),
       plot = FigS3,
       units="mm",  
       width=60, 
       height=70, 
       scale=3, 
       dpi=300)
