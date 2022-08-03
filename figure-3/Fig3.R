
# This script produces Figure 3 of the main text.

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
## Panel A

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

# Eventually will have wbc data too, so load it here to
# make your color bar

colz = scales::hue_pal()(length(unique(pan.dat$order))-1)
colz=c(colz, "red")

names(colz) <- c(unique(pan.dat$order)[unique(pan.dat$order)!="Chiroptera"], "Chiroptera")

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


# And plot with predictions:

pan.dat.sum <- cbind.data.frame(log10mass_g=seq(log10(min(pan.dat$mass_g)), max(log10(pan.dat$mass_g)), length=10), order = "Chiroptera")
pan.dat.sum$predict_lifespan <- 10^(predict.gam(m1, newdata = pan.dat.sum, exclude = "s(order)"))
pan.dat.sum$mass_g <- 10^(pan.dat.sum$log10mass_g)
pan.dat.sum = pan.dat.sum[pan.dat.sum$mass_g==min(pan.dat.sum$mass_g, na.rm=T) | pan.dat.sum$mass_g==max(pan.dat.sum$mass_g, na.rm=T),]

pAleg <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,   fill=order), size =3, pch=21) +  
  scale_fill_manual(values=colz) +
  geom_line(data = pan.dat.sum, aes(x=mass_g, y=predict_lifespan), size=1) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  xlab("mass (g)") + ylab("max lifespan (yrs)")
pAleg

# Grab leg
pleg <- cowplot::get_legend(pAleg)

pA <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,  fill=order), size =3, pch=21, show.legend = F) +  
  scale_fill_manual(values=colz) +
  geom_line(data = pan.dat.sum, aes(x=mass_g, y=predict_lifespan), size=1) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  xlab("mass (g)") + ylab("max lifespan (yrs)")
pA


######################################################
######################################################
## Panel B


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


pB <- ggplot(wbc.dat) + 
  geom_point(aes(x=BMR_W_g, y=neutro_conc, fill=order), size =3, pch=21, show.legend = F) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  geom_line(data = wbc.pred, aes(x=BMR_W_g , y=predict_neut), size=1) +
  scale_fill_manual(values=colz2) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=16),
        axis.text = element_text(size = 12),
        plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  xlab("mass-specific BMR (W/g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
pB




######################################################
######################################################
## Panel C



