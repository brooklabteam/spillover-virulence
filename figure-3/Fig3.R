
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

names(colz) <- c(sort(unique(pan.dat$order)[unique(pan.dat$order)!="Chiroptera"]), "Chiroptera")

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

pA <- ggplot(data=pan.dat) + 
  geom_point(aes(x=mass_g, y=max_lifespan_yrs,  fill=order), size =3, pch=21)+#, show.legend = F) +  
  scale_fill_manual(values=colz) +
  geom_line(data = pan.dat.sum, aes(x=mass_g, y=predict_lifespan), size=1) +
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18),
        legend.position = c(.77,.18),
        legend.title = element_blank(),
        legend.key.size = unit(.3, "lines"),
        legend.background = element_rect(color="black"),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5,.65,.1,.15), "lines")) +
  xlab("mass (g)") + ylab("max lifespan (yrs)")
pA


######################################################
######################################################
## Panel B


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
length(wbc.dat$mass_g[!is.na(wbc.dat$BMR_W)]) #144 with BMR


# Slim your colors down to only those in this dataset
colz2 = sort(colz[unique(wbc.dat$order)])

wbc.dat$order <- as.factor(wbc.dat$order)

# plot with predictions
# faux gam first
wbc.dat$ln10BMR_W_g <- log10(wbc.dat$mass_g)/wbc.dat$BMR_W
wbc.dat$ln10BMR_W_g <- log10(wbc.dat$BMR_W_g)


m2a <- gam(ln10neutro~s(ln10BMR_W_g, bs="tp") +
             s(order, bs="re"),
           data=wbc.dat)
summary(m2a) #51.8%; n =144

wbc.pred = cbind.data.frame(ln10BMR_W_g=seq(min(wbc.dat$ln10BMR_W_g, na.rm=T),  max(wbc.dat$ln10BMR_W_g, na.rm=T), length=2),  order = "Chiroptera")
wbc.pred$predict_neut <- 10^(predict.gam(m2a, newdata = wbc.pred, exclude = "s(order)"))
wbc.pred$BMR_W_g <- 10^(wbc.pred$ln10BMR_W_g)


pB <- ggplot(wbc.dat) + 
  geom_point(aes(x=BMR_W_g, y=neutro_conc, fill=order), size =3, pch=21, show.legend = F) +  
  theme_bw() + scale_y_log10() + scale_x_log10(labels=scales::comma) +
  geom_line(data = wbc.pred, aes(x=BMR_W_g , y=predict_neut), size=1) +
  scale_fill_manual(values=colz) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=18),
        axis.text = element_text(size = 14),
        plot.margin = unit(c(.5,2,.1,1.6), "lines")) +
  xlab("mass-specific BMR (W/g)") + ylab(bquote("neutrophil concentrations ("~10^9~"cells/L)"))
pB


######################################################
######################################################
## Panel C (phylogeny)

# Load the timetree from the phylo-tree subfolder
tree <- ape::read.tree(file = paste0(homewd,"phylo-tree/Timetree_ReservoirRepresentatives.nwk"))

# Rename for the icons we want to use
tree$tip.label
tree$tip.label <- sub(pattern = " ", replacement = "_", x=tree$tip.label)
tree$tip.label[tree$tip.label=="Homo_sapiens"] <-  "Pan"
tree$tip.label[tree$tip.label=="Rattus_rattus"] <-  "Mus_musculus_domesticus"
tree$tip.label[tree$tip.label=="Antilocapra_americana"] <-  "Sus_scrofa"
tree$tip.label[tree$tip.label=="Phascolarctos_cinereus"] <-  "Macropus_rufus"
tree$tip.label[tree$tip.label=="Caenolestes_sangay"]<- "Caenolestes_convelatus"
tree$tip.label[tree$tip.label=="Sarcophilus_harrisii"] <- "Dasyurus_viverrinus"
tree$tip.label[tree$tip.label=="Oryctolagus_cuniculus"] <- "Ochotona_princeps"
tree$tip.label[tree$tip.label=="Tupaia_glis"] <- "Dermoptera"


# Call the phylopics
tree$tip.label[tree$tip.label=="Rhyncholestes_raphanurus"] <- "Paucituberculata"
tree$tip.label[tree$tip.label=="Suncus_murinus"] <- "Atopogale_cubana"
tree$tip.label[tree$tip.label=="Manis_javanica"] <- "Smutsia_gigantea"
tree$tip.label[tree$tip.label=="Acerodon_celebensis"] <- "Pteropus_medius"
tree$tip.label[tree$tip.label=="Procavia_capensis"] <- "Hyracoidea"
tree$tip.label[tree$tip.label=="Dugong_dugon"] <- "Dugongidae"
tree$tip.label[tree$tip.label=="Elephantulus_myurus"] <- "Macroscelididae"
tree$tip.label[tree$tip.label=="Equus_asinus"] <- "Equus_ferus_caballus"
d <- ggimage::phylopic_uid(tree$tip.label)


tree.dat <- cbind.data.frame(species=tree$tip.label)



# Load the order meta-data and phylogenetic distance and combine with tree
order.dat <- read.csv(file=paste0(homewd, "/phylo-tree/Timetree_ReservoirMapping.csv"), header = T, stringsAsFactors = F)
head(order.dat)
order.dat$species
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
#order.dat$species[order.dat$species=="Homo_sapiens"] <-  "Macaca_mulatta"
order.dat$species[order.dat$species=="Rattus_rattus"] <-  "Mus_musculus_domesticus"
order.dat$species[order.dat$species=="Antilocapra_americana"] <-  "Sus_scrofa"
order.dat$species[order.dat$species=="Phascolarctos_cinereus"] <-  "Macropus_rufus"
order.dat$species[order.dat$species=="Caenolestes_sangay"]<- "Caenolestes_convelatus"
order.dat$species[order.dat$species=="Sarcophilus_harrisii"] <- "Dasyurus_viverrinus"
order.dat$species[order.dat$species=="Oryctolagus_cuniculus"] <- "Ochotona_princeps"
order.dat$species[order.dat$species=="Tupaia_glis"] <- "Dermoptera"

tree.dat$order <- tree.dat$phylo_dist <- NA

for(i in 1:length(order.dat$order)){
  tree.dat$order[tree.dat$species==order.dat$species[i]] <- order.dat$order[i]
  tree.dat$phylo_dist[tree.dat$species==order.dat$species[i]] <- order.dat$phylo_dist[i]
}
tree.dat$color= NA
for(i in 1:length(colz)){
  tree.dat$color[tree.dat$order==names(colz)[i]] <- colz[i]
}


d$phylo_dist = tree.dat$phylo_dist
d$name <- tree.dat$order
d$order <- tree.dat$order
d$order <- as.factor(d$order)
d$phylo_dist <- as.numeric(as.character(d$phylo_dist))
d$color <- tree.dat$color
tree$tip.label <- tree.dat$order

d$phylo_lab = round(d$phylo_dist, 0)
names(d)[names(d)=="name"] <- "tip_label"

#tree$tip.label[is.na(tree$tip.label)] <- "Paucituberculata"
#d$color[is.na(d$tip_label)] <- colz[names(colz)=="Paucituberculata"]
d$order <- as.character(d$order)
#d$order[is.na(d$tip_label)] <- "Paucituberculata"
d$order <- as.factor(d$order)
#d$tip_label[is.na(d$tip_label)] <- "Paucituberculata"

pC1 <- ggtree(tree, size=1)   %<+% d + 
  aes(color=phylo_dist) +
  scale_x_reverse() +
  scale_color_viridis_c(name="phylogenetic distance\nfrom primates (Myr)", direction=-1, 
                        guide = guide_legend(direction = "vertical",title.position = "top")) +
  geom_tiplab(aes(label=label), color="black",offset = -65, size=4) + 
  # scale_color_gradient(low="black", high="red",  name=bquote(eta)) +
  ggnewscale::new_scale_color()+  theme_bw() +
  geom_tiplab(aes(image=uid, color=order),  geom="phylopic",offset = -10, size=.04) +
  scale_color_manual(values=colz, guide="none") + 
  theme(legend.position = c(.85,.8), 
        legend.direction = "vertical", legend.title = element_text(size=12),
        legend.text = element_text(size=12), legend.key.size = unit(c(.8), "cm"),
        plot.margin = unit(c(.3,.3,9.8,3), "lines"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) #+
#geom_text(aes(x=branch, label=phylo_lab), size=5, color="black", vjust=-0.4) 

# 
# 
# # Now flip branches into the correct order
node1 <- MRCA(tree, which(tree$tip.label=="Primates"), which(tree$tip.label=="Lagomorpha"))
node2 <- MRCA(tree, which(tree$tip.label=="Carnivora"), which(tree$tip.label=="Eulipotyphla"))

pC <- flip(pC1, node1, node2) 
# pB3 <- flip(pB2, 7, 8) 
# pB4 <- flip(pB3, 2, 12) 
# pB5 <- flip(pB4, 6, 13) 



######################################################
######################################################
## Panel D - compare alpha with literature


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

# This is a part of Fig. 3 in the main text
pDa <- ggplot(data=subset(plot.dat, tolerance!="complete"))  +  geom_hline(aes(yintercept=0), size=.2) +
  geom_errorbar(aes(x=order, ymin=alpha_lci, ymax=alpha_uci, color=order),  width=0, linetype=3, show.legend = F) +
  geom_point(aes(order, alpha, fill=order, size=N, shape=source)) + 
  theme_bw() +
  scale_color_manual(values=colz, guide="none") +
  scale_fill_manual(values=colz, guide="none") +
  scale_shape_manual(values=shapez, guide="none") +
  #facet_grid(source~., scales = "free_y") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=18), 
        legend.direction = "horizontal", legend.position = c(.74,.95),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(size=14, vjust=.1, hjust=-.2, angle=90),
        plot.margin = unit(c(.3,1,.8,1), "cm")) + 
  ylab(bquote("relative spillover virulence,"~alpha[S])) + 
  scale_y_continuous(breaks=c(-1,-.5, 0, .5, 1), labels=c(1,.5, 0, .5, 1)) +
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") + 
  geom_phylopic(data=pic.df, aes(x=order, y = -1.3, image=uid, color=order), size=.05)

pD <- pDa + geom_text(x=20, y=0, label="      From nested model       From zoonotic literature   ", angle=270, nudge_y = 2, size=6) + 
  coord_cartesian(ylim=c(-1.1,1.1), clip = "off") 



### And combine

Fig3top <- cowplot::plot_grid(pA, pB, ncol = 2, nrow = 1, labels = c("A", "B"), label_size = 22, rel_widths = c(1,1.1))
Fig3bottom <- cowplot::plot_grid(pC, pD, ncol = 2, nrow = 1, labels = c("C", "D"), label_size = 22, rel_widths = c(1,1.1))

Fig3 <- cowplot::plot_grid(Fig3top, Fig3bottom, ncol = 1, nrow = 2, rel_heights = c(1,1.25))

ggsave(file = paste0(homewd,"/main-figs/Fig3.png"),
       plot = Fig3,
       units="mm",  
       width=120, 
       height=120, 
       scale=3, 
       dpi=300)

#and pdf for submission

ggsave(file = paste0(homewd,"/main-figs/Fig3.pdf"),
       plot = Fig3,
       units="mm",  
       width=120, 
       height=120, 
       scale=3, 
       dpi=300)
