rm(list=ls())

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(deSolve)
library(RColorBrewer)
library(scales)
library(mgcv)
library(stringr)
library(lmerTest)

#setwd

homewd = "/Users/carabrook/Developer/spillover-virulence/"
subwd = "source"
setwd(paste0(homewd, subwd))


# This script is pulled from Guth et al. 2021 and estimates order-specific 
# virulence for zoonoses derived from disparate mammalian reservoirs

############################
############################
## First, model CFR by order as in Guth et al. 2021

## Load data
dat <- read.csv(file="stringent_data_guth2021.csv", header = T, stringsAsFactors = F)

head(dat)

unique(dat$hOrder)
dat <- subset(dat, hOrder!="AVES")# & hOrder!="PILOSA")
## convert all categorical predictors to factors
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- str_to_title(dat$hOrder)
dat$hOrder <- as.factor(dat$hOrder)
dat$Tr.primary <- as.factor(dat$Tr.primary)
dat$vFamily <- as.factor(dat$vFamily)
dat$phylo_dist <- as.numeric(dat$phylo_dist)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)

# Convert spill_type to binary
dat$spill_type[dat$spill_type=="bridged"] <- 1
dat$spill_type[dat$spill_type=="direct"] <- 0
dat$spill_type <- as.factor(dat$spill_type)

# Remove viruses with only one human case recorded in history
dat_rare <- dat[dat$tot_cases==1, ]
dat_suf <- dat[dat$tot_cases>1, ]

# Final set of model variables
dat.mort <- dat_suf %>% distinct(SppName_ICTV_MSL2018b, CFR_avg, hOrder, spill_type, duration_symptoms_avg_days, .keep_all = TRUE)
length(unique(dat.mort$SppName_ICTV_MSL2018b))

# How many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 3 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

# Convert CFRs to proportions to model as beta distribution
dat.mort$CFR <- dat.mort$CFR_avg/100

# Transformation for extreme 0 and 1 values
dat.mort$CFR <- (dat.mort$CFR*(86-1) + 0.5)/86

# Global model w/ global CFR estimates
gam_mort <- gam(CFR ~
                  s(vFamily, bs = 're'  ) +
                  #s(phylo_dist, k=5, bs="tp") +
                  #s(ReservoirNspecies, k=7, bs="tp") +
                  s(hOrder, bs = 're'  ) +
                  #s(VirusSppPublicationCount, k=7, bs="tp") +
                  s(spill_type, bs = 're' )  +
                  s(IsVectorBorne, bs = 're'),
                data = dat.mort,
                select=FALSE,
                family = betar(link = "logit"))
summary(gam_mort)

unique(dat.mort$hOrder)#8 orders - Pilosa dropped
# Build empty dataset to house predictions
# Include all predictor variables in final global model from Guth et al 
gam.dat <- ddply(dat, .(hOrder), summarise, 
                 ReservoirPublicationCount = median(ReservoirPublicationCount), 
                 Nobs=length(hOrder))

# Remove birds and Pilosa
gam.dat <- subset(gam.dat, hOrder!="Aves" & hOrder!="Pilosa" )
gam.dat$hOrder <- as.factor(gam.dat$hOrder)
gam.dat$spill_type <- 0
gam.dat$spill_type <- as.factor(gam.dat$spill_type)
gam.dat$IsVectorBorne <- 0
gam.dat$IsVectorBorne <- as.factor(gam.dat$IsVectorBorne)

# Fill in a place holder for virus family but then exclude it in the predictions
gam.dat$vFamily <- "Poxviridae"
gam.dat$CFR = NA

#Predict CF by order, excluding viral family
gam.dat$CFR <- predict.gam(gam_mort, newdata=gam.dat, exclude = "s(vFamily)", type = "response")

#Get confidence intervals on those estimates
gam.dat$CFR_lci <- gam.dat$CFR - (1.96*predict(gam_mort, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$CFR_uci <-gam.dat$CFR + (1.96*predict(gam_mort, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$CFR_lci[gam.dat$CFR_lci<0] <- 0
gam.dat$CFR[gam.dat$CFR<0] <- 0
gam.dat$hOrder <- str_to_title(gam.dat$hOrder)

gam.dat <- arrange(gam.dat, desc(CFR))

gam.dat$hOrder <- factor(gam.dat$hOrder, levels=unique(gam.dat$hOrder))

gam.dat$CFR <- 100*gam.dat$CFR
gam.dat$CFR_lci <- 100*gam.dat$CFR_lci
gam.dat$CFR_uci <- 100*gam.dat$CFR_uci


colz= c("Afrosoricida" = "#F8766D",  "Carnivora"="#D89000",  "Cetartiodactyla"="#A3A500", "Chiroptera"="red", "Dasyuromorphia"="#39B600", 
        "Diprotodontia"="#00BF7D", "Eulipotyphla"="#00BFC4", "Peramelemorphia"="#E76BF3", "Perissodactyla"="#9590FF",  "Primates"="#00B0F6" ,  "Rodentia"="#FF62BC")


# And plot
p1 <- ggplot(data=gam.dat) + 
  geom_errorbar(aes(x=hOrder, ymin=CFR_lci, ymax=CFR_uci, color=hOrder),  width=0, linetype=3) +
  geom_point(aes(hOrder, CFR, color=hOrder, size=Nobs)) + 
  theme_bw() +
  #scale_color_manual(values=colz) +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12, angle = 90), legend.title = element_blank()) + 
  ylab(bquote("predicted CFR by order")) #+ coord_cartesian(ylim=c(0,100))

print(p1)


############################
############################
## Now, we would prefer to measure alpha, 
# so need to convert CFR to alpha using duration of infection

#alpha = (CFR*(1-Di*mu))/Di
#assume human mu = (1/(79*365))

head(dat.mort)
dat.mort$alpha = (dat.mort$CFR/dat.mort$duration_symptoms_avg_days)

#alpha is in units of day^-1

# And plot by order
p2 <- ggplot(dat.mort) + 
      geom_boxplot(aes(x=hOrder, y=duration_symptoms_avg_days, color=hOrder)) + 
      scale_y_log10() +scale_color_manual(values=colz) +
      theme(axis.text.x = element_text(angle=45), 
        axis.title.x = element_blank())
print(p2)

p3 <- ggplot(dat.mort) + 
    geom_boxplot(aes(x=hOrder, y=alpha, color=hOrder)) + 
    scale_y_log10() +scale_color_manual(values=colz) +
      theme(axis.text.x = element_text(angle=45), 
            axis.title.x = element_blank())
print(p3)

p4 <- ggplot(dat.mort) + 
  geom_boxplot(aes(x=vFamily, y=CFR, color=vFamily)) + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle=45), 
        axis.title.x = element_blank())
print(p4)


# Now, make the same GAM for alpha
gam_alpha <- gam(alpha ~
                  s(vFamily, bs = 're'  ) +
                  #s(phylo_dist, k=5, bs="tp") +
                  #s(ReservoirNspecies, k=7, bs="tp") +
                  s(hOrder, bs = 're'  ) +
                  s(VirusSppPublicationCount, k=7, bs="tp") +
                  s(spill_type, bs = 're' )  +
                  s(IsVectorBorne, bs = 're'),
                data = dat.mort,
                select=FALSE,
                family = betar(link = "logit"))
summary(gam_alpha)

# And visualize the partial effects
source(paste0(homewd,"source/mollentze-streicker-2020-functions.R"))

# Get partial effects of order on alpha
m1out <- get_partial_effects(fit= gam_alpha, var="hOrder", seWithMean = T)

plot.partial <- function(df, var, ylab1){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  
  
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  #p2 <- ggplot(data=df2, aes(var,  Residual)) +
  #     geom_boxplot(aes(var~Residual))
  
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_crossbar(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                  alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, angle = 45),
          plot.margin = unit(c(.1,.1,.5,.1), "cm"))+
    ylab(ylab1) 
  
  #print(p1)
  
  return(p1)
  
}

# Plot 
plot.partial(m1out, var="hOrder", ylab1 = "partial effect of order on alpha") #primates have significantly high Tw, eulipotyphla significantly low. bats are not sig but right on the edge



# Now, estimate alpha by order
# Build empty dataset to house predictions
# Include all predictor variables in final global model from Guth et al 
gam.dat <- ddply(dat, .(hOrder), summarise, 
                 VirusSppPublicationCount = median(VirusSppPublicationCount), 
                 Nobs=length(hOrder))
gam.dat <- subset(gam.dat, hOrder!="Aves" & hOrder!="Pilosa" )
gam.dat$hOrder <- as.factor(gam.dat$hOrder)
gam.dat$spill_type <- 0
gam.dat$spill_type <- as.factor(gam.dat$spill_type)
gam.dat$IsVectorBorne <- 0
gam.dat$IsVectorBorne <- as.factor(gam.dat$IsVectorBorne)
# Fill in a place holder for virus family but then exclude it in the predictions
gam.dat$vFamily <- "Poxviridae"
gam.dat$alpha = NA

gam.dat$alpha <- predict.gam(gam_alpha, newdata=gam.dat, exclude = "s(vFamily)", type = "response")

gam.dat$alpha_lci <- gam.dat$alpha - (1.96*predict(gam_alpha, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$alpha_uci <-gam.dat$alpha + (1.96*predict(gam_alpha, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$alpha_lci[gam.dat$alpha_lci<0] <- 0
gam.dat$alpha[gam.dat$alpha<0] <- 0
gam.dat$hOrder <- str_to_title(gam.dat$hOrder)

gam.dat <- arrange(gam.dat, desc(alpha))

gam.dat$hOrder <- factor(gam.dat$hOrder, levels=unique(gam.dat$hOrder))

colz= c("Afrosoricida" = "#F8766D",  "Carnivora"="#D89000",  "Cetartiodactyla"="#A3A500", "Chiroptera"="red", "Dasyuromorphia"="#39B600", 
        "Diprotodontia"="#00BF7D", "Eulipotyphla"="#00BFC4", "Peramelemorphia"="#E76BF3", "Perissodactyla"="#9590FF",  "Primates"="#00B0F6" ,  "Rodentia"="#FF62BC")


# And plot it...
p5 <- ggplot(data=gam.dat) + 
  geom_errorbar(aes(x=hOrder, ymin=alpha_lci, ymax=alpha_uci, color=hOrder),  width=0, linetype=3) +
  geom_point(aes(hOrder, alpha, color=hOrder, size=Nobs)) + 
  theme_bw() +
  scale_color_manual(values=colz) +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12, angle = 90), legend.title = element_blank()) + 
  ylab(bquote("predicted"~alpha[S]~"by order ("~days^-1~")")) #+ coord_cartesian(ylim=c(0,100))


p5

save(gam.dat, file=paste0(homewd,"source/gam.dat.Guth.et.al.2021.Rdata"))




# And redo above, excluding rabies
rm(list=ls())

homewd = "/Users/caraebrook/Documents/R/R_repositories/virulence-evolution/spillover-virulence/"
subwd = "source"
setwd(paste0(homewd, subwd))



## Load data
dat <- read.csv(file="stringent_data_guth2021.csv", header = T, stringsAsFactors = F)

head(dat)

unique(dat$hOrder)
dat <- subset(dat, hOrder!="AVES" & hOrder!="PILOSA")
dat <- subset(dat, SppName_ICTV_MSL2018b!="Rabies lyssavirus")
## Convert all categorical predictors to factors
dat$SppName_ICTV_MSL2018b <- as.factor(dat$SppName_ICTV_MSL2018b)
dat$hOrder <- str_to_title(dat$hOrder)
dat$hOrder <- as.factor(dat$hOrder)
dat$Tr.primary <- as.factor(dat$Tr.primary)
dat$vFamily <- as.factor(dat$vFamily)
dat$phylo_dist <- as.numeric(dat$phylo_dist)
dat$IsVectorBorne <- as.factor(dat$IsVectorBorne)

# Convert spill_type to binary
dat$spill_type[dat$spill_type=="bridged"] <- 1
dat$spill_type[dat$spill_type=="direct"] <- 0
dat$spill_type <- as.factor(dat$spill_type)

# Remove viruses with only one human case recorded in history
dat_rare <- dat[dat$tot_cases==1, ]
dat_suf <- dat[dat$tot_cases>1, ]

# Final set of model variables
dat.mort <- dat_suf %>% distinct(SppName_ICTV_MSL2018b, CFR_avg, hOrder, spill_type, duration_symptoms_avg, .keep_all = TRUE)
length(unique(dat.mort$SppName_ICTV_MSL2018b))



# How many viruses have multiple entries?
dup <- dat.mort[duplicated(dat.mort$SppName_ICTV_MSL2018b)==TRUE, ] #only 2 viruses
dup_dat <- dat.mort %>%
  filter(SppName_ICTV_MSL2018b %in% dup$SppName_ICTV_MSL2018b)

# Convert CFRs to proportions to model as beta distribution
dat.mort$CFR <- dat.mort$CFR_avg/100

# Transformation for extreme 0 and 1 values
dat.mort$CFR <- (dat.mort$CFR*(86-1) + 0.5)/86

# We would prefer to measure alpha, so need to convert CFR to alpha using duration infection
#alpha = (CFR*(1-Di*mu))/Di
#assume human mu = (1/(79*365))

head(dat.mort)
dat.mort$alpha = (dat.mort$CFR/dat.mort$duration_symptoms_avg)


colz= c("Afrosoricida" = "#F8766D",  "Carnivora"="#D89000",  "Cetartiodactyla"="#A3A500", "Chiroptera"="red", "Dasyuromorphia"="#39B600", 
        "Diprotodontia"="#00BF7D", "Eulipotyphla"="#00BFC4", "Peramelemorphia"="#E76BF3", "Perissodactyla"="#9590FF",  "Primates"="#00B0F6" ,  "Rodentia"="#FF62BC")


#and plot by order
p2 <- ggplot(dat.mort) + 
  geom_boxplot(aes(x=hOrder, y=duration_symptoms_avg, color=hOrder)) + 
  scale_y_log10() +scale_color_manual(values=colz) +
  theme(axis.text.x = element_text(angle=45), 
        axis.title.x = element_blank())
print(p2)

p3 <- ggplot(dat.mort) + 
  geom_boxplot(aes(x=hOrder, y=alpha, color=hOrder)) + 
  scale_y_log10() +scale_color_manual(values=colz) +
  theme(axis.text.x = element_text(angle=45), 
        axis.title.x = element_blank())
print(p3)

p4 <- ggplot(dat.mort) + 
  geom_boxplot(aes(x=vFamily, y=CFR, color=vFamily)) + 
  scale_y_log10() +
  theme(axis.text.x = element_text(angle=45), 
        axis.title.x = element_blank())
print(p4)


# And gam of alpha
gam_alpha <- gam(alpha ~
                   s(vFamily, bs = 're'  ) +
                   #s(phylo_dist, k=5, bs="tp") +
                   #s(ReservoirNspecies, k=7, bs="tp") +
                   s(hOrder, bs = 're'  ) +
                   s(VirusSppPublicationCount, k=7, bs="tp") +
                   s(spill_type, bs = 're' )  +
                   s(IsVectorBorne, bs = 're'),
                 data = dat.mort,
                 select=FALSE,
                 family = betar(link = "logit"))
summary(gam_alpha)

# And show the partial effects
source(paste0(homewd,"source/mollentze-streicker-2020-functions.R"))

# Make tolerance (Tw) a scaled version of the order effect in this relationship
m1out <- get_partial_effects(fit= gam_alpha, var="hOrder", seWithMean = T)

plot.partial <- function(df, var, ylab1){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  
  
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  #p2 <- ggplot(data=df2, aes(var,  Residual)) +
  #     geom_boxplot(aes(var~Residual))
  
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_crossbar(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                  alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, angle = 45),
          plot.margin = unit(c(.1,.1,.5,.1), "cm"))+
    ylab(ylab1) 
  
  #print(p1)
  
  return(p1)
  
}
plot.partial(m1out, var="hOrder", ylab1 = "partial effect of order on alpha") #primates have significantly high Tw, eulipotyphla significantly low. bats are not sig but right on the edge


# Build empty dataset to house predictions
# Include all predictor variables in final global model from Guth et al 
gam.dat <- ddply(dat, .(hOrder), summarise, 
                 VirusSppPublicationCount = median(VirusSppPublicationCount), 
                 Nobs=length(hOrder))
gam.dat <- subset(gam.dat, hOrder!="Aves" & hOrder!="Pilosa" )
gam.dat$hOrder <- as.factor(gam.dat$hOrder)
gam.dat$spill_type <- 0
gam.dat$spill_type <- as.factor(gam.dat$spill_type)
gam.dat$IsVectorBorne <- 0
gam.dat$IsVectorBorne <- as.factor(gam.dat$IsVectorBorne)

# Fill in a place holder for virus family but then exclude it in the predictions
gam.dat$vFamily <- "Poxviridae"
gam.dat$alpha = NA

gam.dat$alpha <- predict.gam(gam_alpha, newdata=gam.dat, exclude = "s(vFamily)", type = "response")

gam.dat$alpha_lci <- gam.dat$alpha - (1.96*predict(gam_alpha, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$alpha_uci <-gam.dat$alpha + (1.96*predict(gam_alpha, newdata=gam.dat, type="response", se.fit=TRUE)$se)
gam.dat$alpha_lci[gam.dat$alpha_lci<0] <- 0
gam.dat$alpha[gam.dat$alpha<0] <- 0
gam.dat$hOrder <- str_to_title(gam.dat$hOrder)

gam.dat <- arrange(gam.dat, desc(alpha))

gam.dat$hOrder <- factor(gam.dat$hOrder, levels=unique(gam.dat$hOrder))

colz= c("Afrosoricida" = "#F8766D",  "Carnivora"="#D89000",  "Cetartiodactyla"="#A3A500", "Chiroptera"="red", "Dasyuromorphia"="#39B600", 
        "Diprotodontia"="#00BF7D", "Eulipotyphla"="#00BFC4", "Peramelemorphia"="#E76BF3", "Perissodactyla"="#9590FF",  "Primates"="#00B0F6" ,  "Rodentia"="#FF62BC")


# And plot
p5 <- ggplot(data=gam.dat) + 
  geom_errorbar(aes(x=hOrder, ymin=alpha_lci, ymax=alpha_uci, color=hOrder),  width=0, linetype=3) +
  geom_point(aes(hOrder, alpha, color=hOrder, size=Nobs)) + 
  theme_bw() +
  scale_color_manual(values=colz) +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12, angle = 90), legend.title = element_blank()) + 
  ylab(bquote("predicted"~alpha[S]~"by order")) #+ coord_cartesian(ylim=c(0,100))


p5

# Save the version excluding rabies too
save(gam.dat, file=paste0(homewd, "source/SI.gam.dat.no.rabies.Guth.et.al.2021.Rdata"))


