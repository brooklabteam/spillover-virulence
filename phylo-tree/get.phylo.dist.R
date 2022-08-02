rm(list=ls())
library(ape)
library(ggplot2)

#get time to most recent common ancestory by order across several orders
#first load the sequences




homewd =  "/Users/carabrook/Developer/spillover-virulence/"
setwd(homewd)

tree <- ape::read.tree(file = paste0(homewd,"/phylo-tree/Timetree_ReservoirRepresentatives.nwk"))

plot(tree)
dist.matrix= cophenetic(tree)
head(dist.matrix)
hum.dat <- dist.matrix[colnames(dist.matrix)=="Homo_sapiens",]

phylo.dat <- cbind.data.frame(species=names(hum.dat), phylo_dist=hum.dat)
rownames(phylo.dat) <- c()

#and link to order and save
order.dat <- read.csv(file= paste0(homewd, "/phylo-tree/Timetree_ReservoirMapping.csv"), header = T, stringsAsFactors = F)
head(order.dat)
names(order.dat) <- c("order", "species")
order.dat$species <- sub(pattern = " ", replacement = "_", x=order.dat$species)
phylo.dat <- merge(phylo.dat, order.dat, by="species", all.x = T)

phylo.dat = subset(phylo.dat, order!="OUTGROUP")

phylo.dat <- dplyr::select(phylo.dat, order, phylo_dist)
#note that this is divergence by MILLIONS OF YEARS

save(phylo.dat, file = "phylo-tree/phylo.dat.final.Rdata")

