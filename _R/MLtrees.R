library(seqinr)
library(stringr)
library(ggtree)
library(ape)

source("~/Packaging_Type/_R/Function.R")
setwd("~/Desktop/data_souce/")

### raxml tree --------------------------------

## 234 ----------------
c234_ha_raxml <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/234_ha")
geo_region    <- read.csv("./geo.csv")

colnames(geo_region)[1] <-  "taxa"
row.names(geo_region)   <-   NULL
geo_region$taxa         <-   paste0("'", geo_region$taxa, "'")

temp_t <- 
ggtree(c234_ha_raxml, ladderize = FALSE) %<+% geo_region +
  geom_tippoint(aes(fill = region), color = "black", size = 2.5, stroke = 0.5, shape = 21) + 
  scale_fill_manual(values = c("darkorange4", "darkorange", "black", "firebrick1", "chartreuse2") ) + 
  theme(legend.position = "right") 

b2344like <- findtaxa(type     = 1, 
                      tree     = c234_ha_raxml, 
                      targetid = c("EU195400_mallard_Huadong_lk_2005_H5N1_2005.496", 
                                   "KP233704_duck_Hunan_316_2005_H5N1_2005.301",
                                   "HM172100_duck_Jiangxi_80_2005_H5N1_2005.496",
                                   "DQ992838_crested_myna_Hong_Kong_540_2006_H5N1_2006.496",
                                   "DQ992790_duck_Hunan_324_2006_H5N1_2006.496"), 
                      target   = rep("red", 5) )
c234_t <- 
temp_t %<+% b2344like + aes(color = I(colorr)) + geom_treescale(offset = 4)



## 232 ----------------
c232_ha_raxml  <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/232_ha")
geo_region_232 <- read.csv("./geo_232.csv")

colnames(geo_region_232)[1] <-  "taxa"
row.names(geo_region_232)   <-   NULL
geo_region_232$taxa         <-   paste0("'", geo_region_232$taxa, "'")

c232_t <- 
  ggtree(c232_ha_raxml, ladderize = FALSE) %<+% geo_region_232 +
  geom_tippoint(aes(fill = region), color = "black", size = 2.5, stroke = 0.5, shape = 21) + 
  scale_fill_manual(values = c("darkorange4", "darkorange", "black", "forestgreen", "firebrick1", "chartreuse2") ) + 
  theme(legend.position = "right") + geom_treescale(x = 0, y = 20, offset = 4)




## 7 ----------------
c7_ha_raxml  <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/7_ha")
geo_region_7 <- read.csv("./geo_7.csv")

colnames(geo_region_7)[1] <-  "taxa"
row.names(geo_region_7)   <-   NULL
geo_region_7$taxa         <-   paste0("'", geo_region_7$taxa, "'")

c7_t <- 
  ggtree(c7_ha_raxml, ladderize = FALSE) %<+% geo_region_7 +
  geom_tippoint(aes(fill = region), color = "black", size = 2.5, stroke = 0.5, shape = 21) + 
  scale_fill_manual(values = c("darkorange4", "darkorange", "black", "forestgreen", "chartreuse2") ) + 
  theme(legend.position = "right") + 
  geom_treescale(x = 0, y = 20, offset = 1.5)
