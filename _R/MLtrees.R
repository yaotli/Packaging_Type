library(seqinr)
library(stringr)
library(ggtree)
library(ape)

source("~/Packaging_Type/_R/Function.R")
setwd("~/Desktop/data_souce/")

### raxml tree --------------------------------

## 234 ----------------
c234_ha_raxml <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/234_ha")
geo_region    <- read.csv("./tree_info/geo.csv")

colnames(geo_region)[1] <-  "taxa"
row.names(geo_region)   <-   NULL
geo_region$taxa         <-   paste0("'", geo_region$taxa, "'")

# host 

hu_i  <- 
  grep("[A-Z]+", 
       sapply( strsplit(geo_region$taxa, split = "_"), function(x) x[[2]] ) )[-c(7,28)]
pig_i <- 
  grep("swine",
       sapply( strsplit(geo_region$taxa, split = "_"), function(x) x[[2]] ), ignore.case = TRUE)

host        <- rep("avian", dim(geo_region)[1] )
host[hu_i]  <- "Human"
host[pig_i] <- "Swine"

host_i  <- rep(" ", dim(geo_region)[1])
host_i[ which(host == "Human") ] = "+"

geo_region  <- data.frame(geo_region, host, host_i, stringsAsFactors = FALSE)

# plot 

temp_t <- 
ggtree( c234_ha_raxml, ladderize = FALSE ) %<+% geo_region +
geom_tippoint( aes(fill = region), 
               color = "black", size = 2.5, stroke = 0.5, shape = 21 ) + 
scale_fill_manual( values = c("darkorange4", "darkorange", "gray60", "firebrick1", "chartreuse2") ) + 
theme( legend.position = "right" ) + 
geom_text( aes(label = host_i))  # , vjust = 0.35

b2344like <- findtaxa(type     = 1, 
                      tree     = c234_ha_raxml, 
                      targetid = c("EU195400_mallard_Huadong_lk_2005_H5N1_2005.496", 
                                   "KP233704_duck_Hunan_316_2005_H5N1_2005.301",
                                   "HM172100_duck_Jiangxi_80_2005_H5N1_2005.496",
                                   "DQ992838_crested_myna_Hong_Kong_540_2006_H5N1_2006.496",
                                   "DQ992790_duck_Hunan_324_2006_H5N1_2006.496"), 
                      target   = rep("dotted", 5) )
c234_t <- 
temp_t %<+% b2344like + aes(linetype = colorr) + geom_treescale(offset = 4)




## 232 ----------------
c232_ha_raxml  <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/232_ha")
geo_region_232 <- read.csv("./tree_info//geo_232.csv")

colnames(geo_region_232)[1] <-  "taxa"
row.names(geo_region_232)   <-   NULL
geo_region_232$taxa         <-   paste0("'", geo_region_232$taxa, "'")

# host 

hu_i  <- 
  grep("[A-Z]+", 
       sapply( strsplit(geo_region_232$taxa, split = "_"), function(x) x[[2]] ) )[ c(1,6,7,10,12) ]
pig_i <- 
  grep("swine",
       sapply( strsplit(geo_region_232$taxa, split = "_"), function(x) x[[2]] ), ignore.case = TRUE)

host        <- rep("avian", dim(geo_region_232)[1] )
host[hu_i]  <- "Human"
host[pig_i] <- "Swine"

host_i  <- rep(" ", dim(geo_region_232)[1])
host_i[ which(host == "Human") ] = "+"

geo_region_232  <- data.frame(geo_region_232, host, host_i, stringsAsFactors = FALSE)


# plot 
c232_t <- 
  ggtree(c232_ha_raxml, ladderize = FALSE) %<+% geo_region_232 +
  geom_tippoint(aes(fill = region), color = "black", size = 2.5, stroke = 0.5, shape = 21) + 
  scale_fill_manual(values = c("darkorange4", "darkorange", "gray60", "forestgreen", "firebrick1", "chartreuse2") ) + 
  theme(legend.position = "right") + geom_treescale(x = 0, y = 20, offset = 4) + 
  geom_text( aes(label = host_i) )




## 7 ----------------
c7_ha_raxml  <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/7_ha")
geo_region_7 <- read.csv("./tree_info//geo_7.csv")

colnames(geo_region_7)[1] <-  "taxa"
row.names(geo_region_7)   <-   NULL
geo_region_7$taxa         <-   paste0("'", geo_region_7$taxa, "'")

# host 

hu_i  <- 
  grep("[A-Z]+", 
       sapply( strsplit(geo_region_7$taxa, split = "_"), function(x) x[[2]] ) )[1]
pig_i <- 
  grep("swine",
       sapply( strsplit(geo_region_7$taxa, split = "_"), function(x) x[[2]] ), ignore.case = TRUE)

host        <- rep("avian", dim(geo_region_7)[1] )
host[hu_i]  <- "Human"
host[pig_i] <- "Swine"

host_i  <- rep(" ", dim(geo_region_7)[1])
host_i[ which(host == "Human") ] = "+"

geo_region_7  <- data.frame(geo_region_7, host, host_i, stringsAsFactors = FALSE)

# plot
c7_t <- 
  ggtree(c7_ha_raxml, ladderize = FALSE) %<+% geo_region_7 +
  geom_tippoint(aes(fill = region), color = "black", size = 2.5, stroke = 0.5, shape = 21) + 
  scale_fill_manual(values = c("darkorange4", "darkorange", "gray60", "forestgreen", "chartreuse2") ) + 
  theme(legend.position = "right") + 
  geom_treescale(x = 0, y = 20, offset = 1.5) +
  geom_text( aes(label = host_i), color = "white")

