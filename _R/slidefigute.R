require(ggtree)
require(ape)

setwd("/Volumes/EDGE2/LoVE/ReassortSubtype/data_souce/")
source("~/Packaging_Type/_R/Function.R")


## big 2344 tree ----------------

# 2014 --
trefile  <- read.nexus("Nov2017_hana/tree/rmd_pH5_7326/rmd_pH5_gsgd_p.tre")
targetid <- paste0( c("_2014.", "_2015.", "_2016.", "_2017."), collapse = "|" )
target   <- "#ff7f0e"


tree.d                         <- fortify(trefile)
tree.d[, ncol(tree.d) + 1]     <- gsub("'", "", tree.d$label)
colnames(tree.d)[ncol(tree.d)] <- "taxaid"

shapetaxa <- data.frame( node   = c( 1:length( tree.d$isTip ) ), shapee = NA)
shapetaxa$shapee[  grep(  targetid, tolower(tree.d$taxaid) ) ] <- target

ggtree(trefile, color = "gray", alpha = 0.6, size = 1.1) %<+% shapetaxa + 
  geom_tippoint( aes(color = I(shapee)), alpha = 0.7, size = 1) 
#+ geom_text(aes(label = node), size = 0.5)

# subtype

shapetaxa <- data.frame( shapetaxa, Nx = NA)

targetna <- c( "_H5N6_", "_H5N2_", "_H5N8_", "_H5N5_", "_H5N3_")
colna    <- c( "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

for( i in 1: length(targetna) )
{
    shapetaxa$Nx[ grep( targetna[i], tree.d$taxaid ) ] <- colna[i]
}

ls.c2344 <- read.table( "tree/rmd_pH5_7326/ls_c2344.txt", header = FALSE, stringsAsFactors = FALSE)
shapetaxa$Nx[ match( grep( "_H5N1_", ls.c2344[,1], ignore.case = TRUE, value = TRUE), tree.d$taxaid ) ] <- "black"


p = 
ggtree(trefile, color = "gray", alpha = 0.6, size = 1.1) %<+% shapetaxa + 
  geom_tippoint( aes(color = I(Nx)), alpha = 0.7, size = 2.5, stroke = 0)

# clade 2344
viewClade(p, node = 5008 )


## trees 234 ----------------

trefile <- "./Clade_allh5/234/c234_202.tre"
rawtre  <- read.tree(trefile)
tredata <- fortify( rawtre )
tredata <- data.frame(tredata, geo = NA, stringsAsFactors = FALSE)

tredata$geo[ 1:202 ]       <- geoID( tredata$label )[1:202]
tredata$geo[ c(192, 195) ] <- c("cnN", "cnE")
tredata$geo[ which(tredata$geo == "cnNW" ) ] <- "cnN"

table(tredata$geo)

col.geo <- c("#8c564b", "#2ca02c", "#7f7f7f", "#d62728", "#1f77b4", "#17becf")
#c("#8c564b", "#ff7f0e", "#7f7f7f", "#d62728", "#76ee00", "#17becf")
#76ee00
ml.234 <- 
ggtree( rawtre, right = TRUE, ladderize = FALSE) %<+% tredata +
  geom_tippoint( aes(fill = geo ), 
                 color = "black", size = 2.5, stroke = 0.5, shape = 21 ) + 
  scale_fill_manual( values = col.geo ) + 
  theme( legend.position = "right" )
flip(ml.234, 299, 339)



b_tre.234    <- "~/Desktop/b_geo/234/1012/c234_202_10101012-ann.tre"
rawbeast.234 <- read.beast( b_tre.234 )
tredata.234  <-  fortify( rawbeast.234 )

g <- 
ggtree( rawbeast.234, right = TRUE, size = 1, mrsd = "2011-11-18") + 
  aes( color = geo ) + 
  scale_color_manual( values = col.geo ) + 
  scale_fill_manual( values = col.geo ) +
  theme_tree2( axis.text.x = element_text( size = 16  ), legend.position = c(0,0), legend.justification = c(0,0)) + 
  scale_y_continuous( expand = c(0,5) ) + 
  geom_tippoint(aes(fill = geo), shape = 21, color = "black", stroke = 0.5 ) + 
  scale_x_continuous( breaks = seq(2002.5, 2012.5, by = 2), 
                      labels = seq(2002, 2012, by = 2) )  + 
  theme( axis.ticks  = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 18))


rectdf <- data.frame( xstart = seq( 2002, 2011, 2), 
                      xend   = seq( 2003, 2012, 2))
g + geom_rect(data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
            fill = "gray", alpha = 0.2, inherit.aes=FALSE)

ggsave("234.pdf", height = 5, width = 4)  
  

## trees 232 ----------------

  
trefile <- "./Clade_allh5/232/c232_149_e.tre"
rawtre  <- read.nexus(trefile)
tredata <- fortify( rawtre )
tredata <- data.frame(tredata, geo = NA, stringsAsFactors = FALSE)

tredata$geo[ 1:149 ]       <- geoID( tredata$label )[1:149]
tredata$geo[ c(3, 6, 9) ] <- c( "cnC", "cnNW", "cnS")
tredata$geo[ which(tredata$geo == "cnNW" ) ] <- "cnN"


col.geo <- c("#8c564b", "#ff7f0e", "#7f7f7f", "#d62728", "#76ee00", "#e377c2","#17becf")
#76ee00

ggtree( rawtre, right = TRUE) %<+% tredata +
  geom_tippoint( aes(fill = geo ), 
                 color = "black", size = 2.5, stroke = 0.5, shape = 21 ) + 
  scale_fill_manual( values = col.geo ) + 
  theme( legend.position = "right" )


b_tre.232    <- "~/Desktop/b_geo/232/c232_149_1007a_ann.trees"
rawbeast.232 <- read.beast( b_tre.232 )
tredata.232  <-  fortify( rawbeast.232 )

k <- 
  ggtree( rawbeast.232, right = TRUE, size = 0.6,mrsd = "2011-12-30") + 
  aes( color = geo ) + 
  scale_color_manual( values = col.geo ) + 
  scale_fill_manual( values = col.geo ) +
  theme_tree2( axis.text.x = element_text( size = 16  ), legend.position = c(0,0), legend.justification = c(0,0)) + 
  scale_y_continuous( expand = c(0,5) ) + 
  geom_tippoint(aes(fill = geo), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2014.5, by = 2), 
                      labels = seq(2002, 2014, by = 2) )  + 
  theme( axis.ticks  = element_blank() )


rectdf <- data.frame( xstart = seq( 2002, 2013, 2), 
                      xend   = seq( 2003, 2014, 2))
k + geom_rect(data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.2, inherit.aes=FALSE)

ggsave("232.pdf", height = 5, width = 4)  


### re-pH5 ML tree ----------------------------------

# 234
c234_ha_raxml <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/234_ha")
geo_region    <- read.csv("./tree_info/geo.csv")

colnames(geo_region)[1] <-  "taxa"
row.names(geo_region)   <-   NULL
geo_region$taxa         <-   paste0("'", geo_region$taxa, "'")

temp_t <- 
  ggtree( c234_ha_raxml, ladderize = FALSE, size = 1.1) %<+% geo_region +
  geom_tippoint( aes(fill = region), 
                 color = "black", size = 3, stroke = 0.5, shape = 21 ) + 
  scale_fill_manual( values = c("#8c564b", "#2ca02c", "#7f7f7f", "#d62728", "#1f77b4") ) 


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


# 232
c232_ha_raxml  <- read.tree("Clade/sampled_clade_tree_CNHK/sampled_besttree/232_ha")
geo_region_232 <- read.csv("./tree_info//geo_232.csv")

colnames(geo_region_232)[1] <-  "taxa"
row.names(geo_region_232)   <-   NULL
geo_region_232$taxa         <-   paste0("'", geo_region_232$taxa, "'")

geo_region_232$region[ which(geo_region_232$region == "NW" ) ] = "N"

c232_t <- 
  ggtree(c232_ha_raxml, ladderize = FALSE, size = 1.1 ) %<+% geo_region_232 +
  geom_tippoint(aes(fill = region), color = "black", size = 3, stroke = 0.5, shape = 21) + 
  scale_fill_manual( values = c("#8c564b", "#2ca02c", "#7f7f7f", "#d62728", "#1f77b4") ) +
  geom_treescale(x = 0, y = 20, offset = 4) + theme(legend.position = "right") 


### big c234 and c232 ML tree ------------------------------------

# c232
trefile.c232 <- "./c232/rmd/tree/ph5_c232_1007_rmdp_phy_phyml/ph5_c232_1007_rmdp_phy_phyml_e1201.tre"
tre.c232 <- read.nexus( trefile.c232 )

key_nCN.c232 <- 
paste0( grep( "China|Hong_Kong", as.data.frame( table(taxaInfo( trefile.c232, useTree = TRUE)[2]), stringsAsFactors = FALSE)$Var1, value = TRUE, invert = TRUE ),
        collapse = "|")

colordf.c232 <- 
findtaxa( type = 1, tree = tre.c232, targetid = key_nCN.c232, target = "red")

ggtree( tre.c232, ladderize = FALSE ) %<+% colordf.c232 + aes( color = I(colorr) ) +
  geom_treescale(width = 0.005, offset = -15 )

# c234
trefile.c234 <- "./c234/rmd/tree/ph5_c234_1699_rmdp_phy_phyml/ph5_c234_1699_rmdp_phy_phyml_e1201.tre"
tre.c234 <- read.nexus( trefile.c234 )

key_nCN.c234 <- 
  paste0( grep( "China|Hong_Kong", as.data.frame( table(taxaInfo( trefile.c234, useTree = TRUE)[2]), stringsAsFactors = FALSE)$Var1, value = TRUE, invert = TRUE ),
          collapse = "|")

colordf.c234 <- 
  findtaxa( type = 1, tree = tre.c234, targetid = key_nCN.c234, target = "red")

ggtree( tre.c234, ladderize = FALSE ) %<+% colordf.c234 + aes( color = I(colorr) ) + 
  geom_treescale(width = 0.005, offset = -25 )


### c232 + c234 vs other in China (line plot) --------------------------------

pri_tre_data    <- taxaInfo( useTree = TRUE, file = "Nov2017_hana/tree/rmd_pH5_7326/rmd_pH5_gsgd_e0108.tre")

df_pri_tre_data <- data.frame( time  = floor(pri_tre_data[[4]]), 
                               clade = pri_tre_data[[7]], 
                               geo   = pri_tre_data[[2]],
                               stringsAsFactors = FALSE)

df_pri_tre_data <- df_pri_tre_data[ which( df_pri_tre_data$geo == "China" | df_pri_tre_data$geo == "Hong_Kong"), ]

df_pri_tre_data$clade[ is.na(df_pri_tre_data$clade) ] <- "Other"
df_pri_tre_data$clade <- ifelse( df_pri_tre_data$clade == "Other", "Other", "c23")
  
  
df_pri_tre_data = as.data.frame( prop.table( table( df_pri_tre_data[, c(1,2)] ), margin = 1) )


ggplot( data = df_pri_tre_data ) + 
  geom_line( aes(x = time, y = Freq, color = clade, group = clade), size = 2) + 
  theme_bw() + xlab("") + ylab("% Sequences") +
  scale_x_discrete( breaks = seq(1996, 2016, by = 4), labels = seq(1996, 2016, by = 4) ) + 
  scale_color_manual( values = c( pyCol("red"), pyCol("gray") )  ) +
  theme( panel.grid.minor   = element_blank(), 
         panel.grid.major   = element_blank(), 
         axis.title       = element_text(size = 14),
         axis.ticks.x = element_blank(),
         axis.text    = element_text(size = 10),
         legend.position = "none"  )






