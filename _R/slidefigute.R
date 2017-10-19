require(ggtree)
require(ape)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")


## big 2344 tree ----------------

# 2014 --
trefile  <- read.tree("./Tree/allh5_GsGDlike_GsGDonly_e1011")
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

p = 
ggtree(trefile, color = "gray", alpha = 0.6, size = 1.1) %<+% shapetaxa + 
  geom_tippoint( aes(color = I(Nx)), alpha = 0.7, size = 4, stoke = 0)

# clade 2344
viewClade(p, node = 6347 )


## trees 234 ----------------

trefile <- "./Clade_allh5/234/c234_202.tre"
rawtre  <- read.tree(trefile)
tredata <- fortify( rawtre )
tredata <- data.frame(tredata, geo = NA, stringsAsFactors = FALSE)

tredata$geo[ 1:202 ]       <- geoID( tredata$label )[1:202]
tredata$geo[ c(192, 195) ] <- c("cnN", "cnE")
tredata$geo[ which(tredata$geo == "cnNW" ) ] <- "cnN"

table(tredata$geo)

col.geo <- c("#8c564b", "#ff7f0e", "#7f7f7f", "#d62728", "#76ee00", "#17becf")
#76ee00
ml.234 <- 
ggtree( rawtre, right = TRUE, ladderize = FALSE) %<+% tredata +
  geom_tippoint( aes(fill = geo ), 
                 color = "black", size = 2.5, stroke = 0.5, shape = 21 ) + 
  scale_fill_manual( values = col.geo ) + 
  theme( legend.position = "right" )
flip(ml.234, 299, 339)



b_tre.234    <- "~/Desktop/b_geo/234/c234_202_1011-ann.trees"
rawbeast.234 <- read.beast( b_tre.234 )
tredata.234  <-  fortify( rawbeast.234 )

g <- 
ggtree( rawbeast.234, right = TRUE, size = 0.6,mrsd = "2011-11-18") + 
  aes( color = geo ) + 
  scale_color_manual( values = col.geo ) + 
  scale_fill_manual( values = col.geo ) +
  theme_tree2( axis.text.x = element_text( size = 16  ), legend.position = c(0,0), legend.justification = c(0,0)) + 
  scale_y_continuous( expand = c(0,5) ) + 
  geom_tippoint(aes(fill = geo), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2012.5, by = 2), 
                      labels = seq(2002, 2012, by = 2) )  + 
  theme( axis.ticks  = element_blank() )


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











