source("functions.R")

require(ggtree) 
require(tidyverse)
require(ggpubr)


c234_geo.tre <- "./bigmcc/results/20180805_geo_big_234_h5-1-anno.tre"
c232_geo.tre <- "./bigmcc/results/20180805_geo_big_232_h5-1-anno.tree"
c7_geo.tre   <- "./bigmcc/results/20180717_eco_big_7_h5-anno.tre"


# c234 tre

rawbeast.c234 <- read.beast( c234_geo.tre )

t1 <- 
ggtree( rawbeast.c234, right = TRUE, size = 0.5, mrsd = "2017-04-18" ) + 
  aes( color = states ) + 
  geom_tippoint(size = 0.5) +
  theme_tree2( axis.text.x = element_text(), 
               legend.title = element_blank(),
               legend.position = c(0.1,0.5) ) + 
  scale_x_continuous( breaks = seq(2002, 2017, by = 2), 
                      labels = seq(2002, 2017, by = 2),
                      limit  = c(2001.5, 2017.5) )  + 
  scale_y_continuous( expand = c(0,1) ) + 
  scale_color_manual( values =  c( "#d0694a",
                                   "#7f6e85",
                                   "#cc79a7",
                                   "#77bedb",
                                   "#ccc197",
                                   "#e1c72f" ) )

recdf <- data.frame( xstart = seq(2002, 2017, 2),
                     xend   = seq(2003, 2018, 2) )

t11 <- t1 + geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                       fill = "gray", alpha = 0.2, inherit.aes = FALSE) + 
       ggtitle( "clade 2.3.4" )


# c232 tre

rawbeast.c232 <- read.beast( c232_geo.tre )

t2 <- 
  ggtree( rawbeast.c232, right = TRUE, size = 0.5, mrsd = "2017-01-01" ) + 
  aes( color = states ) + 
  geom_tippoint(size = 0.5) +
  theme_tree2( axis.text.x = element_text(), 
               legend.title = element_blank(),
               legend.position = c(0.1,0.5) ) + 
  scale_x_continuous( breaks = seq(2002, 2017, by = 2), 
                      labels = seq(2002, 2017, by = 2),
                      limit  = c(2001.5, 2017.5) ) + 
  scale_y_continuous( expand = c(0,1) ) + 
  scale_color_manual( values =  c( "#48a365",
                                   "#7f6e85",
                                   "#cc79a7",
                                   "#77bedb",
                                   "#ccc197",
                                   "#e1c72f" ) )

recdf <- data.frame( xstart = seq(2002, 2017, 2),
                     xend   = seq(2003, 2018, 2) )

t22 <- t2 + geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                       fill = "gray", alpha = 0.2, inherit.aes = FALSE) + 
  ggtitle( "clade 2.3.2" )



# c7 tree

rawbeast.c7 <- read.beast( c7_geo.tre )

t3 <- 
  ggtree( rawbeast.c7, right = TRUE, size = 0.5, mrsd = "2015-01-15" ) + 
  aes( color = states ) + 
  geom_tippoint(size = 0.5) +
  theme_tree2( axis.text.x = element_text(), 
               legend.title = element_blank(),
               legend.position = c(0.1,0.5) ) + 
  scale_x_continuous( breaks = seq(2002, 2017, by = 2), 
                      labels = seq(2002, 2017, by = 2),
                      limit  = c(2001.5, 2017.5) ) + 
  
  scale_y_continuous( expand = c(0,1) ) + 
  scale_color_manual( values =  c( "#7f6e85","#ccc197" ) )

recdf <- data.frame( xstart = seq(2002, 2017, 2),
                     xend   = seq(2003, 2018, 2) )

t33 <- t3 + geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                       fill = "gray", alpha = 0.2, inherit.aes = FALSE) + 
  ggtitle( "clade 7" )

# combine

ggarrange( t11, t22, ggarrange( t33, nrow = 2, heights = c(0.6, 1) ), 
           ncol = 3, labels = c("A", "B", "C")) #a5


