source("./functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(ggpubr)
require(tidyverse)
require(HDInterval)


## figures TREE ----------------

eco_tre.232  <- "./eco/results/201807_eco_all/20180709_eco_2014_232_h5-anno.tre"
eco_tre.234  <- "./eco/results/201807_eco_all/20180709_eco_all_234_h5nx-anno.tre"
eco_tre.7    <- "./eco/results/201807_eco_all/20180723_eco_all_7_h5nx-anno.tre"
rawbeast.232 <- read.beast( eco_tre.232 )
rawbeast.234 <- read.beast( eco_tre.234 )
rawbeast.232 <- read.beast( eco_tre.232 )
rawbeast.7   <- read.beast( eco_tre.7 )


g1 <- 
  ggtree( rawbeast.232, right = TRUE, size = 0.6, mrsd = "2013-11-18") + 
  aes( color = states ) +
  scale_color_manual( values = c( "#8c564b", "#17becf" )  ) + 
  theme_tree2( axis.text.x = element_text() ) + 
  scale_x_continuous( breaks = seq(2002,2013, by = 2), 
                      labels = seq(2002,2013, by = 2), 
                      limit = c( 2002, 2014) )  +
  scale_y_continuous( expand = c(0, 1) ) 
  # geom_tippoint(aes(fill = eco), shape = 21, color = "black") + 

rectdf <- data.frame( xstart = seq( 2002, 2013, 2), 
                      xend   = seq( 2003, 2014, 2))

g11 <- g1 + geom_rect( data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                       fill = "gray", alpha = 0.2, inherit.aes=FALSE) +
            ggtitle("clade 2.3.2")

# zoom in 
#ggtree( rawbeast.232, right = TRUE, size = 1, mrsd = "2016-04-12") + 
#  geom_range(range='height_0.95_HPD', size = 1) +
#  geom_tippoint() +
#  aes( color = states ) +
#  scale_color_manual( values = c( "#8c564b", "#17becf" )  ) + 
#  theme_tree2( axis.text.x = element_text( size = 12  ), 
#               panel.grid.major.x = element_line(color = "gray")) + 
#  scale_y_continuous( expand = c(0,1), limits = c(80, 180) ) + 
  # ecom_tippoint(aes(fill = eco), shape = 21, color = "black") + 
#  scale_x_continuous( limits = c(2004, 2008) )  


g2 <- 
  ggtree( rawbeast.234, right = TRUE, size = 0.6, mrsd = "2013-12-18") + 
  aes( color = states ) +
  scale_color_manual( values = c( "#8c564b", "#17becf" )  ) + 
  theme_tree2( axis.text.x = element_text() ) +
  scale_x_continuous( breaks = seq(2004,2013, by = 2), 
                      labels = seq(2004,2013, by = 2), 
                      limit = c( 2003, 2014) )  +
  scale_y_continuous( expand = c(0, 1) ) 
  #ecom_tippoint(aes(fill = eco), shape = 21, color = "black") + 


rectdf <- data.frame( xstart = seq( 2002, 2014, 2), 
                      xend   = seq( 2003, 2015, 2))

g22 <- g2 + geom_rect( data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
                       fill = "gray", alpha = 0.2, inherit.aes=FALSE) + 
            ggtitle("clade 2.3.4")

g22 <- flip(g22, 384, 243)

# clade 7
g3 <- 
  ggtree( rawbeast.7, right = TRUE, size = 0.6, mrsd = "2013-12-21") + 
  aes( color = states ) +
  scale_color_manual( values = c( "#8c564b", "#17becf" )  ) + 
  theme_tree2( axis.text.x = element_text() ) +
  scale_x_continuous( breaks = seq(2002,2013, by = 2), 
                      labels = seq(2002,2013, by = 2), 
                      limit = c( 2002, 2014) )  +
  scale_y_continuous( expand = c(0, 1) ) 


rectdf <- data.frame( xstart = seq( 2002, 2014, 2), 
                      xend   = seq( 2003, 2015, 2))

g33 <- g3 + geom_rect( data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
                       fill = "gray", alpha = 0.2, inherit.aes=FALSE) + 
            ggtitle("clade 7")






## figures TRUNK ----------------


trunk.232 = read.table( "./eco/PACT/skylines.20180709_eco_2014_232_h5", header = T, stringsAsFactors = F )
trunk.234 = read.table( "./eco/PACT/skylines.20180709_eco_all_234_h5nx", header = T, stringsAsFactors = F )
trunk.7   = read.table( "./eco/PACT/skylines.20180723_eco_all_7_h5nx", header = T, stringsAsFactors = F )


t1 <- 
  ggplot( data = trunk.232, aes( x = time, y = mean )) + 
  geom_area( aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2004,2013, by = 2), 
                      labels = seq(2004,2013, by = 2), 
                      limits = c( 2002, 2014) ) +
  
  scale_y_continuous( expand = c(0,0.005) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.border     = element_blank(), 
         #axis.title       = element_text(size = 12),
         #axis.text.y      = element_text(size = 8),
         #axis.text.x      = element_text(size = 10),
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )



t2 <- 
  ggplot( data = trunk.234, aes( x = time, y = mean )) + 
  geom_area( size = 1, aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2004,2013, by = 2), 
                      labels = seq(2004,2013, by = 2),
                      limit = c( 2003, 2014) ) +
  
  scale_y_continuous( expand = c(0, 0.005) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.border     = element_blank(), 
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )


# clade 7

t3 <- 
  ggplot( data = trunk.7, aes( x = time, y = mean )) + 
  geom_area( aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2004,2013, by = 2), 
                      labels = seq(2004,2013, by = 2), 
                      limits = c( 2002, 2014) ) +
  
  scale_y_continuous( expand = c(0,0.005) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         panel.border     = element_blank(), 
         #axis.title       = element_text(size = 12),
         #axis.text.y      = element_text(size = 8),
         #axis.text.x      = element_text(size = 10),
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )



ggarrange( g11, g22, g33,  t1, t2, t3, ncol = 3, nrow = 2, heights = c(3,1), align = "v", 
           labels = c("A", "B", "C", "", "", "")) #5*9 landscape


## figures Persistence ----------------

per_232 <- read.table("./eco/PACT/persis.20180709_eco_2014_232_h5", header = TRUE)
per_234 <- read.table("./eco/PACT/persis.20180709_eco_all_234_h5nx", header = TRUE)
stater  <- read.table("./eco/results/eco_2012_stater", header = TRUE)


persis.df <- data.frame( clade = c("232", "234"), 
                          mean  = c( per_232$mean[3], per_234$mean[3] ),
                          lower = c( per_232$lower[3], per_234$lower[3] ), 
                          upper = c( per_232$upper[3], per_234$upper[3] ), 
                          data  = "persis" )
state.r  <- data.frame( clade = c("232", "234"), 
                        mean  = c( median(stater$X232), median(stater$X234) ),
                        lower = c( hdi(stater$X232, 0.95)[[1]], hdi(stater$X234, 0.95)[[1]] ), 
                        upper = c( hdi(stater$X232, 0.95)[[2]], hdi(stater$X234, 0.95)[[2]] ), 
                        data  = "stater" )

persis.df = rbind( persis.df, state.r )


p1 <- 
  persis.df %>%
  filter( data == "persis" ) %>%
  select( clade, mean, lower, upper ) %>%
  ggplot( aes( x = clade, y = mean, color = clade ) ) +
  geom_point( size = 5 )  +
  geom_errorbar( aes( x = clade, ymin = lower, ymax = upper ), width = 0, size = 1) +
  xlab("") + ylab("Persistence in wild states (year)") +
  ggtitle(" ") +
  scale_color_manual( values = pyCol( c("green", "red")  ) ) + 
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         axis.ticks.y     = element_blank(),
         axis.text.y      = element_text(size = 10),
         axis.text.x      = element_text(size = 8),
         axis.title.x     = element_text(size = 8),
         legend.position  = "none") 

# 234
p2 <- 
  persis.df %>%
  filter( data == "stater" ) %>%
  select( clade, mean, lower, upper ) %>%
  ggplot( aes( x = clade, y = mean, color = clade ) ) +
  geom_point( size = 5 )  +
  geom_errorbar( aes( x = clade, ymin = lower, ymax = upper ), width = 0, size = 1) +
  xlab("") + ylab("State transition rate") +
  ggtitle(" ") +
  scale_y_continuous( limits = c(0,1) ) +
  scale_color_manual( values = pyCol( c("green", "red")  ) ) + 
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         axis.ticks.y     = element_blank(),
         axis.text.y      = element_text(size = 10),
         axis.text.x      = element_text(size = 8),
         axis.title.x     = element_text(size = 8),
         legend.position  = "none") 


ggarrange( p2, p1, labels = c("A", "B"), widths = c(1, 0.9) ) #3*5