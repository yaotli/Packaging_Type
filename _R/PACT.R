library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)
library(tidyr)


source("~/Packaging_Type/_R/Function.R")
setwd("/Volumes/EDGE2/LoVE/ReassortSubtype/BEAST/")


a = read.table( "~/PACT-master/out.232.0102.skylines", header = T, stringsAsFactors = F )

p.a <- 
ggplot( data = a, aes( x = time, y = mean )) + 
  geom_area( aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2003,2015, by = 2), limits = c(2004.1, 2016), 
                      expand = c(0, 0) ) +
  scale_y_continuous( expand = c(0,0) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks       = element_blank(), 
         panel.border     = element_blank(), 
         axis.title       = element_text(size = 24),
         axis.text.y      = element_text(size = 12),
         axis.text.x      = element_text(size = 20),
         #axis.text.x = element_blank(),
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )


b = read.table( "~/PACT-master/out.234.0102-0.1.skylines", header = T, stringsAsFactors = F )

p.b <- 
ggplot( data = b, aes( x = time, y = mean )) + 
  geom_area( size = 1, aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2003,2015, by = 2), limits = c(2004.1, 2016), 
                      expand = c(0, 0) ) +
  scale_y_continuous( expand = c(0,0) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks       = element_blank(), 
         panel.border     = element_blank(), 
         axis.title       = element_text(size = 24),
         axis.text.y      = element_text(size = 12),
         axis.text.x      = element_text(size = 20),
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )

multiplot(p.b, p.a, ncol = 1)


### persistence 0122 --------------------------------

PACT_232 <- read.table("~/PACT-master/out.232.0120.persis.stats", header = TRUE)
PACT_234 <- read.table("~/PACT-master/out.234.0120.persis.stats", header = TRUE)

# 232
g1 <- ggplot(data = PACT_232[-1,], aes( x= statistic, y = mean, color = statistic)) +
  geom_point(size = 5) +
  geom_errorbar( aes( ymin = lower, ymax = upper ), width = 0.01, size = 1) +
  coord_flip() + 
  xlab("") + ylab("Persistence (year)") +
  scale_color_manual( values = c( pyCol( c( "brown", "green", "gray", "red", "blue" ) ) ) ) + 
  
  scale_x_discrete( labels = c("Central China", "Eastern China", "Northern China", "Southern China", 
                               "Southwest China") ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y     = element_blank(),
         panel.border     = element_blank(),
         axis.text.y      = element_text(size = 16),
         axis.text.x      = element_text(size = 15),
         axis.title.x     = element_text(size = 16),
         legend.position  = "none") 

# 234
g2 <- ggplot(data = PACT_234[-1,], aes( x= statistic, y = mean, color = statistic)) +
  geom_point(size = 5) +
  geom_errorbar( aes( ymin = lower, ymax = upper ), width = 0.01, size = 1) +
  coord_flip() + 
  xlab("") + ylab("Persistence (year)") +
  scale_color_manual( values = c( pyCol( c( "brown", "green", "red", "blue" ) ) ) ) + 
  
  scale_x_discrete( labels = c("Central China", "Eastern China", "Southern China", 
                               "Southwest China"  ) ) +
  scale_y_continuous( limits = c(0,3) ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y     = element_blank(),
         panel.border     = element_blank(),
         axis.text.y      = element_text(size = 16),
         axis.text.x      = element_text(size = 15),
         axis.title.x     = element_text(size = 16),
         legend.position  = "none") 

multiplot(g2, g1, ncol = 1) #5*6


### geo-trunk 0122 --------------------------------

a = read.table( "~/PACT-master/out.232.0120.skylines", header = T, stringsAsFactors = F )

p.a <- 
  ggplot( data = a, aes( x = time, y = mean )) + 
  geom_area( aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2003,2015, by = 2), limits = c(2004.1, 2016), 
                      expand = c(0, 0) ) +
  scale_y_continuous( expand = c(0,0) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks       = element_blank(), 
         panel.border     = element_blank(), 
         axis.title       = element_text(size = 24),
         axis.text.y      = element_text(size = 12),
         axis.text.x      = element_text(size = 20),
         #axis.text.x = element_blank(),
         legend.position  = "none") + 
  scale_fill_manual( values = c( pyCol( c( "brown", "green", "gray", "red", "blue" ) ) ) ) 
  

b = read.table( "~/PACT-master/out.234.0120.skylines", header = T, stringsAsFactors = F )

p.b <- 
  ggplot( data = b, aes( x = time, y = mean )) + 
  geom_area( size = 1, aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2003,2015, by = 2), limits = c(2004.1, 2016), 
                      expand = c(0, 0) ) +
  scale_y_continuous( expand = c(0,0) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks       = element_blank(), 
         panel.border     = element_blank(), 
         axis.title       = element_text(size = 24),
         axis.text.y      = element_text(size = 12),
         axis.text.x      = element_text(size = 20),
         legend.position  = "none") + 
  scale_fill_manual( values = c( pyCol( c( "brown", "green", "red", "blue" ) ) ) ) 

multiplot(p.b, p.a, ncol = 1) #4*6
