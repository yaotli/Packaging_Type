source("functions.R")

require(ggplot2)
require(ggpubr)
require(tidyverse)


# readin 

test = read.table( "~/Desktop/t2", header = TRUE)
test = test[ ,c(1,2,4)]

rawtb <- read.csv( "./dynamics/dynamics.csv", header = TRUE, stringsAsFactors = FALSE )


# violin plot 

v1 <- 
test %>% 
  gather( method, value, sr:ur ) %>% 
  ggplot( aes( x = method, y = 2013.989 - value  ) ) + 
  geom_violin( ) + 
  theme_bw() +
  stat_summary( fun.data = mean_sdl, geom = "pointrange" ) + 
  scale_y_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) + 
  scale_x_discrete( labels = c( "uncorrelated", "strict" ), breaks = c( "ur", "sr") ) +
  coord_flip( ylim = c(2003, 2012) ) 
  
  
# skygrowth 

l1 <- 
rawtb %>%
  filter( type == "mcmc" ) %>%
  filter( note == 234 )    %>%
  select( e, time, note )   %>%
  ggplot(  ) + theme_bw() + 
  geom_line( aes( x = time, y = e ), size = 2) + 
  xlab("") + 
  coord_cartesian( ylim = c(-6, 9), xlim = c(2003, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  theme( axis.title.x = element_blank(), 
         axis.text.x  = element_blank(),
         panel.grid.minor.y = element_blank())




# combine 

ggarrange( l1, v1, nrow = 2,  heights = c(1,2), align = "v")






