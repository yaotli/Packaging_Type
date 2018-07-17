source("functions.R")

require(ggplot2)
require(ggpubr)
require(skygrowth)
require(ape)
require(tidyverse)


# raw input -------- 
# 1 .tre
c232_D <- data.frame( read.table("./eco/results/201807_grid_eco/table/20180702_grid_eco_232_h5_D-com", sep = "\t", header = T), Eco = "Domestic", Clade = "232" )
c232_W <- data.frame( read.table("./eco/results/201807_grid_eco/table/20180702_grid_eco_232_h5_W-com", sep = "\t", header = T), Eco = "Wild", Clade = "232" )
c234_D <- data.frame( read.table("./eco/results/201807_grid_eco/table/20180706_grid_eco_234_h5_D-com", sep = "\t", header = T), Eco = "Domestic", Clade = "234" )
c234_W <- data.frame( read.table("./eco/results/201807_grid_eco/table/20180706_grid_eco_234_h5_W-com", sep = "\t", header = T), Eco = "Wild", Clade = "234" )


tb = list( c232_D, c232_W, c234_D, c234_W ) 
tb = do.call( rbind, tb )


# 2 .nwk 

c232_D_tre <- "./eco/results/201807_grid_eco/20180702_grid_eco_232_h5_D-anno.nwk"
c232_W_tre <- "./eco/results/201807_grid_eco/20180702_grid_eco_232_h5_W-anno.nwk"
c234_D_tre <- "./eco/results/201807_grid_eco/20180706_grid_eco_234_h5_D-anno.nwk"
c234_W_tre <- "./eco/results/201807_grid_eco/20180706_grid_eco_234_h5_W-anno.nwk"

t_c232_D <- 2011.989
t_c232_W <- 2011.948
t_c234_D <- 2011.953
t_c234_W <- 2009.496



# skygrowth --------
tb.g = skygrowth_df( ls.nwk = c( c232_D_tre, c232_W_tre, c234_D_tre, c234_W_tre ), 
                     ls.t   = c( t_c232_D, t_c232_W, t_c234_D, t_c234_W ), 
                     name   = c( "232D", "232W", "234D", "234W"), n.res = 50) 



# figures --------

p1 <- 
ggplot( data = tb ) + theme_bw() +
  facet_grid( Eco~. ) + 
  geom_line( aes( x = Time, y = log(Median), color = Clade ), size = 1.5) + 
  geom_ribbon( aes( x = Time, ymin = log(Lower), ymax = log(Upper), group = Clade), alpha = 0.1) +
  coord_cartesian( ylim = c(-1.5, 6), xlim = c(2003, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  #scale_linetype_manual(  values =  c( "solid", "dotted") ) +
  xlab( "Year" ) + ylab( "Effective population size" ) + 
  scale_color_manual( values = pyCol( c( "green", "red" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none", 
         strip.text = element_blank() ) +
  geom_vline( xintercept = 2007, linetype = "dotted") 


# growth rate figure

tb.g[, 7] <- ifelse( startsWith( as.character(tb.g[, 5]) , prefix = "232" ),  "232", "234" )
colnames(tb.g)[7] = "Clade"
tb.g[, 8] <- ifelse( endsWith( as.character(tb.g[, 5]) , suffix = "D" ),  "D", "W" )
colnames(tb.g)[8] = "Eco"

p2 <- 
  tb.g %>%
  filter( type == "rate" ) %>% 
  select( lb, e, ub, time, note, Clade, Eco) %>%
  ggplot() + theme_bw() +
  facet_grid( Eco~.) +
  geom_line( aes( x = time, y = e, color = Clade) , size = 1.5) + 
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub, group = Clade), alpha = 0.1) + 
  coord_cartesian( ylim = c(-6, 8), xlim = c(2003, 2012) ) +
  # scale_linetype_manual(  values =  c( "solid", "dotted") ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "Year" ) + ylab( "Growth rate" ) + 
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         strip.text = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" )  +
  geom_vline( xintercept = 2007, linetype = "dotted")


ggarrange( p1, p2, ncol = 2 )


