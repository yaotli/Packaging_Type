source("functions.R")

require(ggplot2)
require(ggpubr)
require(tidyverse)

# raw input -------- 
# 1 .tre
h5_232 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180629_grid_2012_232_h5-com", sep = "\t", header = T), Gene = "H5", Clade = "232" )
h5_234 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180629_grid_all_234_h5-com", sep = "\t", header = T), Gene = "H5", Clade = "234" )
h5_7   <- data.frame( read.table("./dynamics/results/201807_grid/table/20180706_grid_2012_7_h5-com", sep = "\t", header = T), Gene = "H5", Clade = "7" )
n1_232 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180629_grid_2012_232_n1-com", sep = "\t", header = T), Gene = "N1", Clade = "232" )
n1_234 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180629_grid_all_234_n1-com", sep = "\t", header = T), Gene = "N1", Clade = "234" )
n1_7   <- data.frame( read.table("./dynamics/results/201807_grid/table/20180706_grid_2012_7_n1-com", sep = "\t", header = T), Gene = "N1", Clade = "7" )


tb = list( h5_232, h5_234, h5_7, n1_232, n1_234, n1_7 ) 
tb = do.call( rbind, tb )

# 2 .nwk 

c232_nwk <- "./dynamics/results/201807_grid/20180629_grid_2012_232_h5-anno.nwk"
c234_nwk <- "./dynamics/results/201807_grid/20180629_grid_all_234_h5-anno.nwk"
c7_nwk   <- "./dynamics/results/201807_grid/20180706_grid_2012_7_h5-anno.nwk"
  
t_c232_D <- 2011.989
t_c234_D <- 2011.953
t_7_D    <- 2011.86


# skygrowth --------
tb.g = skygrowth_df( ls.nwk = c( c232_nwk, c234_nwk, c7_nwk ), 
                     ls.t   = c( t_c232_D, t_c234_D, t_7_D ), 
                     name   = c( "232", "234", "7" ), n.res = 50) 

# write.csv( tb.g, "./dynamics/dynamics.csv")


# figures --------
f1 <- 
tb %>%
  filter( Gene == "H5" ) %>%
  select( Clade, Median, Time, Upper, Lower ) %>%
  ggplot(  ) + theme_bw() +
  facet_grid( Clade~. ) + 
  geom_line( aes( x = Time, y = log(Median), color = Clade), size = 2) + 
  geom_ribbon( aes( x = Time, ymin = log(Lower), ymax = log(Upper)), alpha = 0.1) + 
  coord_cartesian( ylim = c(-1.6, 5.5), xlim = c(2003, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "Year" ) + ylab( "Effective population size" ) + 
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none", 
         strip.text = element_blank() ) +
  geom_vline( xintercept = 2007, linetype = "dotted")

  

# growth rate figure

dat_text <- data.frame( label = c("2.3.2", "2.3.4", "      7"),
                        note = c( "232", "234", "7" ) )

f2 <- 
tb.g %>%
  filter( type == "rate" ) %>% 
  select( lb, e, ub, time, note ) %>%
  ggplot() + theme_bw() +
  facet_grid( note~.) +
  geom_line( aes( x = time, y = e, color = note), size = 2) + 
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1) + 
  coord_cartesian( ylim = c(-6, 9), xlim = c(2003, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "Year" ) + ylab( "Growth rate" ) + 
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         strip.text = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" )  +
  geom_text( data = dat_text, mapping = aes(x = -Inf, y = -Inf, label = label ), 
             hjust = -7.5, vjust = -14) +
  geom_vline( xintercept = 2007, linetype = "dotted")



ggarrange( f1, f2, ncol = 2 )


