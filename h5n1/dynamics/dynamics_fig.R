source("functions.R")

require(ggplot2)
require(ggpubr)
require(tidyverse)

# raw input -------- 
# 1 .tre
h5_232 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180723_grid_2014_232_h5-com", sep = "\t", header = T), Gene = "H5", Clade = "232" )
h5_234 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180629_grid_all_234_h5-com", sep = "\t", header = T), Gene = "H5", Clade = "234" )
h5_7   <- data.frame( read.table("./dynamics/results/201807_grid/table/20180719_grid_2014_7_h5-com", sep = "\t", header = T), Gene = "H5", Clade = "7" )
n1_232 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180723_grid_2014_232_n1-com", sep = "\t", header = T), Gene = "N1", Clade = "232" )
n1_234 <- data.frame( read.table("./dynamics/results/201807_grid/table/20180629_grid_all_234_n1-com", sep = "\t", header = T), Gene = "N1", Clade = "234" )
n1_7   <- data.frame( read.table("./dynamics/results/201807_grid/table/20180723_grid_2014_7_n1-com", sep = "\t", header = T), Gene = "N1", Clade = "7" )


tb = list( h5_232, h5_234, h5_7, n1_232, n1_234, n1_7 ) 
tb = do.call( rbind, tb )

# 2 .nwk 

c232_nwk    <- "./dynamics/results/201807_grid/20180723_grid_2014_232_h5-anno.nwk"
c234_nwk    <- "./dynamics/results/201807_grid/20180629_grid_all_234_h5-anno.nwk"
c7_nwk      <- "./dynamics/results/201807_grid/20180719_grid_2014_7_h5-anno.nwk"

c232.na_nwk <- "./dynamics/results/201807_grid/20180723_grid_2014_232_n1-anno.nwk"
c234.na_nwk <- "./dynamics/results/201807_grid/20180629_grid_all_234_n1-anno.nwk"
c7.na_nwk   <- "./dynamics/results/201807_grid/20180723_grid_2014_7_n1-anno.nwk"

  
t_c232_D <- 2013.879
t_c234_D <- 2011.953
t_7_D    <- 2013.153


# skygrowth --------
tb.g = skygrowth_df( ls.nwk = c( c232_nwk, c234_nwk, c7_nwk ), 
                     ls.t   = c( t_c232_D, t_c234_D, t_7_D ), 
                     name   = c( "232", "234", "7" ), n.res = 50) 


tb.g.na = skygrowth_df( ls.nwk = c( c232.na_nwk, c234.na_nwk, c7.na_nwk ), 
                        ls.t   = c( t_c232_D, t_c234_D, t_7_D ), 
                        name   = c( "232", "234", "7" ), n.res = 50) 


# write.csv( tb.g, "./dynamics/dynamics.csv")


# HA figures --------

# beast 

f1 <- 
tb %>%
  filter( Gene == "H5" ) %>%
  filter( 2003 <= Time ) %>% 
  filter( Time <= 2012 ) %>% 
  select( Clade, Median, Time, Upper, Lower ) %>%
  ggplot(  ) + theme_bw() +
  facet_grid( Clade~. ) + 
  geom_line( aes( x = Time, y = Median, color = Clade), size = 2) +
  geom_ribbon( aes( x = Time, ymin = Lower, ymax = Upper ), alpha = 0.1) +
  coord_cartesian( ylim = c( 0.1, 1000 ), xlim = c(2003, 2012) ) +
  
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))  ) +
  annotation_logticks() +
  
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "" ) + ylab( "Effective population size" ) + 
  ggtitle("Skygrid") +
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none", 
         strip.text = element_blank() ) +
  geom_vline( xintercept = 2007, linetype = "dotted")

  
# skygrowth 

f2 <- 
tb.g %>%
  filter( type == "mcmc" ) %>% 
  filter( 2003 <= time ) %>% 
  filter( time <= 2012 ) %>% 
  select( lb, e, ub, time, note ) %>%
  ggplot() + theme_bw() +
  facet_grid( note~.) +
  geom_line( aes( x = time, y = e, color = note), size = 2) + 
  coord_cartesian( ylim = c( 0.1, 1000 ), xlim = c(2003, 2012) ) +
  geom_ribbon( aes( x = time, ymin = lb , ymax = ub  ), alpha = 0.1) + 
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))  ) +
  annotation_logticks() +  
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "" ) + ylab( "Effective population size" ) + 
  ggtitle("Skygrowth") +
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         strip.text = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" )  +
  geom_vline( xintercept = 2007, linetype = "dotted")

# growth rate 

dat_text <- data.frame( label = c("2.3.2", "2.3.4", "      7"),
                        note = c( "232", "234", "7" ) )

f3 <- 
tb.g %>%
  filter( type == "rate" ) %>% 
  filter( 2003 <= time ) %>% 
  filter( time <= 2012 ) %>% 
  select( lb, e, ub, time, note ) %>%
  ggplot() + theme_bw() +
  facet_grid( note~.) +
  geom_line( aes( x = time, y = e, color = note), size = 2) + 
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1) + 
  coord_cartesian( ylim = c(-6, 8), xlim = c(2003, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "" ) + ylab( "Growth rate" ) +
  ggtitle(" ") +
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         strip.text = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" )  +
  geom_text( data = dat_text, mapping = aes(x = -Inf, y = -Inf, label = label ), 
             hjust = -2, vjust = -8) +
  geom_vline( xintercept = 2007, linetype = "dotted")



ggarrange( f1, f2, f3, ncol = 3, labels = c("A", "B", "C")) #4.5*8 inch landcape 


# NA figures --------

# beast 

n1 <- 
  tb %>%
  filter( Gene == "N1" ) %>%
  filter( 2003 <= Time ) %>% 
  filter( Time <= 2012 ) %>% 
  select( Clade, Median, Time, Upper, Lower ) %>%
  ggplot(  ) + theme_bw() +
  facet_grid( Clade~. ) + 
  geom_line( aes( x = Time, y = Median, color = Clade), size = 2) + 
  geom_ribbon( aes( x = Time, ymin = Lower, ymax = Upper ), alpha = 0.1) + 
  coord_cartesian( ylim = c( 0.1, 1000 ), xlim = c(2003, 2012) ) +
  
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))  ) +
  annotation_logticks() +
  
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "" ) + ylab( "Effective population size" ) + 
  ggtitle("Skygrid") +
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none", 
         strip.text = element_blank() ) +
  geom_vline( xintercept = 2007, linetype = "dotted")


# skygrowth 

n2 <- 
  tb.g.na %>%
  filter( type == "mcmc" ) %>% 
  filter( 2003 <= time ) %>% 
  filter( time <= 2012 ) %>% 
  select( lb, e, ub, time, note ) %>%
  ggplot() + theme_bw() +
  facet_grid( note~.) +
  geom_line( aes( x = time, y =  e , color = note), size = 2) + 
  geom_ribbon( aes( x = time, ymin = lb , ymax = ub  ), alpha = 0.1) + 
  coord_cartesian( ylim = c( 0.1, 1000 ), xlim = c(2003, 2012) ) +
  geom_ribbon( aes( x = time, ymin = lb , ymax = ub  ), alpha = 0.1) + 
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))  ) +
  annotation_logticks() +  
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "" ) + ylab( "Effective population size" ) + 
  ggtitle("Skygrowth") +
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         strip.text = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" )  +
  geom_vline( xintercept = 2007, linetype = "dotted")

# growth rate 

dat_text <- data.frame( label = c("2.3.2", "2.3.4", "      7"),
                        note = c( "232", "234", "7" ) )

n3 <- 
  tb.g.na %>%
  filter( 2003 <= time ) %>% 
  filter( time <= 2012 ) %>% 
  filter( type == "rate" ) %>% 
  select( lb, e, ub, time, note ) %>%
  ggplot() + theme_bw() +
  facet_grid( note~.) +
  geom_line( aes( x = time, y = e, color = note), size = 2) + 
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1) + 
  coord_cartesian( ylim = c(-6, 8), xlim = c(2002.5, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "" ) + ylab( "Growth rate" ) + 
  ggtitle(" ")+
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         strip.text = element_blank(),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" )  +
  geom_text( data = dat_text, mapping = aes(x = -Inf, y = -Inf, label = label ), 
             hjust = -7.5, vjust = -14) +
  geom_vline( xintercept = 2007, linetype = "dotted")



ggarrange( n1, n2, n3, ncol = 3, labels = c("A", "B", "C")) #4.5*8 inch landcape 

