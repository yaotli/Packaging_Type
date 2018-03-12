source("functions.R")

library(ggplot2)
require(ggpubr)


c232_ha_ride     <- read.table("./dynamics/ride_all/result/ride_all_232_h5pika-com", sep = "\t", header = T)
c232_na_ride     <- read.table("./dynamics/ride_all/result/ride_all_232_n1pika-com", sep = "\t", header = T)
c232_ha_ride2012 <- read.table("./dynamics/ride_2012/result/ride_2012_232_h5pika-com", sep = "\t", header = T)
c232_na_ride2012 <- read.table("./dynamics/ride_2012/result/ride_2012_232_n1pika-com", sep = "\t", header = T)

c234_ha_ride <- read.table("./dynamics/ride_all/result/ride_all_234_h5-com", sep = "\t", header = T)
c234_na_ride <- read.table("./dynamics/ride_all/result/ride_all_234_n1-com", sep = "\t", header = T)

c232_h5_D_ride <- read.table("./eco/eco_group/result/eco_group_232_h5pika_D", sep = "\t", header = T)
c232_h5_W_ride <- read.table("./eco/eco_group/result/eco_group_232_h5pika_W", sep = "\t", header = T)
c234_h5_D_ride <- read.table("./eco/eco_group/result/eco_group_234_h5_D", sep = "\t", header = T)
c234_h5_W_ride <- read.table("./eco/eco_group/result/eco_group_234_h5_Wsrd", sep = "\t", header = T)

c2344_ha_ride <- read.table("./Nx/result/h5nx_83-com", sep = "\t", header = T)


r1 <- 
  ggplot() + theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 12), 
         axis.title=element_text(size = 15) ) +
  
  scale_x_continuous(breaks = seq(2004, 2012, by=1) ) +
  xlab("") + ylab("Population size") +
  
  geom_line( data = c234_ha_ride, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_ha_ride2012, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  geom_ribbon( data = c234_ha_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_ha_ride2012, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 


r2 <- 
  ggplot() + theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 12), 
         axis.title=element_text(size = 15) ) +
  
  scale_y_continuous( breaks = c(0,2,4), labels = c(0,2,4) ) +
  scale_x_continuous(breaks = seq(2004, 2012, by=1) ) +
  xlab("") + ylab("Population size") +
  
  geom_line( data = c234_na_ride, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_na_ride2012, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  geom_ribbon( data = c234_na_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_na_ride2012, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 

r3 <- 
ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 12), 
         axis.title=element_text(size = 15) ) +
  
  scale_x_continuous(breaks = seq(2004, 2012, by=1), limits = c(2003.4, 2012)) +
  scale_y_continuous(limits = c(-1.5, 5) ) +
  
  xlab("") + ylab("Population size") +
  geom_line( data = c232_h5_D_ride, 
             aes(x = Time, y = log(Median) ), color = pyCol("brown"), size = 2) +
  geom_line( data = c232_h5_W_ride, 
             aes(x = Time, y = log(Median) ), color = pyCol("cyan"), size = 2) +
  geom_ribbon( data = c232_h5_D_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = pyCol("brown"), alpha = 0.1) + 
  geom_ribbon( data = c232_h5_W_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = pyCol("cyan"), alpha = 0.1) 

r4 <- 
  ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 12), 
         axis.title=element_text(size = 15) ) +
  
  scale_x_continuous(breaks = seq(2004, 2012, by=1), limits = c(2003.4, 2012)) +
  scale_y_continuous(limits = c(-1.5, 5) ) +
  
  xlab("") + ylab("Population size") +
  geom_line( data = c234_h5_D_ride, 
             aes(x = Time, y = log(Median) ), color = pyCol("brown"), size = 2) +
  geom_line( data = c234_h5_W_ride, 
             aes(x = Time, y = log(Median) ), color = pyCol("cyan"), size = 2) +
  geom_ribbon( data = c234_h5_D_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = pyCol("brown"), alpha = 0.1) + 
  geom_ribbon( data = c234_h5_W_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = pyCol("cyan"), alpha = 0.1) 

ggarrange( r1, r3, r2, r4, ncol = 2, nrow = 2, heights = c(1,1) ) #10*6


rall <-
  ggplot() + theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 12), 
         axis.title=element_text(size = 15) ) +
  
  scale_x_continuous(breaks = seq(2004, 2016, by=1) ) +
    xlab("") + ylab("Population size") +
  geom_line( data = c234_ha_ride, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  #geom_line( data = c232_ha_ride, 
  #           aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  geom_line( data = c2344_ha_ride, 
             aes(x = Time, y = log(Median) ), color = pyCol("blue"), size = 2) +
  geom_ribbon( data = c234_ha_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  #geom_ribbon( data = c232_ha_ride, 
  #             aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) +
  geom_ribbon( data = c2344_ha_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = pyCol("blue"), alpha = 0.1) +

  geom_vline( xintercept =  2013.989-6.5819, color = pyCol("blue"), size = 1) +
  geom_vline( xintercept =  2013.989-6.0049, color = pyCol("blue"), size = 0.5, linetype = "dashed") +
  geom_vline( xintercept =  2013.989-7.3194, color = pyCol("blue"), size = 0.5, linetype = "dashed") 
  
  #geom_point( aes( x = (2013.989-6.5819), y = 0 ), size = 4, color = "blue" ) + 
  #geom_errorbarh( aes( y    = 0,  x = (2013.989-6.5819),
  #                     xmax = 2013.989-6.0049, 
  #                     xmin = 2013.989-7.3194 ), color = "blue", height = 0, size = 1)

  
  
  
  

