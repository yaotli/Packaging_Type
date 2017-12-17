library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)
library(tidyr)

source("~/Packaging_Type/_R/Function.R")
setwd("/Volumes/EDGE2/LoVE/ReassortSubtype/BEAST/")


### skygrid_0810 --------------------------------

## 234 HA ----------------
c234_ha_mcc  <- read.beast("./skygrid_0810/result_0810/clade234_h5_ucrl_grid_100g_6e7.trees.tre")
c234_ha_grid <- read.table("./skygrid_0810/result_0810/clade234_h5_grid.csv", sep = "\t", header = T)

g234_ha <- 
ggtree(c234_ha_mcc, mrsd = "2012-12-21")  + theme_tree2() + 
  scale_x_continuous(breaks = seq(2002, 2013, by=1), limit = c(2002, 2013)) + 
  theme_tree2(panel.grid.major.x = element_line(color = "gray", size = 0.2),
              panel.grid.minor.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank() )

tip2344lile <- findtaxa(type     = 0,
                        tree     = c234_ha_mcc, 
                        targetid = c("EU195400_mallard_Huadong_lk_2005_H5N1_2005.496", 
                                   "KP233704_duck_Hunan_316_2005_H5N1_2005.301",
                                   "HM172100_duck_Jiangxi_80_2005_H5N1_2005.496",
                                   "DQ992838_crested_myna_Hong_Kong_540_2006_H5N1_2006.496",
                                   "DQ992790_duck_Hunan_324_2006_H5N1_2006.496"), 
                        target   = rep("red", 5) )

g234_ha <- g234_ha %<+% tip2344lile + 
  geom_tippoint(aes(color = I(shapee)), size = 1, alpha = 0.8)

g234_ha_grid <- 
  g234_ha + 
  geom_line(data = c234_ha_grid, 
            aes( x = Time, y = log(Median)*10 + 50 ), 
            color = "gray", alpha = 0.4, size = 2 ) + 
  geom_ribbon(data = c234_ha_grid, 
              aes(x = Time, ymin = log(Lower)*10 + 50, ymax = log(Upper)*10 + 50),  
              inherit.aes=FALSE, fill = "gray", alpha = 0.2) + ggtitle("HA - Clade 234*") + 
  geom_vline(xintercept = 2008, size = 1, color = "black", linetype = "dotted") +
  geom_point(x = 2007.044, y = 0, size = 0.8, color = "blue") + 
  geom_errorbarh( aes(y = 0, 
                     xmax = 2007.849, 
                     xmin = 2005.558), color = "blue")

## MRCA of Nx virus (rate_0818) ----------------

rate_0818_log <- data.frame( method = c("strict_exp", "strict_skyride", "relaxed_exp", "relaxed_skyride"),
                             median = c(3.0269, 3.0468, 2.9611, 3.0233),
                             hpd_l  = c(2.2218, 2.2994, 1.7829, 1.9933),
                             hpd_u  = c(4.5126, 4.2755, 10.0257, 5.2335), stringsAsFactors = FALSE)

ggplot(data = rate_0818_log[-3,], aes( x= method, y = 2010.071-median, color = method)) +
  geom_point(size = 4) +
  geom_errorbar( aes( ymin = 2010.071-hpd_l, ymax = 2010.071-hpd_u), width = 0.1) +
  coord_flip() + 
  xlab("") + ylab("") +
  scale_color_manual(values = c("black", "blue", "black") ) + 
  scale_y_continuous(breaks = seq(2000,2008, by = 1)) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         text = element_text(size = 28, face = "bold"), 
         legend.position = "none") 
  

## 234 NA ----------------

c234_na_mcc  <- read.beast("./skygrid_0810/result_0810/clade234_n1_ucrl_grid_100g_6e7.trees.tre")
c234_na_grid <- read.table("./skygrid_0810/result_0810/clade234_n1_grid.csv", sep = "\t", header = T)

g234_na <- 
  ggtree(c234_na_mcc, mrsd = "2012-12-21")  + theme_tree2() + 
  scale_x_continuous(breaks = seq(2002, 2013, by=1), limit = c(2002, 2013)) + 
  theme_tree2(panel.grid.major.x = element_line(color = "gray", size = 0.2),
              panel.grid.minor.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank() )

tip2344lile <- findtaxa(type     = 0,
                        tree     = c234_na_mcc, 
                        targetid = c("EU195402_mallard_Huadong_lk_2005_H5N1_2005.496", 
                                     "KP233706_duck_Hunan_316_2005_H5N1_2005.301",
                                     "HM172200_duck_Jiangxi_80_2005_H5N1_2005.496",
                                     "EF124202_crested_myna_Hong_Kong_540_2006_H5N1_2006.496",
                                     "EF124253_duck_Hunan_324_2006_H5N1_2006.496"), 
                        target   = rep("red", 5) )

g234_na <- g234_na %<+% tip2344lile + 
  geom_tippoint(aes(color = I(shapee)), size = 1, alpha = 0.8)

g234_na_grid <- 
  g234_na + 
  geom_line(data = c234_na_grid, 
            aes( x = Time, y = log(Median)*10 + 50 ), 
            color = "gray", alpha = 0.4, size = 2 ) + 
  geom_ribbon(data = c234_na_grid, 
              aes(x = Time, ymin = log(Lower)*10 + 50, ymax = log(Upper)*10 + 50),  
              inherit.aes=FALSE, fill = "gray", alpha = 0.2) + ggtitle("NA - Clade 234*") 

# multiplot(g234_ha_grid, g234_na_grid, ncol = 1)


## 232 HA ----------------

c232_ha_mcc  <- read.beast("./skygrid_0810/result_0810/clade232_h5_ucrl_grid_100g_6e7.trees.tre")
c232_ha_grid <- read.table("./skygrid_0810/result_0810/clade232_h5_grid.csv", sep = "\t", header = T)

g232_ha <- 
  ggtree(c232_ha_mcc, mrsd = "2013-11-18")  + theme_tree2() + 
  scale_x_continuous(breaks = seq(1998, 2013, by=1)) + 
  theme_tree2(panel.grid.major.x = element_line(color = "gray", size = 0.2),
              panel.grid.minor.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank() )


g232_ha_grid <- 
  g232_ha + 
  geom_line(data = c232_ha_grid, 
            aes( x = Time, y = log(Median)*10 + 50 ), 
            color = "gray", alpha = 0.4, size = 2 ) + 
  geom_ribbon(data = c232_ha_grid, 
              aes(x = Time, ymin = log(Lower)*10 + 50, ymax = log(Upper)*10 + 50),  
              inherit.aes=FALSE, fill = "gray", alpha = 0.2) + ggtitle("HA - Clade 232")



## 232 NA ----------------

c232_na_mcc  <- read.beast("./skygrid_0810/result_0810/clade232_n1_ucrl_grid_100g_6e7.trees.tre")
c232_na_grid <- read.table("./skygrid_0810/result_0810/clade232_n1_grid.csv", sep = "\t", header = T)

g232_na <- 
  ggtree(c232_na_mcc, mrsd = "2013-11-18")  + theme_tree2() + 
  scale_x_continuous(breaks = seq(1998, 2013, by=1), limit = c(1998, 2014)) + 
  theme_tree2(panel.grid.major.x = element_line(color = "gray", size = 0.2),
              panel.grid.minor.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank() )


g232_na_grid <- 
  g232_na + 
  geom_line(data = c232_na_grid, 
            aes( x = Time, y = log(Median)*10 + 50 ), 
            color = "gray", alpha = 0.4, size = 2 ) + 
  geom_ribbon(data = c232_na_grid, 
              aes(x = Time, ymin = log(Lower)*10 + 50, ymax = log(Upper)*10 + 50),  
              inherit.aes=FALSE, fill = "gray", alpha = 0.2) + ggtitle("NA - Clade 232") +
  scale_y_continuous(limit = c(0, 140))





## 7 HA ----------------

c7_ha_mcc  <- read.beast("./skygrid_0810/result_0810/clade7_h5_ucrl_grid_100g_6e7.trees.tre")
c7_ha_grid <- read.table("./skygrid_0810/result_0810/clade7_h5_grid.csv", sep = "\t", header = T)

g7_ha <- 
  ggtree(c7_ha_mcc, mrsd = "2013-02-26")  + theme_tree2() + 
  scale_x_continuous(breaks = seq(1998, 2013, by=1), limit = c(1998, 2014)) + 
  theme_tree2(panel.grid.major.x = element_line(color = "gray", size = 0.2),
              panel.grid.minor.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank() )


g7_ha_grid <- 
  g7_ha + 
  geom_line(data = c7_ha_grid, 
            aes( x = Time, y = log(Median)*10 ), 
            color = "gray", alpha = 0.4, size = 2 ) + 
  geom_ribbon(data = c7_ha_grid, 
              aes(x = Time, ymin = log(Median)*10, ymax = log(Upper)*10 ),  
              inherit.aes=FALSE, fill = "gray", alpha = 0.2) + ggtitle("HA - Clade 7") + 
  geom_vline(xintercept = 2009, size = 1, color = "black", linetype = "dotted")


## 7 NA ----------------

c7_na_mcc  <- read.beast("./skygrid_0810/result_0810/clade7_n1_ucrl_grid_100g_6e7.trees.tre")
c7_na_grid <- read.table("./skygrid_0810/result_0810/clade7_n1_grid.csv", sep = "\t", header = T)

g7_na <- 
  ggtree(c7_na_mcc, mrsd = "2013-02-26")  + theme_tree2() + 
  scale_x_continuous(breaks = seq(1998, 2013, by=1), limit = c(1998, 2014)) + 
  theme_tree2(panel.grid.major.x = element_line(color = "gray", size = 0.2),
              panel.grid.minor.x = element_line(color = "grey", size = 0.1),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank() )


g7_na_grid <- 
  g7_na + 
  geom_line(data = c7_na_grid, 
            aes( x = Time, y = log(Median)*10), 
            color = "gray", alpha = 0.4, size = 2 ) + 
  geom_ribbon(data = c7_na_grid, 
              aes(x = Time, ymin = log(Median)*10, ymax = log(Upper)*10),  
              inherit.aes=FALSE, fill = "gray", alpha = 0.2) + ggtitle("NA - Clade 7") 

  
  
# combine 
multiplot(g234_ha_grid, g232_ha_grid,
          g7_ha_grid, g234_na_grid,
          g232_na_grid, g7_na_grid, ncol = 2)




### mean rate from BEAST --------------------------------

# HA 
ha_meanRate <- data.frame( t(read.table("./skygrid_0810/result_0810/ha_meanRate.csv", 
                          sep = "\t", header = TRUE) ), stringsAsFactors = FALSE)

colnames(ha_meanRate) <- ha_meanRate[1, ][1:10]
ha_meanRate           <- ha_meanRate[-1,]
ha_meanRate$mean      <- as.numeric(ha_meanRate$mean)

Upper = as.numeric( unlist( strsplit(gsub("\\[|\\]", "", ha_meanRate$`95% HPD Interval`), ", "))[c(2,4,6)] )
Lower = as.numeric( unlist( strsplit(gsub("\\[|\\]", "", ha_meanRate$`95% HPD Interval`), ", "))[c(1,3,5)] )

ha_meanRate <- data.frame(ha_meanRate, Upper, Lower)

ggplot(data = ha_meanRate, aes(x=c("232", "234", "7"), y = mean)) +
  geom_point() +
  geom_errorbar( aes(ymin = Lower, ymax = Upper), width = 0.1) +
  xlab("") + ylab("meanRate") +
  theme_bw() + ggtitle("HA substitution rate")



# NA
na_meanRate <- data.frame( t( read.table("./skygrid_0810/result_0810/na_meanRate.csv", 
                                           sep = "\t", header = TRUE) ), stringsAsFactors = FALSE)

colnames(na_meanRate) <- na_meanRate[1, ][1:10]
na_meanRate           <- na_meanRate[-1,]
na_meanRate$mean      <- as.numeric(na_meanRate$mean)

Upper = as.numeric( unlist( strsplit(gsub("\\[|\\]", "", na_meanRate$`95% HPD Interval`), ", "))[c(2,4,6)] )
Lower = as.numeric( unlist( strsplit(gsub("\\[|\\]", "", na_meanRate$`95% HPD Interval`), ", "))[c(1,3,5)] )

na_meanRate <- data.frame(na_meanRate, Upper, Lower)

ggplot(data = na_meanRate, aes(x=c("232", "234", "7"), y = mean)) +
  geom_point() +
  geom_errorbar( aes(ymin = Lower, ymax = Upper), width = 0.1) +
  xlab("") + ylab("meanRate") +
  theme_bw() + ggtitle("NA substitution rate")


### prepare combined tree --------------------------------

c("#FF4040", "#66CD00", "#000000")
c("#4F94CD")
rectdf <- data.frame( xstart = seq(2001, 2014, 2), xend = seq(2002, 2014, 2))

ggtree(c234_ha_mcc, mrsd = "2012-12-21", color = "#FF4040", size = 1)  + theme_tree2() + 
  scale_x_continuous(breaks = seq(2001, 2014, by=1), limit = c(2001, 2014)) + 
  geom_rect(data = rectdf, 
            aes(xmin = xstart, xmax = xend, 
                ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2, inherit.aes=FALSE)

ggtree(c232_ha_mcc, mrsd = "2013-11-18", color = "#66CD00", size = 1)  + theme_tree2() + 
  scale_x_continuous(breaks = seq(2001, 2014, by=1), limit = c(2001, 2014)) + 
  geom_rect(data = rectdf, 
            aes(xmin = xstart, xmax = xend, 
                ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2, inherit.aes=FALSE)

ggtree(c7_ha_mcc, mrsd = "2013-02-26", size = 1)  + theme_tree2() + 
  scale_x_continuous(breaks = seq(2001, 2014, by=1), limit = c(2001, 2014)) + 
  geom_rect(data = rectdf, 
            aes(xmin = xstart, xmax = xend, 
                ymin = -Inf, ymax = Inf), fill = "gray", alpha = 0.2, inherit.aes=FALSE)

  
### combined grid / ride --------------------------------

c234_ha_grid <- read.table("./skygrid_0810/result_0810/clade234_h5_grid.csv", sep = "\t", header = T)
c234_na_grid <- read.table("./skygrid_0810/result_0810/clade234_n1_grid.csv", sep = "\t", header = T)

c232_ha_grid <- read.table("./skygrid_0810/result_0810/clade232_h5_grid.csv", sep = "\t", header = T)
c232_na_grid <- read.table("./skygrid_0810/result_0810/clade232_n1_grid.csv", sep = "\t", header = T)

c7_ha_grid   <- read.table("./skygrid_0810/result_0810/clade7_h5_grid.csv", sep = "\t", header = T)
c7_na_grid   <- read.table("./skygrid_0810/result_0810/clade7_n1_grid.csv", sep = "\t", header = T)

N5_8_ride    <- read.table("./rate_0818/result/strict_gmrf_n5", sep = "\t", header = T) 


# HA 
combined_ride_ha <- 
ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 18), 
         axis.title=element_text(size = 16, face="bold") ) +
  
  scale_x_continuous(breaks = seq(2004, 2014, by=1) ) +
  
  xlab("") + ylab("Relative population size") +
  
  geom_line( data = c234_ha_grid, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_ha_grid, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  geom_line( data = N5_8_ride, 
             aes(x = Time, y = log(Median) ), color = "darkblue", size = 1) +
  
  #geom_line( data = c7_ha_grid, 
  #           aes(x = Time, y = log(Median) ), color = "gray", size = 1.5) +
  

  geom_ribbon( data = c234_ha_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_ha_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) +
  
  geom_ribbon( data = N5_8_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "blue", alpha = 0.1) 
  
  #geom_ribbon( data = c7_ha_grid, 
  #             aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "gray", alpha = 0.1) 
  

# NA
combined_ride_na <- 
ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 18), 
         axis.title=element_text(size = 16, face="bold") ) +
  
  scale_x_continuous(breaks = seq(2004, 2014, by=1) ) +
  scale_y_continuous(breaks = seq(-3, 5, by=2) ) +
  
  xlab("") + ylab("Relative population size") +
  
  geom_line( data = c234_na_grid, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_na_grid, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  
  #geom_line( data = c7_na_grid, 
  #           aes(x = Time, y = log(Median) ), color = "gray", size = 1.5) +
  
  
  geom_ribbon( data = c234_na_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_na_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 
  
  #geom_ribbon( data = c7_na_grid, 
  #             aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "gray", alpha = 0.1) 



multiplot(combined_ride_ha, combined_ride_na, ncol = 1)



### mean rate dist --------------------------------

rate_dist_234ha <- read.table("./skygrid_0810/result_0810/ratedist_ha_234.csv", sep = "\t", header = T)
rate_dist_232ha <- read.table("./skygrid_0810/result_0810/ratedist_ha_232.csv", sep = "\t", header = T)
rate_dist_7ha   <- read.table("./skygrid_0810/result_0810/ratedist_ha_7.csv", sep = "\t", header = T)
rate_dist_n5    <- read.table("./skygrid_0810/result_0810/ratedist_ha_n5_strict_gmrf.csv", sep = "\t", header = T)

rate_dist_234n1 <- read.table("./skygrid_0810/result_0810/ratedist_n1_234.csv", sep = "\t", header = T)
rate_dist_232n1 <- read.table("./skygrid_0810/result_0810/ratedist_n1_232.csv", sep = "\t", header = T)
rate_dist_7n1   <- read.table("./skygrid_0810/result_0810/ratedist_n1_7.csv", sep = "\t", header = T)

# HA 

ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y = element_blank(), 
         panel.border = element_blank(),
         axis.text.y = element_blank(), 
         axis.text.x = element_text(size = 20)) +
  
  scale_x_continuous(limit = c(0., 0.008)) +
  
  xlab("") + 
  
  geom_ribbon( data = rate_dist_7ha, 
               aes(x = clade7_h5, ymin = 0 , ymax = meanRate ), 
               fill = "gray", alpha = 0.4) + 
    
  geom_ribbon( data = rate_dist_232ha, 
               aes(x = clade232_h5, ymin = 0 , ymax = meanRate ), 
               fill = "#2ca02c", alpha = 0.4) + 
  
  geom_ribbon( data = rate_dist_234ha, 
               aes(x = clade234_h5, ymin =  0, ymax = meanRate ), 
               fill = "#d62728", alpha = 0.4) + 
  
  geom_ribbon( data = rate_dist_n5, 
               aes(x = h5_n5, ymin =  0, ymax = meanRate ), 
               fill = "blue", alpha = 0.1) 
  
# NA 

ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y = element_blank(), 
         panel.border = element_blank(),
         axis.text.y = element_blank(), 
         axis.text.x = element_text(size = 20)) +
  
  scale_x_continuous(limit = c(0., 0.008)) +
  
  xlab("") + 
  
  geom_ribbon( data = rate_dist_7n1, 
               aes(x = clade7_n1, ymin = 0 , ymax = meanRate ), 
               fill = "gray", alpha = 0.4) + 
  
  geom_ribbon( data = rate_dist_232n1, 
               aes(x = clade232_n1, ymin = 0 , ymax = meanRate ), 
               fill = "#2ca02c", alpha = 0.4) + 
  
  geom_ribbon( data = rate_dist_234n1, 
               aes(x = clade234_n1, ymin =  0, ymax = meanRate ), 
               fill = "#d62728", alpha = 0.4)
  
  
### more grid / ride / rate plot --------------------------------

c232_ha_grid <- read.table("./grid_1107/result/232_h5_grid_1107.csv", sep = "\t", header = T)
c234_ha_grid <- read.table("./grid_1107/result/234_h5_grid_1107.csv", sep = "\t", header = T)
c232_na_grid <- read.table("./grid_1107/result/232_n1_grid_1107.csv", sep = "\t", header = T)
c234_na_grid <- read.table("./grid_1107/result/234_n1_grid_1107.csv", sep = "\t", header = T)

c232_ha_ride <- read.table("./ride_1205/result/ride_232_h5_1205", sep = "\t", header = T)
c234_ha_ride <- read.table("./ride_1205/result/ride_234_h5_1205", sep = "\t", header = T)
c232_na_ride <- read.table("./ride_1205/result/ride_232_n1_1205", sep = "\t", header = T)
c234_na_ride <- read.table("./ride_1205/result/ride_234_n1_1205", sep = "\t", header = T)


# grid - HA

combined_grid_ha <- 
  ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 18), 
         axis.title=element_text(size = 16, face = "bold") ) +
  
  scale_x_continuous(breaks = seq(2004, 2014, by=1) ) +
  
  xlab("") + ylab("Population size") +
  
  geom_line( data = c234_ha_grid, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_ha_grid, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  
  geom_ribbon( data = c234_ha_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_ha_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 
  
 
# grid - NA 

combined_grid_na <- 
  ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 18), 
         axis.title=element_text(size = 16, face = "bold") ) +
  
  scale_x_continuous(breaks = seq(2004, 2014, by=1) ) +
  scale_y_continuous(breaks = seq(-3, 5, by=2) ) +
  
  xlab("") + ylab("Population size") +
  
  geom_line( data = c234_na_grid, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_na_grid, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +

  geom_ribbon( data = c234_na_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_na_grid, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 

multiplot( combined_grid_ha, combined_grid_na, ncol = 1)

# ride - HA 

combined_ride_ha <- 
  ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text(size = 14), 
         axis.title=element_text(size = 16, face = "bold") ) +
  
  scale_x_continuous(breaks = seq(2004, 2016, by=1) ) +
  
  xlab("") + ylab("Population size") +
  
  geom_line( data = c234_ha_ride, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_ha_ride, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  
  geom_ribbon( data = c234_ha_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_ha_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 



# ride - NA

combined_ride_na <- 
  ggplot() + theme_bw() + 
  
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         axis.text.x = element_text( size = 14), 
         axis.title=element_text(size = 16, face ="bold") ) +
  
  scale_x_continuous(breaks = seq(2004, 2016, by=1) ) +
  scale_y_continuous(breaks = seq(-3, 5, by=2) ) +
  
  xlab("") + ylab("Population size") +
  
  geom_line( data = c234_na_ride, 
             aes(x = Time, y = log(Median) ), color = "#d62728", size = 2) +
  geom_line( data = c232_na_ride, 
             aes(x = Time, y = log(Median) ), color = "#2ca02c", size = 2) +
  
  geom_ribbon( data = c234_na_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#d62728", alpha = 0.1) + 
  geom_ribbon( data = c232_na_ride, 
               aes(x = Time, ymin = log(Lower) , ymax = log(Upper) ), fill = "#2ca02c", alpha = 0.1) 

multiplot( combined_ride_ha, combined_ride_na, ncol = 1)


# rate 

rate_h5n1 <- data.frame( clade  = c( rep( c( rep("232", 2), rep("234", 2) ), 2 ) ),
                         gene   = c( rep( c("HA", "NA"), 4) ),
                         method = c( rep("ride", 4), rep("grid", 4) ),
                         median = c( 5.3323E-3, 4.9161E-3, 3.8127E-3, 3.7575E-3, 5.6302E-3, 4.8016E-3, 4.0174E-3, 3.8179E-3),
                         hpd_l  = c( 4.586E-3, 4.1466E-3, 3.1808E-3, 3.1799E-3, 4.8245E-3, 3.8923E-3, 3.5E-3, 3.2948E-3 ),
                         hpd_u  = c( 6.0966E-3, 5.7407E-3, 4.4604E-3, 4.3539E-3, 6.5103E-3, 5.7716E-3, 4.6053E-3, 4.3723E-3 ), 
                         stringsAsFactors = FALSE)

ggplot(data = rate_h5n1, aes(x = method, y = median, color = clade )) +
  geom_point( size = 4, position = position_dodge(0.8)) +
  geom_errorbar( aes(ymin = hpd_l, ymax = hpd_u), width = 0.01, size = 1.2, position = position_dodge(0.8)) +
  facet_wrap(~gene, ncol = 2) + 
  xlab("") + ylab("meanRate") +
  scale_color_manual(values = c( "#00BFC4", "#F8766D" )) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 20), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        legend.text = element_text(size = 12), 
        legend.title = element_blank(), 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 18, face = "bold") )


### geo_1112 --------------------------------

# 234
b_tre.234    <- "/Volumes/EDGE\ 2/LoVE/ReassortSubtype/b_geo/234/1112/234_h5_1112_ann.tre"
rawbeast.234 <- read.beast( b_tre.234 )
tredata.234  <-  fortify( rawbeast.234 )

g <- 
  ggtree( rawbeast.234, right = TRUE, size = 1, mrsd = "2011-12-18") + 
  aes( color = geo ) + 
  scale_color_manual( values = c("#8c564b", "#2ca02c", "#d62728", "#1f77b4", "#17becf") ) + 
  scale_fill_manual( values = c("#8c564b", "#2ca02c", "#d62728", "#1f77b4", "#17becf") ) +
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


# 232

b_tre.232    <- "/Volumes/EDGE\ 2/LoVE/ReassortSubtype/b_geo/232/1112/232_h5_1112_ann.tre"
rawbeast.232 <- read.beast( b_tre.232 )
tredata.232  <-  fortify( rawbeast.232 )

k <- 
  ggtree( rawbeast.232, right = TRUE, size = 1,mrsd = "2011-12-30") + 
  aes( color = geo ) + 
  scale_color_manual( values = c("#8c564b", "#2ca02c", "#7f7f7f", "#d62728", "#1f77b4", "#bcbd22","#17becf" ) ) + 
  scale_fill_manual( values = c("#8c564b", "#2ca02c", "#7f7f7f", "#d62728", "#1f77b4", "#bcbd22","#17becf" ) ) +
  theme_tree2( axis.text.x = element_text( size = 16  ), legend.position = c(0,0), legend.justification = c(0,0)) + 
  scale_y_continuous( expand = c(0,5) ) + 
  geom_tippoint(aes(fill = geo), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2014.5, by = 2), 
                      labels = seq(2002, 2014, by = 2) )  + 
  theme( axis.ticks  = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 18))

rectdf <- data.frame( xstart = seq( 2002, 2013, 2), 
                      xend   = seq( 2003, 2014, 2))
k + geom_rect(data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
              fill = "gray", alpha = 0.2, inherit.aes=FALSE)


# persistence 

PACT_232 <- read.table("/Volumes/EDGE\ 2/LoVE/ReassortSubtype/b_geo/out.232.stats", header = TRUE)
PACT_234 <- read.table("/Volumes/EDGE\ 2/LoVE/ReassortSubtype/b_geo/out.234.stats", header = TRUE)

# 232
ggplot(data = PACT_232[-1,], aes( x= statistic, y = mean, color = statistic)) +
  geom_point(size = 5) +
  geom_errorbar( aes( ymin = lower, ymax = upper ), width = 0.01, size = 1) +
  coord_flip() + 
  xlab("") + ylab("Persistence") +
  scale_color_manual( values = c("#8c564b", "#2ca02c", "#7f7f7f", "#d62728", "#1f77b4", "#bcbd22","#17becf" ) ) + 
  
  scale_x_discrete( labels = c("Central China", "Eastern China", "Northern China", "Southern China", 
                               "SouthWest China", "Northern Asia", "Southeast Asia") ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y     = element_blank(),
         panel.border     = element_blank(),
         axis.text.y      = element_text(size = 16),
         axis.text.x      = element_text(size = 15),
         axis.title.x     = element_text(size = 18),
         legend.position  = "none") 
  
# 234
ggplot(data = PACT_234[-1,], aes( x= statistic, y = mean, color = statistic)) +
  geom_point(size = 5) +
  geom_errorbar( aes( ymin = lower, ymax = upper ), width = 0.01, size = 1) +
  coord_flip() + 
  xlab("") + ylab("Persistence") +
  scale_color_manual( values = c("#8c564b", "#2ca02c", "#d62728", "#1f77b4", "#17becf") ) + 
  
  scale_x_discrete( labels = c("Central China", "Eastern China", "Southern China", 
                               "SouthWest China",  "Southeast Asia") ) +
  scale_y_continuous( limits = c(0,3) ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y     = element_blank(),
         panel.border     = element_blank(),
         axis.text.y      = element_text(size = 16),
         axis.text.x      = element_text(size = 15),
         axis.title.x     = element_text(size = 18),
         legend.position  = "none") 


### MCC tree 12/10 --------------------------------

# 232 ha 
mcc_tre.232  <- "ride_1205/result/com/232_h5_1205-anno.tree"
rawbeast.232 <- read.beast( mcc_tre.232 )

k2 <- 
  ggtree( rawbeast.232, right = TRUE, size = 0.4, mrsd = "2016-4-12") + 
  scale_y_continuous( expand = c(0,1) ) + 
  # geom_tippoint(aes(fill = geo), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2016.5, by = 2), 
                      labels = seq(2002, 2016, by = 2) )  + 
  theme( axis.ticks  = element_blank(), 
         axis.text.x = element_blank() )

rectdf <- data.frame( xstart = seq( 2002, 2017, 2), 
                      xend   = seq( 2003, 2017, 2))
k22 <- k2 + geom_rect(data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                      fill = "gray", alpha = 0.2, inherit.aes=FALSE)


# 234 ha
mcc_tre.234  <- "ride_1205/result/com/234_h5_1205-anno.tree"
rawbeast.234 <- read.beast( mcc_tre.234 )

k1 <- 
  ggtree( rawbeast.234, right = TRUE, size = 0.4, mrsd = "2011-12-15") + 
  theme_tree2( axis.text.x = element_text( size = 16  ) ) + 
  scale_y_continuous( expand = c(0,1) ) + 
  # geom_tippoint(aes(fill = geo), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2016.5, by = 2), 
                      labels = seq(2002, 2016, by = 2) )  + 
  theme( axis.ticks.x  = element_blank() ) + 
  
  geom_point(x = 2010.071-3.0057, y = 152, size = 2, color = "blue") + 
  geom_errorbarh( aes(y = 152, 
                      xmax = 2010.071-2.2808, 
                      xmin = 2010.071-4.1949), color = "blue", size = 1, height = 0)


rectdf <- data.frame( xstart = seq( 2002, 2017, 2), 
                      xend   = seq( 2003, 2017, 2))

k11 <- k1 + geom_rect(data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
                      fill = "gray", alpha = 0.2, inherit.aes=FALSE)


# Nx 
mcc_tre.Nx  <- "Nx_1206/results/h5_1206_S_con-anno.tree"
rawbeast.Nx <- read.beast( mcc_tre.Nx )
tredata.Nx  <- fortify( rawbeast.Nx )

ggtree( rawbeast.Nx, right = TRUE, mrsd = "2010-01-27" ) + theme_tree2() + 
  


multiplot(k11, k22, ncol = 1)


## MRCA of Nx ----------------

#h5
h5_1206_mrca0 <- read.table("Nx_1206/results/h5_1206", sep = "\t", stringsAsFactors = FALSE, 
                           row.names = NULL, header = T)
h5_1206_mrca0 = h5_1206_mrca0[,-1]

h5_1206_mrca <- data.frame( Method = c( "strict_constant", "strict_skyride", "UCLN_constant", "UCLN_skyride"), 
                            median = as.numeric(h5_1206_mrca0[5,][1:4]), 
                            hpd_l  = as.numeric( str_match( h5_1206_mrca0[8,][1:4], "\\[([0-9.]+)")[,2] ),
                            hpd_u  = as.numeric( str_match( h5_1206_mrca0[8,][1:4], "([0-9.]+)\\]")[,2] ) )

ggplot( data = h5_1206_mrca, aes( x= Method, y = 2010.071-median, color = Method )) +
  geom_point(size = 4) +
  geom_errorbar( aes( ymin = 2010.071-hpd_l, ymax = 2010.071-hpd_u), width = 0, size = 1) +
  coord_flip() + 
  xlab("") + ylab("") +
  scale_color_manual(values = c( "blue", "black","black", "black") ) + 
  scale_y_continuous(breaks = seq(2000,2008, by = 1), limits = c(2005,2008.5)) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         text = element_text(size = 20, face = "bold"), 
         legend.position = "none") 


#n5
n5_1206_mrca0 <- read.table("Nx_1206/results/n5_1206", sep = "\t", stringsAsFactors = FALSE, 
                            row.names = NULL, header = T)
n5_1206_mrca0 = n5_1206_mrca0[,-1]

n5_1206_mrca <- data.frame( Method = c( "strict_constant", "strict_skyride", "UCLN_constant", "UCLN_skyride"), 
                            median = as.numeric(n5_1206_mrca0[5,][1:4]), 
                            hpd_l  = as.numeric( str_match( n5_1206_mrca0[8,][1:4], "\\[([0-9.]+)")[,2] ),
                            hpd_u  = as.numeric( str_match( n5_1206_mrca0[8,][1:4], "([0-9.]+)\\]")[,2] ) )

ggplot( data = n5_1206_mrca, aes( x= Method, y = 2010.071-median, color = Method )) +
  geom_point(size = 4) +
  geom_errorbar( aes( ymin = 2010.071-hpd_l, ymax = 2010.071-hpd_u), width = 0, size = 1) +
  coord_flip() + 
  xlab("") + ylab("") +
  scale_color_manual(values = c( "blue", "black","black", "black") ) + 
  # scale_y_continuous(breaks = seq(2000,2008, by = 1) ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         text = element_text(size = 20, face = "bold"), 
         legend.position = "none") 
