library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)

source("~/Packaging_Type/_R/Function.R")
setwd("/Volumes/EDGE 2/LoVE/ReassortSubtype/Ne/")


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

  #geom_rect(data = rects, 
  #          aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), 
  #          alpha = 0.4) 
  

## MRCA of Nx virus (rate_0818) ----------------

rate_0818_log <- data.frame( method = c("strict_exp", "strict_skyride", "relaxed_exp", "relaxed_skyride"),
                             median = c(3.0269, 3.0468, 2.9611, 3.0233),
                             hpd_l  = c(2.2218, 2.2994, 1.7829, 1.9933),
                             hpd_u  = c(4.5126, 4.2755, 10.0257, 5.2335), stringsAsFactors = FALSE)

ggplot(data = rate_0818_log, aes( x= method, y = 2010.071-median )) +
  geom_point() +
  geom_errorbar( aes( ymin = 2010.071-hpd_l, ymax = 2010.071-hpd_u), width = 0.1) +
  coord_flip() + 
  xlab("") + ylab("") +
  scale_y_continuous(breaks = seq(2000,2008, by = 1)) +
  theme_bw() + ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), panel.grid.major.y = element_blank() ) 
  

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
ha_meanRate <- as.data.frame( t(read.table("./skygrid_0810/result_0810/ha_meanRate.csv", 
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
na_meanRate <- as.data.frame( t(read.table("./skygrid_0810/result_0810/na_meanRate.csv", 
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




