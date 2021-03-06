# data preparation: seqPrep.R

library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

pool_csv <- "./pool_df.csv"
gsgdtree <- "./Tree/h5_GsGD"
n1tree   <- "./Tree/N1_pool"

h5_GsGD_table_dir <- "./h5_GsGD_table.csv"
n1_table_dir      <- "./n1_table.csv"   


### read-in table --------------------------------

h5_GsGD_table <- read.csv( h5_GsGD_table_dir, header = TRUE, stringsAsFactors = FALSE)
n1_table      <- read.csv( n1_table_dir, header = TRUE, stringsAsFactors = FALSE)


# stratified by geo

h5_GsGD_table_cnhk <- h5_GsGD_table[which(h5_GsGD_table$geo == "China" | h5_GsGD_table$geo == "Hong_Kong" ), ]
n1_table_cnhk      <- n1_table[which(n1_table$geo == "China" | n1_table$geo == "Hong_Kong" ), ]



### Figuree --------------------------------

ggplot( data = h5_GsGD_table, 
        aes(x = time, y = distance) ) + 
  geom_point(color = "gray") +
  geom_point( data = h5_GsGD_table[which(h5_GsGD_table$subtype == "N1"), ],
              aes(x = time, y = distance) ) +
  theme_bw() + 
  facet_wrap(~ geo)


## clade234_HA_China ----------------

clade234_HA_cnhk <- 
  
  ggplot( data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade234_h5_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$clade234 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$s_clade234_h5_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  # geom_vline(xintercept = 2008, size = 1, color = "black") + 
  ggtitle("clade234_HA_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  xlab("") + ylab("Root-to-tip divergence")



## clade7_HA_China ----------------

clade7_HA_cnhk <- 
  
  ggplot( data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade7_h5_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$clade7 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$s_clade7_h5_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  # geom_vline(xintercept = 2006, size = 1, color = "black") + 
  ggtitle("clade7_HA_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  xlab("") + ylab("Root-to-tip divergence")



## clade232_HA_China ----------------

clade232_HA_cnhk <- 
  
  ggplot( data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade232_h5_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$clade232 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$s_clade232_h5_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  # geom_vline(xintercept = 2012, size = 1, color = "black") + 
  ggtitle("clade232_HA_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  xlab("") + ylab("Root-to-tip divergence") 
  # + geom_text(aes(label=name), size = 1)
  



## clade234_N1_China ----------------

clade234_N1_cnhk <- 
  
  ggplot( data = n1_table_cnhk[ which(n1_table_cnhk$s_clade234_n1_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = n1_table_cnhk[ which(n1_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$clade234 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$s_clade234_n1_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  # geom_vline(xintercept = 2008, size = 1, color = "black") + 
  ggtitle("clade234_N1_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  scale_y_continuous(limits = c(0.05, 0.2) ) + 
  xlab("") + ylab("Root-to-tip divergence")





## clade7_NA_China ----------------

clade7_N1_cnhk <- 
  
  ggplot( data = n1_table_cnhk[ which(n1_table_cnhk$s_clade7_n1_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = n1_table_cnhk[ which(n1_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$clade7 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$s_clade7_n1_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  # geom_vline(xintercept = 2006, size = 1, color = "black") + 
  ggtitle("clade7_N1_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  scale_y_continuous(limits = c(0.05, 0.2) ) + 
  xlab("") + ylab("Root-to-tip divergence")



## clade232_NA_China ----------------

clade232_N1_cnhk <- 
  
  ggplot( data = n1_table_cnhk[ which(n1_table_cnhk$s_clade232_n1_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = n1_table_cnhk[ which(n1_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$clade232 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$s_clade232_n1_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  # geom_vline(xintercept = 2012, size = 1, color = "black") + 
  ggtitle("clade232_N1_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  scale_y_continuous(limits = c(0.05, 0.2) ) + 
  xlab("") + ylab("Root-to-tip divergence")


## combine ----------------

multiplot(clade234_HA_cnhk, clade7_HA_cnhk, clade232_HA_cnhk, 
          clade234_N1_cnhk, clade7_N1_cnhk, clade232_N1_cnhk, ncol = 2)


### rate comparison --------------------------------

N5_8 <- c("CY091643_duck_Guangdong_wy24_2008_H5N5_2008.430", 
          "CY091627_duck_Guangdong_wy11_2008_H5N5_2008.430",
          "GU727661_duck_Eastern_China_031_2009_H5N5_2009.038",
          "GU727653_duck_Eastern_China_008_2008_H5N5_2008.956",
          "JQ973670_goose_Shandong_k1204_2009_H5N5_2009.956",
          "CY091635_duck_Guangdong_wy19_2008_H5N5_2008.430",
          "JQ973678_quail_Jiangsu_k0104_2010_H5N5_2010.019",
          "JQ973686_goose_Guangdong_k0103_2010_H5N5_2010.071" )


h5_GsGD_table_cnhk_N5_8 <- h5_GsGD_table_cnhk[match(N5_8, h5_GsGD_table_cnhk$name), ]


ggplot() + 
  theme_bw() + 
  xlab("") + ylab("Root-to-tip divergence (sub./site)") + 
  
  theme(panel.grid.minor = element_blank(), 
       panel.grid.major = element_blank(), 
       axis.ticks = element_blank(), 
       panel.border = element_blank(),
       axis.text.x = element_text(size = 18), 
       axis.title.y = element_text(size = 18)) +
  
  scale_x_continuous( limits = c(2000, 2016), 
                      breaks = seq(2000, 2016, by = 2), 
                      labels = seq(2000, 2016, by = 2)) +
  
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, size = 3, stroke = 0) + 
  
  stat_smooth(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
              aes(x = time, y = distance), method = "glm", 
              colour = "gray", size = 2, alpha = 0.8) +

  
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade234_h5_CNHK == 1), ], 
             aes(x = time, y = distance),
             color = "#d62728", alpha = 0.5, size = 3, stroke = 0) + 
  
  stat_smooth(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade234_h5_CNHK == 1), ], 
              aes(x = time, y = distance), method = "glm", 
              colour = "#d62728", size = 2, se = FALSE, alpha = 0.8) +

  
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade232_h5_CNHK == 1), ], 
             aes(x = time, y = distance),
             color = "#2ca02c", alpha = 0.5, size = 3, stroke = 0) + 
    
  stat_smooth(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade232_h5_CNHK == 1 & h5_GsGD_table_cnhk$time < 2014), ], 
              aes(x = time, y = distance), method = "glm", 
              colour = "#2ca02c", size = 2, se = FALSE, alpha = 0.8) +
  

  geom_point(data = h5_GsGD_table_cnhk_N5_8, 
             aes(x = time, y = distance),
             color = "#1f77b4", size = 3, stroke = 0) 
  
  #stat_smooth(data = h5_GsGD_table_cnhk_N5_8, 
  #            aes(x = time, y = distance), method = "glm", 
  #           colour = "#1f77b4", size = 2, se = FALSE)   


  
  
  









