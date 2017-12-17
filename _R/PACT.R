library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)
library(tidyr)


source("~/Packaging_Type/_R/Function.R")
setwd("/Volumes/EDGE2/LoVE/ReassortSubtype/BEAST/")


a = read.table( "~/PACT-master/out.232.skylines", header = T, stringsAsFactors = F )

ggplot( data = a, aes( x = time, y = mean )) + 
  geom_area( aes(color = statistic, fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous(breaks = seq(2004,2016, by = 2), limits = c(2004, 2016)) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y =  element_blank(), 
         panel.border  = element_blank(), 
         text = element_text(size = 15, face = "bold"), 
         legend.position = "none") 


b = read.table( "~/PACT-master/out.234.skylines", header = T, stringsAsFactors = F )

ggplot( data = b, aes( x = time, y = mean )) + 
  geom_area( aes(color = statistic, fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous(breaks = seq(2004,2016, by = 2), limits = c(2004, 2016)) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y =  element_blank(), 
         panel.border  = element_blank(), 
         text = element_text(size = 15, face = "bold"), 
         legend.position = "none") 
