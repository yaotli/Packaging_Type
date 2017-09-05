# This file aims to construct function - subfastaSeq

library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

seq_dir <- paste0("./Clade/sampled_clade_tree_CNHK/sampled_fas/", 
                  list.files( "./Clade/sampled_clade_tree_CNHK/sampled_fas/" ) )

## host ----------------

seq_name_combined = 
as.vector( do.call(c, sapply( seq_dir, function(x) attributes( read.fasta(x) )$names ) ) )

nohuman <- 
  "avian|bird|wildbird|chicken|duck|dove|pigeon|mallard|goose|environment|water|teal|hawk|magpie|munia|myna|kestrel|peregrine|crow|sparrow|mesia|gull|egret|swan|shrike|buzzard|heron|Ph|quail|pheasant|grebe|starling|swallow|white_eye|swine|tiger"

humanID <- grep( nohuman, seq_name_combined, ignore.case = TRUE, value = TRUE, invert = TRUE)

## geo ----------------

temp1.geo <- 
grep("[A-Z]+", 
     sapply( strsplit( humanID, split = "_"), function(x) x[[2]] ), value = TRUE)

temp2.geo <- 
unique( grep("[A-Z]+", sapply( strsplit( animalID, split = "_"), function(x) x[[3]] ), value = TRUE) )

geo_combined <- unique( c( hu.geo, a.geo ) )[-14]
geo_combined[1]  = "Hong_Kong"
geo_combined[28] = "Eastern_China"
geo_combined[30] = "North_China"

geo_combined <- c(geo_combined, "Heilongjiang", "GD", "Qinghai")

geo_term     <- paste0(geo_combined, collapse = "|")

# check
# grep(geo_term, seq_name_combined, invert = TRUE, value = TRUE, ignore.case = TRUE)

geo.N <- c("Shanxi", "Hebei", "Beijing", "Jilin", "Sheny", 
            "Liaoning", "Heilongjiang", "North_China")
geo_N <- paste0( geo.N, collapse = "|")

geo.C <- c("Hunan", "Hubei", "Henan")
geo_C <- paste0( geo.C, collapse = "|")

geo.E <- c("Shandong", "Jiangsu", "Huadong", "Eastern_China", "Fujian", 
          "Anhui", "Shanghai", "Zhejiang", "Jiangxi", "Nanchang")
geo_E <- paste0( geo.E, collapse = "|")

geo.S <- c("Hong_Kong", "Shantou", "Guangdong", "GD", "ST", 
          "Shenzhen", "Guangzhou")
geo_S <- paste0( geo.S, collapse = "|")

geo.SW <- c("Guizhou", "Guangxi", "Yunnan", "Guiyang", "Tibet", 
            "Sichuan", "Chongqing")
geo_SW <- paste0( geo.SW, collapse = "|")

geo.NW <- c("Ningxia", "Xinjiang", "Gansu", "Qinghai", "Shaanxi")
geo_NW <- paste0( geo.NW, collapse = "|")







