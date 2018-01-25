library(seqinr)
library(dplyr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/Nov2017_hana/")
source("~/Packaging_Type/_R/Function.R")

subtreseq( list_filedir = "./c2344/ls_c2344_1385.txt", 
           seq_filedir  = "./tree/trim_pH5_7326_3.fasta" )


# rmDup

rmDup( fasfile = "./c2344/pH5_c2344_1385.fasta", rmdup = TRUE )

