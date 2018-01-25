library(seqinr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/Nov2017_hana/")
source("~/Packaging_Type/_R/Function.R")

### based on previous selection H5N5 (n = 8) --------------------------------

oldh5id    <- tagExtra( "~/Desktop/data_souce/Nx/c2344_21_e1206.tre" )$id[ which( tagExtra( "~/Desktop/data_souce/Nx/c2344_21_e1206.tre" )$tag == "ff0000" ) ]
oldh5id_ac <- str_match( oldh5id, "^[A-Z0-9]+" )[,1]

leafEx( "c234/raw/pH5_c234_2429.fasta", oldh5id )

n5id_ac <- read.csv("R_process/pH5NA.csv", header = TRUE, stringsAsFactors = FALSE)$ac.na[ match( oldh5id_ac, read.csv("R_process/pH5NA.csv", header = TRUE, stringsAsFactors = FALSE)$ac.ha ) ]
subfastaSeq( AC = TRUE, AC_list = n5id_ac, filedir = "./R_process/pNA_7138.fasta")


### all available H5Nx before 2013/2014 --------------------------------

rmDup( "c234/raw/pH5_c234_2429.fasta", year = c(2001, 2014), sero = "H5N1", 
       sero.invert = TRUE, rmdup = TRUE)