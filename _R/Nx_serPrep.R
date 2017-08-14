library(seqinr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

h5selected <- "./Nx/h5n5.txt"
HAseqfile  <- "./align_trim/trim_h5_pool_lth1683.fasta"

## extract sequences ----------------

subtreseq(seq_filedir = HAseqfile, list_filedir = h5selected)

# check by raxml tree
# ./raxml_AVX2 -f a -p 666 -s s_h5_n5_2344.fasta -x 616 -#autoMRE -m GTRGAMMAI --HKY85 -n h5n5
# 
# check substitution model
# java -jar /Volumes/EDGE\ 2/Apps/jmodeltest-2.1.10/jModelTest.jar -d /Users/yaosmacbook/Desktop/data_souce/Nx/s_h5_n5_2344.fas -f -i -g 4 -s 5 -BIC -a -S NNI

