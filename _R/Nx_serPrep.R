library(seqinr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

h5_raw_ncbi <- "./raw/H5_N_7412_20170815.fasta"

### seq name cleaning --------------------------------

# readin
ha_n_seq    <- fastaEx(h5_raw_ncbi)$seq
ha_n_id     <- fastaEx(h5_raw_ncbi)$id

# deal with _A/
ha_n_id[ grep( "_A/", ha_n_id, invert = TRUE) ] <-
  str_replace( string      = ha_n_id[grep( "_A/", ha_n_id, invert = TRUE) ], 
               pattern     = "([A-Z]{1,2})([0-9]{5,6})_", 
               replacement = paste0( 
                 str_match(ha_n_id[grep( "_A/", ha_n_id, invert = TRUE) ], 
                           "([A-Z]{1,2})([0-9]{5,6})_")[,1], "A/") )

ha_n_id <- gsub("A/", "", ha_n_id)

# export
write.fasta(ha_n_seq, ha_n_id, "h5_n_nonpaired.fasta")

# cleanID
cleanID("./Nx/h5_n_nonpaired.fasta")

# curate
# delete = 571 | remain = 6841
curateSeq(maxamb = 150, minseq = 1500, mode = 1, filedir = "./Nx/cleanID_h5_n_nonpaired.fasta")

# align the seq by MAFFT
# SHELL:
# 
# cd ~/Desktop/data_souce/Nx
# mafft ./curateSeq-1_cleanID_h5_pool.fasta > ./align_trim/align_h5_pool.fasta
#
# manually trim in BioEdit (remove stop codon & ~, lth = 1683)


## extract sequences ----------------


# check by raxml tree
# ./raxml_AVX2 -f a -p 666 -s h5_n5.fasta.fasta -x 616 -#autoMRE -m GTRGAMMAI --HKY85 -n h5n5
# 
# check substitution model
# java -jar /Volumes/EDGE\ 2/Apps/jmodeltest-2.1.10/jModelTest.jar -d /Users/yaosmacbook/Desktop/data_souce/Nx/h5_n5_phy -f -i -g 4 -s 5 -BIC -a -S SPR

