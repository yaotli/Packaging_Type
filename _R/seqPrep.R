# prepare sequences for major HA-NA comparative study
# NOTE: 1. use underline (_) replace blank in ID from NCBI
#       2. check >_A (replace with >A_) in ID from GISAID

library(seqinr)
library(stringr)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

# gisaid
file_ha_g <- "./pH5_G_1856_20170608.fasta" 
file_na_g <- "./pNA_G_1855_20170608.fasta"
csv_g     <- "./pH5NA_G_1855_20170608.csv"


# ncbi
file_ha_n <- "./pH5_N_5483_20170608.fasta" 
file_na_n <- "./pNA_N_5338_20170608.fasta"
csv_n_ha  <- "./pH5_N_5483_20170608.csv"
csv_n_na  <- "./pNA_N_5338_20170608.csv"


### gisaid --------------------------------

## remove sequence duplicated isolatedID ----------------

# HA

ha_g_seq    <- keepLongSeq( fastaEx(file_ha_g)$seq, 
                         fastaEx(file_ha_g)$id )$seq

ha_g_id     <- keepLongSeq( fastaEx(file_ha_g)$seq, 
                         fastaEx(file_ha_g)$id)$id

# deal with assession number 

ha_g_id_ID  <- gsub("_ISL_", "", 
                    str_match(ha_g_id, "EPI_ISL_([0-9]+)" )[, 1])

ha_g_id_epi <- str_replace(gsub("_EPI_ISL_([0-9]+)", "", ha_g_id), 
                           "A/", paste0(ha_g_id_ID, "_") ) 


# NA 

nu_g_seq    <- fastaEx(file_na_g)$seq
nu_g_id     <- fastaEx(file_na_g)$id

# deal with assession number 

nu_g_id_ID  <- gsub("_ISL_", "", 
                   str_match(nu_g_id, "EPI_ISL_([0-9]+)" )[, 1])

nu_g_id_epi <- str_replace(gsub("_EPI_ISL_([0-9]+)", "", nu_g_id), 
                           "A/", paste0(nu_g_id_ID, "_") ) 



## prepare the datasheet

### ncbi --------------------------------

# HA 

ha_n_seq    <- fastaEx(file_ha_n)$seq
ha_n_id     <- fastaEx(file_ha_n)$id

# id start without "A/"
ha_n_id[ which( startsWith(prefix = "A/", x = ha_n_id) == FALSE) ] <- 
  paste0("A/", ha_n_id[ which( startsWith(prefix = "A/", x = ha_n_id) == FALSE) ])

ha_n_id_epi <- 
  str_replace(string      = gsub("(_[A-Z]{1,2})([0-9]{5,6})", "", ha_n_id), 
              pattern     = "A/", 
              replacement = paste0( str_match( ha_n_id, "([A-Z]{1,2})([0-9]{5,6})" )[, 1], "_") )
  
# NA 

nu_n_seq <- fastaEx(file_na_n)$seq
nu_n_id  <- fastaEx(file_na_n)$id

nu_n_id[ which( startsWith(prefix = "A/", x = nu_n_id) == FALSE) ] <- 
  paste0("A/", nu_n_id[ which( startsWith(prefix = "A/", x = nu_n_id) == FALSE) ])

nu_n_id_epi <- 
  str_replace(string      = gsub("(_[A-Z]{1,2})([0-9]{5,6})", "", nu_n_id), 
              pattern     = "A/", 
              replacement = paste0( str_match( nu_n_id, "([A-Z]{1,2})([0-9]{5,6})" )[, 1], "_") )


## read-in .csv files

## keep paired HA-NA and prepare datasheet










