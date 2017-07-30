# prepare sequences for major HA-NA comparative study
# NOTE: 1. use underline (_) replace blank in ID from NCBI
#       2. check >_A (replace with >A_) in ID from GISAID

library(seqinr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

# gisaid
file_ha_g <- "./raw/pH5_G_1856_20170608.fasta" 
file_na_g <- "./raw/pNA_G_1855_20170608.fasta"
csv_g     <- "./raw/pH5NA_G_1855_20170608.csv"


# ncbi
file_ha_n <- "./raw/pH5_N_5592_20170612.fasta" 
file_na_n <- "./raw/pNA_N_5560_20170612.fasta"
csv_n_ha  <- "./raw/pH5_N_5592_20170612.csv"
csv_n_na  <- "./raw/pNA_N_5560_20170612.csv"


### gisaid --------------------------------

## HA ----------------

# HA file has data with duplicated assession no. (identical id)
# check by - 
# grep("EPI_ISL_169351", fastaEx(file_ha_g)$id)
# make HA and NA .fas have same # of sequences data

ha_g_seq    <- keepLongSeq( fastaEx(file_ha_g)$seq, 
                         fastaEx(file_ha_g)$id )$seq

ha_g_id     <- keepLongSeq( fastaEx(file_ha_g)$seq, 
                         fastaEx(file_ha_g)$id)$id

# move accession no. from back to front

ha_g_ac     <- gsub("_ISL_", "", 
                    str_match(ha_g_id, "EPI_ISL_([0-9]+)" )[, 1])

ha_g_id_epi <- str_replace( string      = gsub( "_EPI_ISL_([0-9]+)", "", ha_g_id), 
                            pattern     = "A/", 
                            replacement = paste0(ha_g_ac, "_") ) 

ha_g_id_ac0 <- gsub("_EPI_ISL_([0-9]+)", "", ha_g_id)


# duplicated name and remove name with "RG"
# noDupName_g, n = 1849

noDupName_g <- keepLongSeq(ha_g_seq, ha_g_id_ac0, showRemain = TRUE)

ha_g_seq    <- ha_g_seq[ noDupName_g ]
ha_g_id_epi <- ha_g_id_epi[  noDupName_g ]

# rg_ha_g, n = 2, resulting file, n = 1847

rg_ha_g     <- grep(pattern = "RG", ha_g_id_epi)

ha_g_seq    <- ha_g_seq[ seq(1, length(ha_g_seq))[- rg_ha_g ] ]
ha_g_id_epi <- ha_g_id_epi[  seq(1, length(ha_g_id_epi))[- rg_ha_g ] ]


## NA ----------------

nu_g_seq    <- fastaEx(file_na_g)$seq
nu_g_id     <- fastaEx(file_na_g)$id

# move accession no. from back to front

nu_g_ac     <- gsub("_ISL_", "", 
                   str_match(nu_g_id, "EPI_ISL_([0-9]+)" )[, 1])

nu_g_id_epi <- str_replace( string = gsub("_EPI_ISL_([0-9]+)", "", nu_g_id), 
                            pattern = "A/", 
                            replacement = paste0(nu_g_ac, "_") ) 

nu_g_id_ac0 <- gsub("_EPI_ISL_([0-9]+)", "", nu_g_id)


# duplicated name and remove name with "RG"
# resulting in 1847 data as HA

noDupName_g_nu <- match( ha_g_ac[noDupName_g], nu_g_ac )

nu_g_seq    <- nu_g_seq[ noDupName_g_nu ]
nu_g_id_epi <- nu_g_id_epi[  noDupName_g_nu ]


rg_nu_g     <- grep(pattern = "RG", nu_g_id_epi)

nu_g_seq    <- nu_g_seq[ seq(1, length(nu_g_seq))[- rg_nu_g  ] ]
nu_g_id_epi <- nu_g_id_epi[  seq(1, length(nu_g_id_epi))[- rg_nu_g  ] ]


## prepare the datasheet ----------------

sheet_g_file      <- read.csv(csv_g, header = TRUE, stringsAsFactors = FALSE)

sheet_g_assession <- gsub("_ISL_", "", sheet_g_file$Isolate_Id)
sheet_g_subtype   <- gsub("A / H5", "", sheet_g_file$Subtype)


# geo

sheet_g_geo       <- gsub(pattern = " ", "_", sheet_g_file$Location)
sheet_g_geo       <- gsub("_$", "", 
                          str_match(sheet_g_geo, "([A-Za-z_]+)_/_([A-Za-z_]+)")[,3])

# combine 

sheet_g           <- data.frame(ac      = sheet_g_assession, 
                                name    = sheet_g_file$Isolate_Name, 
                                geo     = sheet_g_geo, 
                                subtype = sheet_g_subtype,
                                segment = 46,
                                pac     = sheet_g_assession, 
                                stringsAsFactors = FALSE)


# check accession no. (ac)
# TRUE %in% is.na( match( str_match(ha_g_id_epi, "EPI([0-9]+)")[,1], sheet_g$ac ))
# TRUE %in% is.na( match( str_match(nu_g_id_epi, "EPI([0-9]+)")[,1], sheet_g$ac ))
# TRUE %in% is.na( match( str_match(nu_g_id_epi, "EPI([0-9]+)")[,1], str_match(ha_g_id_epi, "EPI([0-9]+)")[,1] ))

### ncbi --------------------------------

# readin

ha_n_seq    <- fastaEx(file_ha_n)$seq
ha_n_id     <- fastaEx(file_ha_n)$id

nu_n_seq    <- fastaEx(file_na_n)$seq
nu_n_id     <- fastaEx(file_na_n)$id

# id start without "A/"

ha_n_id[ grep( "_A/", ha_n_id, invert = TRUE) ] <-
  str_replace( string      = ha_n_id[grep( "_A/", ha_n_id, invert = TRUE) ], 
               pattern     = "([A-Z]{1,2})([0-9]{5,6})_", 
               replacement = paste0( 
                 str_match(ha_n_id[grep( "_A/", ha_n_id, invert = TRUE) ], 
                           "([A-Z]{1,2})([0-9]{5,6})_")[,1], "A/") )
  
nu_n_id[grep( "_A/", nu_n_id, invert = TRUE) ] <-
  str_replace( string      = nu_n_id[grep( "_A/", nu_n_id, invert = TRUE) ], 
               pattern     = "([A-Z]{1,2})([0-9]{5,6})_", 
               replacement = paste0( 
                 str_match(nu_n_id[grep( "_A/", nu_n_id, invert = TRUE) ], 
                           "([A-Z]{1,2})([0-9]{5,6})_")[,1], "A/") )

sheet_ha_n_file   <- read.csv(csv_n_ha, header = FALSE, stringsAsFactors = FALSE)
sheet_nu_n_file   <- read.csv(csv_n_na, header = FALSE, stringsAsFactors = FALSE)


# some virus name does not appear in sheet_nu_n_file$V9
# check - 
# sheet_ha_n_file$V9[ which( sheet_ha_n_file$V9 %in% "Influenza A virus" == TRUE) ]
# ha_n_id[ which( sheet_ha_n_file$V9 %in% "Influenza A virus" == TRUE)  ]

lossname_ha <- which( sheet_ha_n_file$V9 %in% "Influenza A virus" == TRUE)

sheet_ha_n_file$V9[ lossname_ha ] = 
  str_match( ha_n_id[ lossname_ha ], "(A/[A-Za-z0-9/_-]+)(_H5N)")[,2]


lossname_nu <- which( sheet_nu_n_file$V9 %in% "Influenza A virus" == TRUE)

sheet_nu_n_file$V9[ lossname_nu ] = 
  str_match( nu_n_id[ lossname_nu ], "(A/[A-Za-z0-9/_-]+)(_H5N)")[,2]


## find the discrepancy in HA and NA ----------------

# HA 

# the one that does not appear unique (name)
# find by: setdiff( unique(sheet_ha_n_file$V9), unique(sheet_nu_n_file$V9) )
# result : "Influenza A virus (A/yellow-billed teal/Chile/14/2014(H5Nx))"

# keep unique name (from .csv file)

ha_n_select_remain <- keepLongSeq( ha_n_seq, sheet_ha_n_file$V9, showRemain = TRUE)

# eliminate the no.4926 (A/yellow-billed teal/Chile/14/2014)

HANAdis            <- match( setdiff( unique(sheet_ha_n_file$V9), unique(sheet_nu_n_file$V9) ), 
                             sheet_ha_n_file$V9)

ha_n_select_remain <- ha_n_select_remain[-which(ha_n_select_remain == HANAdis)]
  
# the supposed non-redundant seq

ha_n_select_id     <- ha_n_id[ha_n_select_remain]
ha_n_select_seq    <- ha_n_seq[ha_n_select_remain]

ha_n_select_id_ac0 <- sub("([A-Z]{1,2})([0-9]{5,6})_", "", ha_n_select_id)
ha_n_select_id_min <- gsub("_H5N[0-9]_([0-9]{4})-([0-9-]+)", "", ha_n_select_id_ac0)

ha_n_select_ac     <- str_match( ha_n_select_id, "([A-Z]{1,2})([0-9]{5,6})")[,1]

  
# NA

# keep unique name (from .csv file)

nu_n_select_remain <- keepLongSeq( nu_n_seq, sheet_nu_n_file$V9, showRemain = TRUE )

nu_n_select_id     <- nu_n_id[nu_n_select_remain]
nu_n_select_seq    <- nu_n_seq[nu_n_select_remain]

nu_n_select_id_ac0 <- sub( "([A-Z]{1,2})([0-9]{5,6})_", "", nu_n_select_id)
nu_n_select_id_min <- gsub("_H5N[0-9]_([0-9]{4})-([0-9-]+)", "", nu_n_select_id_ac0)

nu_n_select_ac     <- str_match( nu_n_select_id, "([A-Z]{1,2})([0-9]{5,6})")[,1]


# check: setdiff(ha_n_select_id_ac0, nu_n_select_id_ac0)
# check: setdiff(ha_n_select_id_min, nu_n_select_id_min)
# find : 78 


## unify the HA & NA info ----------------


  dismatch     <- setdiff( ha_n_select_id_ac0, nu_n_select_id_ac0 )
  dismatch_na  <- setdiff( nu_n_select_id_ac0, ha_n_select_id_ac0 )
  dismatch_min <- gsub("_H5N[0-9]_([0-9]{4})-([0-9-]+)", "", dismatch)
  
  # View(cbind(sort(dismatch), sort(dismatch_na), sort(dismatch_min) ) )  

for (i in 1: length(  setdiff(ha_n_select_id_ac0, nu_n_select_id_ac0)  ))
{
  
  ha_i         <- grep( dismatch_min[i], ha_n_select_id_ac0 )
  nu_i         <- grep( dismatch_min[i], nu_n_select_id_ac0 )
  
  appends      <- c( str_match(ha_n_select_id_ac0[ha_i], "_H5N[0-9]_([0-9]{4})-([0-9-]+)")[,1], 
                     str_match(nu_n_select_id_ac0[nu_i], "_H5N[0-9]_([0-9]{4})-([0-9-]+)")[,1] )

  append.i     <- which.max( nchar( appends ) )
    
  append       <- appends[ append.i]
  
  ha_n_select_id[ha_i] <- gsub(pattern     = "_H5N[0-9]_([0-9]{4})-([0-9-]+)", 
                               replacement = append, 
                               x           = ha_n_select_id[ha_i] )
  
  nu_n_select_id[nu_i] <- gsub(pattern     = "_H5N[0-9]_([0-9]{4})-([0-9-]+)", 
                               replacement = append, 
                               x           = nu_n_select_id[nu_i] )
  
  print( c( ha_n_select_id[ha_i],  nu_n_select_id[nu_i] ) )
}


# remove A/
ha_n_select_id <- gsub("A/", "", ha_n_select_id)
nu_n_select_id <- gsub("A/", "", nu_n_select_id)

# both HA, NA sequence, n = 5426


## prepare the datasheet ----------------

sheet_n_raw = rbind( sheet_ha_n_file[, -c(10:15)], sheet_nu_n_file[, -c(10:15)] )

sheet_n_assession <- sheet_n_raw$V1
sheet_n_subtype   <- gsub("H5", "", sheet_n_raw$V5)

sheet_n_geo       <- gsub(" ", "_", sheet_n_raw$V6)
sheet_n_geo[which(sheet_n_geo == "Cote_d'Ivoire")] = "Cote_dIvoire"

sheet_n_name      <- sub("Influenza A virus \\(|", "", 
                         sub("\\)", "", sheet_n_raw$V9))

# ncbi specific

sheet_n_segment   <- as.numeric( gsub(" \\(([A-Z]{2})\\)", "", sheet_n_raw$V4) )
  
sheet_n_NUac      <- 
  nu_n_select_ac[ match( sub( "([A-Z]{1,2})([0-9]{5,6})_", "", ha_n_select_id), 
                         sub( "([A-Z]{1,2})([0-9]{5,6})_", "", nu_n_select_id) ) ]

sheet_n_HAac      <- 
  ha_n_select_ac[ match( sub( "([A-Z]{1,2})([0-9]{5,6})_", "", nu_n_select_id), 
                         sub( "([A-Z]{1,2})([0-9]{5,6})_", "", ha_n_select_id) ) ]

# pair ac

sheet_n_pac = c()
notmatch    = c()

for( k in 1: length( sheet_n_assession) )
{
  if( sheet_n_assession[k] %in% c(nu_n_select_ac, ha_n_select_ac) == TRUE )
  {
    sheet_n_pac[ length(sheet_n_pac) + 1] <- 
      c(sheet_n_NUac, sheet_n_HAac)[
        match( sheet_n_assession[k], c(ha_n_select_ac, nu_n_select_ac) ) ]
      
  }else
  {
    sheet_n_pac[ length(sheet_n_pac) + 1] <- "NA"
    notmatch[ length(notmatch) + 1 ] = k
    
    }
}



# combine 

sheet_n           <- data.frame(ac      = sheet_n_assession, 
                                name    = sheet_n_name, 
                                geo     = sheet_n_geo, 
                                subtype = sheet_n_subtype,
                                segment = sheet_n_segment, 
                                pac     = sheet_n_pac, 
                                stringsAsFactors = FALSE)

# check accession no. (ac)
# TRUE %in% is.na( match( str_match( ha_n_select_ac, "([A-Z]{1,2})([0-9]{5,6})")[,1], sheet_n$ac))
# TRUE %in% is.na( match( str_match(nu_n_select_ac, "([A-Z]{1,2})([0-9]{5,6})")[,1], sheet_n$ac))


### combine files --------------------------------

# 7265 
ha_pool_seq <- c( ha_n_select_seq, ha_g_seq)
ha_pool_id  <- c( ha_n_select_id, ha_g_id_epi)

nu_pool_seq <- c( nu_n_select_seq, nu_g_seq )
nu_pool_id  <- c( nu_n_select_id, nu_g_id_epi)

pool_df     <- rbind(sheet_n, sheet_g)

# export

write.fasta(ha_pool_seq, ha_pool_id, "h5_pool.fasta")
write.fasta(nu_pool_seq, nu_pool_id, "nu_pool.fasta")
write.csv(pool_df, "pool_df.csv")


### sequence cleaning and curate --------------------------------

pool_h5_fasta = "./h5_pool.fasta"
pool_nu_fasta = "./nu_pool.fasta"
pool_csv      = "./pool_df.csv"


## HA ----------------

# cleanID
cleanID(pool_h5_fasta)

# curate file
# 1. maxamb = 150 | 2. minseq = 1500

# delete = 291, remain = 6982

curateSeq(maxamb = 150, minseq = 1500, mode = 1, filedir = "./cleanID_h5_pool.fasta")

# align the seq by MAFFT
# SHELL:
# 
# cd ~/Desktop/data_souce/
# mafft ./curateSeq-1_cleanID_h5_pool.fasta > ./align_trim/align_h5_pool.fasta
#

# manually trim in BioEdit (remove stop codon, lth = 1683)
# FastTree 
#
# ~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./align_trim/trim_h5_pool_lth1683.fasta> ./tree/h5_pool.tre
#
# mannually extract GsGD lineage names (copy taxa: sub_6072)

subtreseq(seq_filedir = "./align_trim/trim_h5_pool_lth1683.fasta", 
          list_filedir = "./Tree/sub_6072.txt")


# re-fasttree
# 
# ~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./tree/sub_6072.fasta> ./tree/h5_GsGD.tre
# ./raxml -f a -p 123 -s sub_6072.fasta -x 616 -#autoMRE -m GTRGAMMAI -n GsGDra

## NA ----------------

# cleanID
cleanID(pool_nu_fasta)

# curate file
# 1. maxamb = 120 | 2. minseq = 1200

# delete: 337, remain: 6936
curateSeq(maxamb = 120, minseq = 1200, mode = 1, filedir = "./cleanID_nu_pool.fasta")

## select N1 virus
# n = 4542

subfastaSeq(subtype = "H5N1", filedir = "./curateSeq-1_cleanID_nu_pool.fasta")

# align by MAFFT 
# 
# mafft ./curateSeq-1_cleanID_nu_pool_H5N1-1000-3000.fasta > ./align_trim/align_nu_pool_H5N1.fasta
#
# trim firstly by trimtool
# manually edit in BioEdit and remove stop codon (lth = 1347)

trimtool(propblank = 0.9, filedir = "./align_trim/align_nu_pool_H5N1.fasta")

# Fasttree 
# ~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./align_trim/trim_nu_pool_H5N1_lth1347.fasta> ./tree/N1_pool.tre


### Clade labeling  --------------------------------


h5_GsGDfile <- read.tree("./Tree/h5_GsGD")


## dirty match for lineage annotation ----------------

# refname derived from the smallref.fasta

refname       <- read.table("./refSeq/ref_name", stringsAsFactors = FALSE)[,1]
refname       <- gsub("/|\\||-|\\{|\\}", "_", gsub(">", "", gsub("A/", "", refname) ) )

refname_split <- strsplit( refname, "_", fixed = TRUE ) 
refname_split <- lapply(refname_split, 
                        function(x)
                          {
                          x_i  <- gsub(pattern = "[0-9]+", "", x)
                          x_ii <- gsub(pattern = "[a-zA-Z]+", "", x)
                          y = c(x_i, x_ii)
                          
                          y = y[!y == ""]
                          
                        }) 

tip_label_min <- gsub("_[0-9]{4}\\.[0-9]{3,4}","", h5_GsGDfile$tip.label)

treeid_split  <- strsplit( gsub("'", "", tip_label_min), "_", fixed = TRUE)
treeid_split  <- lapply(treeid_split, 
                        function(x)
                        {
                          x_i  <- gsub(pattern = "[0-9]+", "", x)
                          x_ii <- gsub(pattern = "[a-zA-Z]+", "", x)
                          y = c(x_i, x_ii)
                          
                          y = y[!y == ""]
                          
                        }) 

dirttymath <- sapply(refname_split, 
                     function(y)
                       {
                       return(
                         which.max(
                           sapply(treeid_split, 
                                  function(x)
                                  {
                                    sum( tolower(y) %in% tolower(x) )
                                  })
                                  )
                             )
                         })

# matchlist <- data.frame( ref = refname, 
#                          all = gsub(pattern = "'", "", h5_GsGDfile$tip.label[ dirttymath ]))



## map back to tree ----------------

h5_GsGDfile_clade <- h5_GsGDfile
h5_GsGDfile_clade$tip.label[-dirttymath] = " "
h5_GsGDfile_clade$tip.label[ dirttymath] = refname

ggtree(h5_GsGDfile, size = 0.2) 
ggtree(h5_GsGDfile_clade, size = 0.2) + geom_tiplab(size = 0.8, col = "RED")

# viewClade(g, node = identify(g)) + geom_tiplab(size = 0.5)

