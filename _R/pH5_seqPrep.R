# .fasta from NCBI: replace blank with underline

library(seqinr)
library(dplyr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")


# gisaid
fas_ha_g <- "./raw/pH5_G_2047_20170922.fasta" 
fas_na_g <- "./raw/pNA_G_2046_20170922.fasta"
csv_g    <- "./raw/pH5NA_G_2046_20170922.csv"


# ncbi
fas_ha_n <- "./raw/pH5_N_5675_20170919.fasta" 
fas_na_n <- "./raw/pNA_N_5643_20170919.fasta"
csv_ha_n <- "./raw/pH5_N_5675_20170919.csv"
csv_na_n <- "./raw/pNA_N_5643_20170919.csv"


### ha-na pairing --------------------------------

# GISAID 

seq.ha.g <- keepLongSeq( fastaEx(fas_ha_g)$seq, fastaEx(fas_ha_g)$id )$seq
id.ha.g  <- keepLongSeq( fastaEx(fas_ha_g)$seq, fastaEx(fas_ha_g)$id )$id
# length( seq.ha.g ) == length( fastaEx( fas_na_g )$seq )

seq.na.g <- fastaEx(fas_na_g)$seq
id.na.g  <- fastaEx(fas_na_g)$id

infolist.ha.g <- idInfo( rawid = id.ha.g, datasource = "g", g.csv = csv_g)
infolist.na.g <- idInfo( rawid = id.na.g, datasource = "g", g.csv = csv_g)

# setdiff( unique( infolist.na.g[[5]] ), unique( infolist.ha.g[[5]] ))
# setdiff( unique( infolist.ha.g[[5]] ), unique( infolist.na.g[[5]] ))

# according to the note column 
infolist.ha.g[[4]][grep("178249", infolist.ha.g[[1]])] = "2014-12_(Day_unknown)"
infolist.na.g[[4]][grep("178249", infolist.na.g[[1]])] = "2014-12_(Day_unknown)"




# NCBI 

id.ha.n <- fastaEx( fas_ha_n )$id
id.na.n <- fastaEx( fas_na_n )$id

# setdiff( unique( idInfo(id.ha.n)[[5]] ), unique( idInfo(id.na.n)[[5]] ) )
# [1] "yellow_billed_teal_Chile_14_2014"

remain = seq(1, length(id.ha.n) )[- grep( "yellow_billed_teal_Chile_14_2014", idInfo(id.ha.n)[[5]] )]

id.ha.n  <- fastaEx( fas_ha_n )$id[ remain ]
seq.ha.n <- fastaEx( fas_ha_n )$seq[ remain ]
seq.na.n <- fastaEx( fas_na_n )$seq

infolist.ha.n <- idInfo( id.ha.n )
infolist.na.n <- idInfo( id.na.n )

# setdiff( unique( infolist.ha.n[[5]] ), unique( infolist.na.n[[5]] ) )
# setdiff( unique( infolist.na.n[[5]] ), unique( infolist.ha.n[[5]] ) )


## combine and remove duplicated strain ----------------

seq.ha <- c( seq.ha.n, seq.ha.g )
seq.na <- c( seq.na.n, seq.na.g )

infolist.ha <- list( c(infolist.ha.n[[1]], infolist.ha.g[[1]] ), 
                     c(infolist.ha.n[[2]], infolist.ha.g[[2]] ),
                     c(infolist.ha.n[[3]], infolist.ha.g[[3]] ),
                     c(infolist.ha.n[[4]], infolist.ha.g[[4]] ),
                     c(infolist.ha.n[[5]], infolist.ha.g[[5]] ), 
                     seq.ha )


infolist.na <- list( c(infolist.na.n[[1]], infolist.na.g[[1]] ), 
                     c(infolist.na.n[[2]], infolist.na.g[[2]] ),
                     c(infolist.na.n[[3]], infolist.na.g[[3]] ),
                     c(infolist.na.n[[4]], infolist.na.g[[4]] ),
                     c(infolist.na.n[[5]], infolist.na.g[[5]] ), 
                     seq.na )


# do selection 
# result n = 7479

s_infolist.ha <- strainSelect( infolist.ha )
s_infolist.na <- strainSelect( infolist.na )

# time

s_infolist.ha[[4]] <- seqDate( s_infolist.ha[[4]] )
s_infolist.na[[4]] <- seqDate( s_infolist.na[[4]] )


write.fasta( sequences = s_infolist.ha[[6]], 
             names     = paste0(s_infolist.ha[[1]], "_",
                                s_infolist.ha[[5]], "_|",
                                s_infolist.ha[[3]], "|_",
                                s_infolist.ha[[2]], "_",
                                s_infolist.ha[[4]]
                                ),
             file.out  = "./pH5_7479.fasta")


write.fasta( sequences = s_infolist.na[[6]], 
             names     = paste0(s_infolist.na[[1]], "_",
                                s_infolist.na[[5]], "_|",
                                s_infolist.na[[3]], "|_",
                                s_infolist.na[[2]], "_",
                                s_infolist.na[[4]]
                                ),
             file.out  = "./pNA_7479.fasta")


## sequence curation and matching ----------------

c_infolist.ha <- seqSelect( minlth = 1500, maxamb = 1, s_infolist.ha, rmdup = FALSE)
c_infolist.na <- seqSelect( minlth = 1200, maxamb = 1, s_infolist.na, rmdup = FALSE)

# shared virus 
# n = 6974
int_hana <- intersect( c_infolist.ha[[5]], c_infolist.na[[5]] )


# check conflicting info
conf = c()

ha.n = c_infolist.ha[[5]][ match(int_hana, c_infolist.ha[[5]]) ]
na.n = c_infolist.na[[5]][ match(int_hana, c_infolist.na[[5]]) ]

ha.g = c_infolist.ha[[3]][ match(int_hana, c_infolist.ha[[5]]) ]
na.g = c_infolist.na[[3]][ match(int_hana, c_infolist.na[[5]]) ]

ha.y = as.numeric(c_infolist.ha[[4]][ match(int_hana, c_infolist.ha[[5]]) ])
na.y = as.numeric(c_infolist.na[[4]][ match(int_hana, c_infolist.na[[5]]) ])

for(i in 1: length(int_hana) )
{
  if( ha.n[i] != na.n[i] ){ conf = c(conf, i) } 
  if( ha.g[i] != na.g[i] ){ conf = c(conf, i) } 
  if( abs(ha.y[i] - na.y[i]) > 0.5 ){ conf = c(conf, i) } 
}

data.frame( a = c(c_infolist.ha[[1]][ match(int_hana, c_infolist.ha[[5]]) ][conf],
                  c_infolist.na[[1]][ match(int_hana, c_infolist.na[[5]]) ][conf]),
            n = c(ha.n[conf], na.n[conf]), 
            g = c(ha.g[conf], na.g[conf]),
            y = c(ha.y[conf], na.y[conf]))  

#a                              n        g        y
#1 AB233319 bar_headed_goose_Mongolia_1_05 Mongolia 2005.534
#2 KU143269       duck_Wenzhou_YHQL22_2014    China 2014.038
#3 KP638516         wildbird_Anhui_82_2005  Vietnam 2005.496
#4 AB239304 bar_headed_goose_Mongolia_1_05    Japan 2005.534
#5 KU143367       duck_Wenzhou_YHQL22_2014    China 2014.953
#6 KP638524         wildbird_Anhui_82_2005    China 2005.496


# unify info
c_infolist.ha[[3]][ match(int_hana, c_infolist.ha[[5]]) ][conf][3] = "China"

int_idx.ha <- match( int_hana, c_infolist.ha[[5]] )
int_idx.na <- match( int_hana, c_infolist.na[[5]] )  

# ha  
paired_infolist_ha <- list()
for( l in 1: length(c_infolist.ha) )
{
  paired_infolist_ha[[ l ]] = c_infolist.ha[[l]][int_idx.ha]
}

# na
paired_infolist_na <- list()
for( l in 1: length(c_infolist.na) )
{
  paired_infolist_na[[ l ]] = c_infolist.na[[l]][int_idx.na]
}


# force NA geo + time infro as HA

paired_infolist_na[[3]] = paired_infolist_ha[[3]]
paired_infolist_na[[4]] = paired_infolist_ha[[4]]

write.fasta( sequences = paired_infolist_ha[[6]], 
             names     = paste0( paired_infolist_ha[[1]], "_",
                                 paired_infolist_ha[[5]], "_|",
                                 paired_infolist_ha[[3]], "|_",
                                 paired_infolist_ha[[2]], "_",
                                 paired_infolist_ha[[4]]
                                ),
             file.out  = "./pH5_6974.fasta")


write.fasta( sequences = paired_infolist_na[[6]], 
             names     = paste0( paired_infolist_na[[1]], "_",
                                 paired_infolist_na[[5]], "_|",
                                 paired_infolist_na[[3]], "|_",
                                 paired_infolist_na[[2]], "_",
                                 paired_infolist_na[[4]]
                                 ),
             file.out  = "./pNA_6974.fasta")

# only N1

subfastaSeq( subtype = "H5N1", filedir = "./pNA_6974.fasta" )


write.csv( file = "pairedH5NA_6974.csv", x =
data.frame( ac.ha = paired_infolist_ha[[1]], 
            name  = paired_infolist_ha[[5]],
            geo   = paired_infolist_ha[[3]],
            sero  = paired_infolist_ha[[2]],
            year  = paired_infolist_ha[[4]],
            ac.na = paired_infolist_na[[1]], 
            stringsAsFactors = FALSE), row.names = FALSE)
  

### alignment, trim and tree --------------------------------

# mafft 

#system("mafft --reorder pH5_6974.fasta > align_pH5_6974.fasta")
#system("mafft --reorder pN1_4428.fasta > align_pN1_4428.fasta")

# trim with BioEdit
# trim with trimtool for NA 

#trimtool( propblank = 0.9, filedir = "./aligntrim_pH5NA/align_pN1_4428.fasta")

# Fasttree

#system("~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./aligntrim_pH5NA/trim3_pH5_6974_lth1683.fasta> ./tree/pH5_6974_lth1683.tre")
#system("~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./aligntrim_pH5NA/trim_pN1_4428_lth1347.fasta> ./tree/pN1_4428_lth1347.tre")


## GsGD ----------------

# manually prepare .txt file containing taxa names of GsGD 

subtreseq( list_filedir = "./Tree/GsGDlist_pH5_6084.txt", 
           seq_filedir  = "./aligntrim_pH5NA/trim3_pH5_6974_lth1683.fasta")


#system("~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./Tree/GsGD_pH5_6084.fasta> ./tree/GsGD_pH5_6084.tre")


### clade annotation and extraction --------------------------------

nwk_GsGD_pH5_6084 <- "./Tree/GsGD_pH5_6084_nwk"
ann_GsGD_pH5_6084 <- "./Tree/GsGD_pH5_6084_clade.tre"
fas_GsGD_pH5_6084 <- "./Tree/GsGD_pH5_6084.fasta"

nwk_pN1_4428      <- "./Tree/pN1_4428_nwk"
fas_pN1_4428      <- "./aligntrim_pH5NA/trim_pN1_4428_lth1347.fasta"

pairedcsv         <- "./pairedH5NA_6974.csv"


# H5
tre_GsGD_pH5_6084 <- read.tree( nwk_GsGD_pH5_6084 )
trefile_GsGD_pH5  <- fortify( tre_GsGD_pH5_6084 )
root.node1        <- length( tre_GsGD_pH5_6084$tip.label ) + 1
Root.dis1         <- dist.nodes( tre_GsGD_pH5_6084 )[ root.node1, 1: (root.node1 - 1) ]
tre.id1           <- gsub("'", "", trefile_GsGD_pH5$label[ 1: root.node1- 1 ] )

tre.id1.a <- str_match( tre.id1, "[A-Z]{1,2}[0-9]{5,6}|EPI[0-9]+" )[,1]
tre.id1.g <- str_match( tre.id1, "\\|([A-Za-z_]+)\\|" )[,2]
tre.id1.s <- str_match( tre.id1, "_(H5N[0-9]{1,2})_" )[,2]
tre.id1.y <- as.numeric( str_match(tre.id1, "_([0-9]{4}.[0-9]{3})$" )[,2] )


# use labeled .tre
tem.clade         <- tagExtra( ann_GsGD_pH5_6084 )

tem.clade$tag[ is.na(tem.clade$tag) ] = "0"
tem.clade$tag     <- gsub("ff3366", "7", 
                          gsub("ff0033", "22", 
                               gsub("cc0033", "232", 
                                    gsub("990033", "234", tem.clade$tag) )))

tretable_GsGD_pH5 <- data.frame( tre.id1.a, tre.id1.g, tre.id1.s, tre.id1.y, 
                                 Distance = Root.dis1, 
                                 Clade    = tem.clade$tag[ match(tre.id1, tem.clade$id) ], 
                                 tre.id1, 
                                 stringsAsFactors = FALSE)

# extract clade 
for (i in 1: 4)
{
  clade <- c("234", "232", "22", "7")
  
  subfastaSeq( AC      = TRUE, 
               no      = clade[i],
               AC_list = tretable_GsGD_pH5$tre.id1.a[ which(tretable_GsGD_pH5$Clade == clade[i] ) ],
               filedir = fas_GsGD_pH5_6084)
}



# extract N1
tre_pN1_4428 <- read.tree( nwk_pN1_4428 )
paired_table <- read.csv( pairedcsv, stringsAsFactors = FALSE )

for (i in 1: 4)
{
  clade <- c("234", "232", "22", "7")
  
  n1_ac <- 
  paired_table[,6][ match( tretable_GsGD_pH5$tre.id1.a[ which( tretable_GsGD_pH5$Clade     == clade[i] & 
                                                               tretable_GsGD_pH5$tre.id1.s == "H5N1") ], 
                           paired_table[,1]) ]
  
  subfastaSeq( AC      = TRUE, 
               no      = clade[i],
               AC_list = n1_ac,
               filedir = fas_pN1_4428)
}

# system( "for f in $(ls ~/Desktop/data_souce/Clade_pH5NA/*.fasta); do ~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <$f> $f.tre ; done" )


### clade-wise --------------------------------

annlth <- function( seqfile, trefile )
{
  seq.s <- getSequence( read.fasta( seqfile ) )
  seq.n <- attributes( read.fasta( seqfile ) )$names
  tre   <- read.tree( trefile )

  seq.s.l <- sapply( seq.s , 
                     function(x)
                     {
                       y = gsub( "-|~", "", c2s(x) )
                       return( length( s2c(y) ) ) 
                     } )

  tre$tip.label <- paste0( tre$tip.label, "___", seq.s.l[ match( tre$tip.label, seq.n ) ] )
  
  write.tree( tre, file = gsub( ".tre", ".lth.tre", trefile ) )
} 

## 234 ----------------

pH5_234            <- "./Clade_pH5NA/234/raw/c234_pH5_2311.fasta"
pN1_234            <- "./Clade_pH5NA/234/raw/c234_pN1_610.fasta"
pH5_234.s1.tre     <- "./Clade_pH5NA/sampled/pH5_c234_s1.tre"
pN1_234.s1.tre     <- "./Clade_pH5NA/sampled/pN1_c234_s1.tre"
pH5_234.s2.tre     <- "./Clade_pH5NA/ml_phyml/c234_pH5_176.phy_phyml_tree.tre"
pN1_234.s2.tre     <- "./Clade_pH5NA/ml_phyml/c234_pN1_164.phy_phyml_tree.tre"
pH5_234.lthann.tre <- "./Clade_pH5NA/ml_phyml/c234_pH5_176.phy_phyml.lth.tree.lth.tre"
pN1_234.lthann.tre <- "./Clade_pH5NA/ml_phyml/c234_pN1_164.phy_phyml.lth.tree.lth.tre"

# s1 to elimicate re-introduction
AC_pH5_234_s1 <- str_match( tagExtra( pH5_234.s1.tre )[, 1][ which( is.na( tagExtra( pH5_234.s1.tre )[, 2] ) ) ], "^[A-Z00-9]+" )[, 1]
AC_pN1_234_s1 <- str_match( tagExtra( pN1_234.s1.tre )[, 1][ which( is.na( tagExtra( pN1_234.s1.tre )[, 2] ) ) ], "^[A-Z00-9]+" )[, 1]

subfastaSeq( AC = TRUE, AC_list = AC_pH5_234_s1, filedir = pH5_234 )
subfastaSeq( AC = TRUE, AC_list = AC_pN1_234_s1, filedir = pN1_234 )

# rmDup
rmDup( fasfile = "./Clade_pH5NA/234/c234_pH5_s2.fasta", year = c(1000, 2012), rmdup = TRUE )
rmDup( fasfile = "./Clade_pH5NA/234/c234_pN1_s2.fasta", year = c(1000, 2012), rmdup = TRUE )

# ann. tree with sequence length 
annlth( seqfile = pH5_234, trefile = pH5_234.s2.tre )
annlth( seqfile = pN1_234, trefile = pN1_234.s2.tre )        

# extract tag taxa
AC_pH5_234_rmLead <-         
  str_match( tagExtra( pH5_234.lthann.tre )[,1][ grep( "ff0000", tagExtra( pH5_234.lthann.tre )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]
AC_pN1_234_rmLead <-         
  str_match( tagExtra( pN1_234.lthann.tre )[,1][ grep( "ff0000", tagExtra( pN1_234.lthann.tre )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]

subfastaSeq( AC = TRUE, AC_list = AC_pH5_234_rmLead, filedir = pH5_234 )
subfastaSeq( AC = TRUE, AC_list = AC_pN1_234_rmLead, filedir = pN1_234 )

pH5.trelist.234 <- taxaInfo( file = "./Clade_pH5NA/ml_phyml_2/c234_pH5_156e.tre", useTree = TRUE, root2tip = TRUE)
pH5.r2t.234     <- pH5.trelist.234[[ 8 ]]
pN1.trelist.234 <- taxaInfo( file = "./Clade_pH5NA/ml_phyml_2/c234_pN1_149e.tre", useTree = TRUE, root2tip = TRUE)
pN1.r2t.234     <- pN1.trelist.234[[ 8 ]]


## 232 ----------------

pH5_232            <- "./Clade_pH5NA/232/raw/c232_pH5_1255.fasta"
pN1_232            <- "./Clade_pH5NA/232/raw/c232_pN1_1244.fasta"
pH5_232.s1.tre     <- "./Clade_pH5NA/sampled/pH5_c232_s1.tre"
pN1_232.s1.tre     <- "./Clade_pH5NA/sampled/pN1_c232_s1.tre"
pH5_232.s2.tre     <- "./Clade_pH5NA/ml_phyml/c232_pH5_126.phy_phyml_tree.tre"
pN1_232.s2.tre     <- "./Clade_pH5NA/ml_phyml/c232_pN1_121.phy_phyml_tree.tre"
pH5_232.lthann.tre <- "./Clade_pH5NA/ml_phyml/c232_pH5_126.phy_phyml.lth.tree.lth.tre"
pN1_232.lthann.tre <- "./Clade_pH5NA/ml_phyml/c232_pN1_121.phy_phyml.lth.tree.lth.tre"

# s1 for re-intro.
AC_pH5_232_s1 <- str_match( tagExtra( pH5_232.s1.tre )[, 1][ which( is.na( tagExtra( pH5_232.s1.tre )[, 2] ) ) ], "^[A-Z00-9]+" )[, 1]
AC_pN1_232_s1 <- str_match( tagExtra( pN1_232.s1.tre )[, 1][ which( is.na( tagExtra( pN1_232.s1.tre )[, 2] ) ) ], "^[A-Z00-9]+" )[, 1]

subfastaSeq( AC = TRUE, AC_list = AC_pH5_232_s1, filedir = pH5_232)
subfastaSeq( AC = TRUE, AC_list = AC_pN1_232_s1, filedir = pN1_232)

# rmDup
rmDup( fasfile = "./Clade_pH5NA/232/c232_pH5_s2.fasta", year = c(1000, 2012), rmdup = TRUE)
rmDup( fasfile = "./Clade_pH5NA/232/c232_pN1_s2.fasta", year = c(1000, 2012), rmdup = TRUE)

# ann. tree with sequence length 
annlth( seqfile = pH5_232, trefile = pH5_232.s2.tre )
annlth( seqfile = pN1_232, trefile = pN1_232.s2.tre )      

AC_pH5_232_rmLead <-         
  str_match( tagExtra( pH5_232.lthann.tre )[,1][ grep( "ff0000", tagExtra( pH5_232.lthann.tre )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]
AC_pN1_232_rmLead <-         
  str_match( tagExtra( pN1_232.lthann.tre )[,1][ grep( "ff0000", tagExtra( pN1_232.lthann.tre )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]

subfastaSeq( AC = TRUE, AC_list = AC_pH5_232_rmLead, filedir = pH5_232)
subfastaSeq( AC = TRUE, AC_list = AC_pN1_232_rmLead, filedir = pN1_232)

pH5.trelist.232 <- taxaInfo( file = "./Clade_pH5NA/ml_phyml_2/c232_pH5_120e.tre", useTree = TRUE, root2tip = TRUE)
pH5.r2t.232     <- pH5.trelist.232[[ 8 ]]
pN1.trelist.232 <- taxaInfo( file = "./Clade_pH5NA/ml_phyml_2/c232_pN1_115e.tre", useTree = TRUE, root2tip = TRUE)
pN1.r2t.232     <- pN1.trelist.232[[ 8 ]]


# 22 ----------------

pH5_22        <- "./Clade_pH5NA/22/raw/c22_pH5_1129.fasta"
pN1_22        <- "./Clade_pH5NA/22/raw/c22_pN1_1129.fasta"
pH5_22.s1.tre <- "./Clade_pH5NA/sampled/pH5_c22_s1.tre"
pN1_22.s1.tre <- "./Clade_pH5NA/sampled/pN1_c22_s1.tre"
pH5_22.s2.tre <- "./Clade_pH5NA/ml_phyml/c22_pH5_46.phy_phyml_tree.tre"
pN1_22.s2.tre <- "./Clade_pH5NA/ml_phyml/c22_pN1_39.phy_phyml_tree.tre"
pH5_22.lthann.tre <- "./Clade_pH5NA/ml_phyml/c22_pH5_46.phy_phyml.lth.tree.lth.tre"
pN1_22.lthann.tre <- "./Clade_pH5NA/ml_phyml/c22_pN1_39.phy_phyml.lth.tree.lth.tre"

# s1 for re-introduction
AC_pH5_22_s1 <- str_match( tagExtra( pH5_22.s1.tre )[, 1][ which( is.na( tagExtra( pH5_22.s1.tre )[, 2] ) ) ], "^[A-Z00-9]+" )[, 1]
AC_pN1_22_s1 <- str_match( tagExtra( pN1_22.s1.tre )[, 1][ which( is.na( tagExtra( pN1_22.s1.tre )[, 2] ) ) ], "^[A-Z00-9]+" )[, 1]

subfastaSeq( AC = TRUE, AC_list = AC_pH5_22_s1, filedir = pH5_22)
subfastaSeq( AC = TRUE, AC_list = AC_pN1_22_s1, filedir = pN1_22)

# rmDup
rmDup( fasfile = "./Clade_pH5NA/22/c22_pH5_s2.fasta", year = c(1000, 2012), rmdup = TRUE)
rmDup( fasfile = "./Clade_pH5NA/22/c22_pN1_s2.fasta", year = c(1000, 2012), rmdup = TRUE)

# ann tree with sequence length 
annlth( seqfile = pH5_22, trefile = pH5_22.s2.tre )
annlth( seqfile = pN1_22, trefile = pN1_22.s2.tre )        

AC_pH5_22_rmLead <-         
  str_match( tagExtra( pH5_22.lthann.tre )[,1][ grep( "ff0000", tagExtra( pH5_22.lthann.tre )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]
AC_pN1_22_rmLead <-         
  str_match( tagExtra( pN1_22.lthann.tre )[,1][ grep( "ff0000", tagExtra( pN1_22.lthann.tre )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]


subfastaSeq( AC = TRUE, AC_list = AC_pH5_22_rmLead, filedir = pH5_22)
subfastaSeq( AC = TRUE, AC_list = AC_pN1_22_rmLead, filedir = pN1_22)

pH5.trelist.22 <- taxaInfo( file = "./Clade_pH5NA/ml_phyml_2/c22_pH5_41e.tre", useTree = TRUE, root2tip = TRUE)
pH5.r2t.22     <- pH5.trelist.22[[ 8 ]]
pN1.trelist.22 <- taxaInfo( file = "./Clade_pH5NA/ml_phyml_2/c22_pN1_39e.tre", useTree = TRUE, root2tip = TRUE)
pN1.r2t.22     <- pN1.trelist.22[[ 8 ]]



# system( "for f in $(ls *phy); do ~/PhyML-3.1/./phyml -i $f -m GTR -v e -a 4; done" )

# 7

## Root-to-tip plot ----------------

ggplot() + 
  theme_bw() + 
  xlab("") + ylab("") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_text(size = 18)) + 
  scale_x_continuous( limit = c(2000, 2017), 
                      breaks = seq(2000, 2016, by = 2), 
                      labels = seq(2000, 2016, by = 2)) +
  
  geom_point( data = tretable_GsGD_pH5[which(tretable_GsGD_pH5$tre.id1.g == "China"|
                                             tretable_GsGD_pH5$tre.id1.g == "Hong_Kong"), ],
              aes( x = tre.id1.y, y = Distance), size = 3, alpha = 0.4, stroke = 0) +
  
  facet_grid( Clade ~.)
  
  
  #geom_point( data = tretable_GsGD_pH5,
  #            aes( x = tre.id1.y, y = Distance), 
  #            color = "gray", size = 3, alpha = 0.1, stroke = 0) + 
  #geom_point( data = tretable_GsGD_pH5[ which(tretable_GsGD_pH5$Clade == "234"), ],
  #            aes( x = tre.id1.y, y = Distance), 
  #            color = "red", size = 3, alpha = 0.4, stroke = 0) + 
  #geom_point( data = tretable_GsGD_pH5[ which(tretable_GsGD_pH5$Clade == "232"), ],
  #            aes( x = tre.id1.y, y = Distance), 
  #           color = "green", size = 3, alpha = 0.4, stroke = 0) + 
  #geom_point( data = tretable_GsGD_pH5[ which(tretable_GsGD_pH5$Clade == "22"), ],
  #            aes( x = tre.id1.y, y = Distance), 
  #            color = "darkorange", size = 3, alpha = 0.4, stroke = 0) + 
  #geom_point( data = tretable_GsGD_pH5[ which(tretable_GsGD_pH5$Clade == "7"), ],
  #            aes( x = tre.id1.y, y = Distance), 
  #            color = "black", size = 3, alpha = 0.4, stroke = 0) 
  


## prepare cleantre for examination ----------------

cleantre_GsGD_pH5_6084           = tre_GsGD_pH5_6084
cleantre_GsGD_pH5_6084$tip.label = 
  sub("[A-Z]{1,2}[0-9]{5,6}_|EPI[0-9]+_", "", 
      cleantre_GsGD_pH5_6084$tip.label)

write.tree(cleantre_GsGD_pH5_6084, "cleantre_GsGD_pH5_6084.tre")

cleantre_pN1_4428           = tre_pN1_4428
cleantre_pN1_4428$tip.label = 
  sub("[A-Z]{1,2}[0-9]{5,6}_|EPI[0-9]+_", "", 
      cleantre_pN1_4428$tip.label)

write.tree(cleantre_pN1_4428, "cleantre_pN1_4428.tre")


## examine reassortment ----------------

cleantxt_pH5       <- fortify( cleantre_GsGD_pH5_6084 )
cleantxt_pH5       <- cleantxt_pH5[ which( cleantxt_pH5$isTip ), ]
cleantxt_pH5$label <- gsub( "'", "", cleantxt_pH5$label )
cleantxt_pH5       <- data.frame( cleantxt_pH5,  gene = "H5", stringsAsFactors = FALSE )

cleantxt_pH5$gene[ which( cleantxt_pH5$y >= 3774 ) ] = "h5_234"
cleantxt_pH5$gene[ which( cleantxt_pH5$y < 3774 & cleantxt_pH5$y >= 2519 ) ] = "h5_232"
cleantxt_pH5$gene[ which( cleantxt_pH5$y < 2519 & cleantxt_pH5$y >= 1390 ) ] = "h5_22"


cleantxt_pN1       <- fortify( cleantre_pN1_4428 )
cleantxt_pN1       <- cleantxt_pN1[ which( cleantxt_pN1$isTip ), ]
cleantxt_pN1$label <- gsub( "'", "", cleantxt_pN1$label )
cleantxt_pN1       <- data.frame( cleantxt_pN1,  gene = "N1", stringsAsFactors = FALSE )


cleantxt <- rbind(cleantxt_pH5, cleantxt_pN1 )
cleantxt <- data.frame( cleantxt, 
                        geo  = str_match( cleantxt$label, "\\|([A-Za-z_]+)\\|" )[, 2],
                        stringsAsFactors = FALSE )


# all
cleantxt %>%
filter( gene == "h5_22" | gene == "N1" ) %>% 
#filter( geo == "China" | geo == "Hong_Kong" ) %>% 
ggplot( aes(x = gene, y = y, group = label) ) + 
  geom_point( size = 0.01 ) +
  geom_line() + 
  theme_bw() + ylab("") + xlab("") + 
  geom_hline( yintercept = 1505, color = "red" ) + 
  geom_hline( yintercept = 340, color = "red" ) +
  geom_hline( yintercept = 2150, color = "red" )

#234
cleantxt_pHA_234 <- cleantxt_pH5[ which( cleantxt_pH5$gene == "h5_234" ), ]
tem.match        <- match( cleantxt_pHA_234$label, cleantxt_pN1$label )
cleantxt_pNA_234 <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
cleantxt_pHA_234 <- cleantxt_pHA_234[ which( cleantxt_pHA_234$label %in% cleantxt_pNA_234$label ), ]
tem.list         <- 
intersect(  cleantxt_pHA_234[ which( cleantxt_pHA_234$y < 4500 ), ]$label, 
            cleantxt_pNA_234[ which( cleantxt_pNA_234$y < 3400 & 
                                     cleantxt_pNA_234$y > 2500 ), ]$label )
ls.234 <- 
grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )


#232
cleantxt_pHA_232 <- cleantxt_pH5[ which( cleantxt_pH5$gene == "h5_232" ), ]
tem.match        <- match( cleantxt_pHA_232$label, cleantxt_pN1$label )
cleantxt_pNA_232 <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
cleantxt_pHA_232 <- cleantxt_pHA_232[ which( cleantxt_pHA_232$label %in% cleantxt_pNA_232$label ), ]
tem.list         <- 
  c( cleantxt_pNA_232[ which( cleantxt_pNA_232$y > 3200 ), ]$label, 
     cleantxt_pNA_232[ which( cleantxt_pNA_232$y < 2800 & 
                              cleantxt_pNA_232$y > 2150 ), ]$label )
ls.232 <- 
  grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )


#22
cleantxt_pHA_22 <- cleantxt_pH5[ which( cleantxt_pH5$gene == "h5_22" ), ]
tem.match       <- match( cleantxt_pHA_22$label, cleantxt_pN1$label )
cleantxt_pNA_22 <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
cleantxt_pHA_22 <- cleantxt_pHA_22[ which( cleantxt_pHA_22$label %in% cleantxt_pNA_22$label ), ]
tem.list        <- cleantxt_pNA_22[ which( cleantxt_pNA_22$y < 1505 & cleantxt_pNA_22$y > 340 ), ]$label 
ls.22 <- 
  grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )


# extract the seq 

ha.id <- fastaEx( fas_GsGD_pH5_6084 )$id
na.id <- fastaEx( fas_pN1_4428 )$id


sampled.ls <- c( "ls.234", "ls.232", "ls.22" )
for( k in 1: 3 )
{
  
  sampled       <- get( sampled.ls[k] )
  tem.strain.ha <- c()
  tem.strain.na <- c()
  
  for( j in 1: length(sampled) )
  {
    tem.strain.ha <- c( tem.strain.ha, grep( sampled[j], ha.id, value = TRUE, fixed = TRUE ) ) 
    tem.strain.na <- c( tem.strain.na, grep( sampled[j], na.id, value = TRUE, fixed = TRUE ) ) 
  }
  
  subfastaSeq( AC = TRUE, filedir = fas_GsGD_pH5_6084, 
               AC_list = str_match( tem.strain.ha, "^[A-Z00-9]+")[, 1], no = sampled.ls[k] )
  subfastaSeq( AC = TRUE, filedir = fas_pN1_4428, 
               AC_list = str_match( tem.strain.na, "^[A-Z00-9]+")[, 1], no = sampled.ls[k] )
  
}

