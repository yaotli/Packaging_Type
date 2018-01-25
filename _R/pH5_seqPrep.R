# .fasta from NCBI: replace blank with underline

library(seqinr)
library(dplyr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/Nov2017_hana/")
source("~/Packaging_Type/_R/Function.R")


# gisaid
fas_ha_g <- "./source/pH5_G_2136_20171124.fasta" 
fas_na_g <- "./source/pNA_G_2135_20171124.fasta"
csv_g    <- "./Pre_process/G_2135_20171124.csv"


# ncbi
fas_ha_n <- "./Pre_process/pH5_N_5753_20171124.fasta" 
fas_na_n <- "./Pre_process/pNA_N_5721_20171124.fasta"


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
# result n = 7626

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
             file.out  = paste0( "./R_process/", "pH5_", length(s_infolist.ha[[4]]), ".fasta") )


write.fasta( sequences = s_infolist.na[[6]], 
             names     = paste0(s_infolist.na[[1]], "_",
                                s_infolist.na[[5]], "_|",
                                s_infolist.na[[3]], "|_",
                                s_infolist.na[[2]], "_",
                                s_infolist.na[[4]]
                                ),
             file.out  = paste0( "./R_process/", "pNA_", length(s_infolist.na[[4]]), ".fasta") )


## sequence curation and matching ----------------

# n = 7326
c_infolist.ha <- seqSelect( minlth = 1500, maxamb = 1, s_infolist.ha, rmdup = FALSE )
# n = 7277
c_infolist.na <- seqSelect( minlth = 1200, maxamb = 1, s_infolist.na, rmdup = FALSE )

# shared virus 
# n = 7138
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

#a                                       n        g        y
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

# export 

write.fasta( sequences = c_infolist.ha[[6]], 
             names     = paste0( c_infolist.ha[[1]], "_",
                                 c_infolist.ha[[5]], "_|",
                                 c_infolist.ha[[3]], "|_",
                                 c_infolist.ha[[2]], "_",
                                 c_infolist.ha[[4]]
             ),
             file.out  = paste0( "./R_process/pH5_", length( c_infolist.ha[[1]] ) ,".fasta")  )

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
             file.out  = paste0( "./R_process/pH5_", length( paired_infolist_ha[[1]] ) ,".fasta")  )


write.fasta( sequences = paired_infolist_na[[6]], 
             names     = paste0( paired_infolist_na[[1]], "_",
                                 paired_infolist_na[[5]], "_|",
                                 paired_infolist_na[[3]], "|_",
                                 paired_infolist_na[[2]], "_",
                                 paired_infolist_na[[4]]
                                 ),
             file.out  = paste0( "./R_process/pNA_", length( paired_infolist_na[[1]] ) ,".fasta") )

# only N1

subfastaSeq( subtype = "H5N1", filedir = "./R_process/pNA_7138.fasta" )


write.csv( file = "pH5NA.csv", x =
data.frame( ac.ha = paired_infolist_ha[[1]], 
            name  = paired_infolist_ha[[5]],
            geo   = paired_infolist_ha[[3]],
            sero  = paired_infolist_ha[[2]],
            year  = paired_infolist_ha[[4]],
            ac.na = paired_infolist_na[[1]], 
            stringsAsFactors = FALSE), row.names = FALSE)
  

### alignment, trim and tree --------------------------------

# mafft 
# system("mafft --reorder in.fasta > out.fasta")

# trim with BioEdit

# trim with trimtool for NA 
# trimtool( propblank = 0.9, filedir = "./align_trim/align_pN1_4468.fasta" )

# Fasttree
# system( "~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <in.fasta> in.tre" )

## GsGD ----------------

ls.gsgd <- tagExtra( "./tree/pH5_7326_e1130.tre" )$id[ !is.na(tagExtra( "./tree/pH5_7326_e1130.tre")$tag) ]
leafEx( "./tree/trim_pH5_7326_3.fasta", ls.gsgd )

# rmd 
rmDup( fasfile = "./tree/trim_pH5_7326_3.fasta", rmdup = TRUE )
ls.rmd.gsgd <- tagExtra( "./tree/rmd_pH5_7326/rmd_pH5_5855_e.tre")$id[ !is.na(tagExtra( "./tree/rmd_pH5_7326/rmd_pH5_5855_e.tre")$tag) ]
leafEx( "./tree/rmd_pH5_7326/rmd_pH5_5855.fasta", ls.rmd.gsgd)

### clade annotation and extraction --------------------------------

ann_6407 <- "./tree/pH5_7326_gsgd/pH5_7326_gsgd_e1130.tre"
pH5_seq  <- "./tree/trim_pH5_7326_3.fasta"
pN1_seq  <- "./align_trim/trim_pN1_4468_2.fasta"
p.table  <- "./R_process/pH5NA.csv"

ann_6407.tag     <- tagExtra( ann_6407 )
ann_6407.tag$tag <- gsub("ff8000", "c234", ann_6407.tag$tag)
ann_6407.tag$tag <- gsub("00ff00", "c232", ann_6407.tag$tag)

pTable <- read.csv( p.table, stringsAsFactors = FALSE)

## 234 ----------------

# rmd mL
ac_all_c234 <- str_match( ann_6407.tag$id[ which( ann_6407.tag$tag == "c234" ) ], "[A-Z0-9]+" )[,1]
ac_h5_c234  <- pTable$ac.ha[ na.omit( match( ac_all_c234, pTable$ac.ha ) ) ]
ac_n1_c234  <- pTable$ac.na[ intersect( na.omit( match( ac_all_c234, pTable$ac.ha ) ), 
                                        which( pTable$sero == "H5N1") ) ]

subfastaSeq( AC = TRUE, filedir = pH5_seq, AC_list = ac_h5_c234 )
subfastaSeq( AC = TRUE, filedir = pN1_seq, AC_list = ac_n1_c234 )

rmDup( fasfile = "./c234/raw/pH5_c234_2429.fasta", rmdup = TRUE )
rmDup( fasfile = "./c234/raw/pN1_c234_607.fasta", rmdup = TRUE )

rmdup_plus( "./c234/rmd/pH5_c234_1699.fasta" )
rmdup_plus( "./c234/rmd/pN1_c234_472.fasta" )

# cn
rmDup( fasfile = "./c234/raw/pH5_c234_2429.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )
rmDup( fasfile = "./c234/raw/pN1_c234_607.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )

reI.c234.ha <- tagExtra( "./filtering/cn/c234/pH5_c234_242.fasta.tre" )$id[ !is.na( tagExtra( "./filtering/cn/c234/pH5_c234_242.fasta.tre" )$tag ) ]
reI.c234.na <- tagExtra( "./filtering/cn/c234/pN1_c234_242.fasta.tre" )$id[ !is.na( tagExtra( "./filtering/cn/c234/pN1_c234_242.fasta.tre" )$tag ) ]

# reassort (most process below)

leafEx( "./c234/raw/pH5_c234_2429.fasta", setdiff( fastaEx( "./filtering/reassort/pH5_c234_rmR.fasta" )$id, reI.c234.ha ) )
leafEx( "./c234/raw/pN1_c234_607.fasta", setdiff( fastaEx( "./filtering/reassort/pN1_c234_rmR.fasta" )$id, reI.c234.na ) )

# pre_tree

rmDup( "./Post_filter/c234/pH5_c234_204.fasta", rmdup = TRUE )
rmDup( "./Post_filter/c234/pN1_c234_204.fasta", rmdup = TRUE )

rmdup_plus( "./Post_filter/c234/rmd/pH5_c234_179.fasta" )
rmdup_plus( "./Post_filter/c234/rmd/pN1_c234_167.fasta" )

# select (keep unique, control period before 2012)

leafEx( "./Post_filter/c234/pH5_c234_204.fasta", tagExtra("./Post_filter/c234/adv_rmd/ml/pH5_c234_160.phy_phyml.tre")$id[ is.na( tagExtra("./Post_filter/c234/adv_rmd/ml/pH5_c234_160.phy_phyml.tre")$tag ) ] )
leafEx( "./Post_filter/c234/pN1_c234_204.fasta", tagExtra("./Post_filter/c234/adv_rmd/ml/pN1_c234_153.phy_phyml.tre")$id[ is.na( tagExtra("./Post_filter/c234/adv_rmd/ml/pN1_c234_153.phy_phyml.tre")$tag ) ] )

# check clockness and prep faslist ---

pH5.trelist.234      <- taxaInfo( "./Post_filter/c234/adv_rmd/select/pH5_c234_158_e1204.tre", useTree = TRUE, root2tip = TRUE)
pN5.r2t.234          <- pH5.trelist.234[[ 8 ]]
pH5.trelist.234[[9]] <- geoID( strings = pH5.trelist.234[[6]], host = TRUE )
pH5.trelist.234[[9]][ which( pH5.trelist.234[[9]] == "Unknown") ] = "ML"

pN1.trelist.234      <- taxaInfo( "./Post_filter/c234/adv_rmd/select/pN1_c234_152_e1204.tre", useTree = TRUE, root2tip = TRUE)
pN1.r2t.234          <- pN1.trelist.234[[ 8 ]]
pN1.trelist.234[[9]] <- geoID( strings = pN1.trelist.234[[6]], host = TRUE )
pN1.trelist.234[[9]][ which( pN1.trelist.234[[9]] == "Unknown") ] = "ML"

# prep for states (eco) - Baye_state.R ---

# ha
pH5.trelist.234[[10]] = pH5.trelist.234[[9]]
pH5.trelist.234[[10]] <- 
  read.csv("./c234/eco_234.csv", stringsAsFactors = FALSE)$states[ 
    match( pH5.trelist.234[[6]], read.csv("./c234/eco_234.csv", stringsAsFactors = FALSE)$name ) ]

pH5.trelist.234[[10]][ is.na(pH5.trelist.234[[10]]) ] <- "W_a"
pH5.trelist.234[[10]] <- ifelse( startsWith(pH5.trelist.234[[10]], "D"), "D", "W")

#na
pN1.trelist.234[[10]] = pN1.trelist.234[[9]]
pN1.trelist.234[[10]] <- 
  read.csv("./c234/eco_234.csv", stringsAsFactors = FALSE)$states[ 
    match( gsub( pattern = "^[A-Z0-9]+", "", pN1.trelist.234[[6]] ), 
           gsub( pattern = "^[A-Z0-9]+", "", read.csv("./c234/eco_234.csv", stringsAsFactors = FALSE)$name ) ) ]

pN1.trelist.234[[10]][ is.na(pN1.trelist.234[[10]]) ] <- c("D_a", "W_a", "D_e", "D_e", "D_m",
                                                          "W_a", "D_a", "D_a", "D_a", "D_a",
                                                          "D_a", "D_a", "D_a", "W_a", "D_a")
pN1.trelist.234[[10]] <- ifelse( startsWith(pN1.trelist.234[[10]], "D"), "D", "W")





# export (for BEAST)
states.234 <- data.frame( id = pH5.trelist.234[[6]], states = pH5.trelist.234[[10]], stringsAsFactors = FALSE)
write.table( x = states.234, file = "eco.234", sep = "\t", quote = FALSE, row.names = FALSE)

# export .fasta
ac.c234_D <- acSearch( pH5.trelist.234, keyword = "D", keyword.dir = 10 )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c234_158.fasta", AC_list =  ac.c234_D, no = "D" )
ac.c234_W <- acSearch( pH5.trelist.234, keyword = "W", keyword.dir = 10 )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c234_158.fasta", AC_list =  ac.c234_W, no = "W" )



# p2
for( i in 1:4 )
{
  tem.y <- seq( 2004, 2011, 1 ) 
  
  y = acSearch( faslist = pH5.trelist.234, range = c( tem.y[ (i*2-1) ], tem.y[ i*2 ] ) )
  subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c234_158.fasta", AC_list =  y, no = tem.y[ i*2 ] )
  
  tem.m <- match( y, str_match( read.table( "./BEAST/eco.234", stringsAsFactors = FALSE, header = TRUE)$id, "[A-Z0-9]+") )
  tem.s <- read.table( "./BEAST/eco.234", stringsAsFactors = FALSE, header = TRUE)$states[ tem.m ]
  
  tem.table <- data.frame( id     = read.table( "./BEAST/eco.234", stringsAsFactors = FALSE, header = TRUE)$id[ tem.m ],
                           states = tem.s, stringsAsFactors = FALSE )
    
  write.table( x = tem.table, file = paste0( "eco.234.", tem.y[ i*2 ]) , sep = "\t", quote = FALSE, row.names = FALSE)
}

# p3

for( i in 1:3 )
{
  tem.y <- c(2002, 2007, 2012)
  
  y = acSearch( faslist = pH5.trelist.234, range = c( tem.y[ i ], (tem.y[ i ]+4) ) )
  subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c234_158.fasta", AC_list =  y, no = tem.y[ i ] )
  
  tem.m <- match( y, str_match( read.table( "./BEAST/eco.234", stringsAsFactors = FALSE, header = TRUE)$id, "[A-Z0-9]+") )
  tem.s <- read.table( "./BEAST/eco.234", stringsAsFactors = FALSE, header = TRUE)$states[ tem.m ]
  
  tem.table <- data.frame( id     = read.table( "./BEAST/eco.234", stringsAsFactors = FALSE, header = TRUE)$id[ tem.m ],
                           states = tem.s, stringsAsFactors = FALSE )
  
  write.table( x = tem.table, file = paste0( "eco.234.", tem.y[ i ]) , sep = "\t", quote = FALSE, row.names = FALSE)
}


# prep for states (geo) - Baye_state.R ---

pH5.trelist.234[[11]] <- 
  pH5_pri_cn_234[[9]]$geo[ match( pH5.trelist.234[[6]], pH5_pri_cn_234[[6]] ) ]
pH5.trelist.234[[11]][ is.na( pH5.trelist.234[[11]] ) ] <- "cnS"

write.table( x    = data.frame( id     = pH5.trelist.234[[6]], 
                                states = pH5.trelist.234[[11]] ),
             file = "./geo.234", sep = "\t", quote = FALSE, row.names = FALSE )


# clade sampling for geo ---

pH5.cladesamplingls.234 <- taxaInfo( file    = "./cladeSampling/234_h5_1205_E.tre", 
                                     useTree = TRUE, root2tip = TRUE)

pH5.cladesamplingls.234[[9]] <- read.table("./BEAST/geo.234", stringsAsFactors = FALSE, header = TRUE)$states[
  match( pH5.cladesamplingls.234[[6]], read.table("./BEAST/geo.234", stringsAsFactors = FALSE, header = TRUE)$id )]
  
cladeSampling( trefile   = "./cladeSampling/234_h5_1205_E.tre", 
               fasfile   = "./cladeSampling/input/pH5_c234_158.fasta", 
               suppList  = TRUE, 
               listinput = pH5.cladesamplingls.234,
               grid      = 1.1, 
               list.x    = c(6, 4, 9), 
               saveFasta = TRUE, showTree = TRUE)

temSample( fasfile    = "./cladeSampling/input/pH5_c234_158_cs.fasta", 
           faslist    = pH5.cladesamplingls.234, 
           list.x     = c(6,4,9), 
           samplelist = list( cnC  = c(2005, 5, 2006, 7, 2007, 3), 
                              cnS  = c(2005, 3, 2006, 4, 2007, 4),
                              cnSW = c(2005, 4, 2006, 7 ) ) )

write.table( x = data.frame( id     = fastaEx("./cladeSampling/pH5_c234_62_ts.fasta")$id,
                             states = pH5.cladesamplingls.234[[9]][ 
                               match( fastaEx("./cladeSampling/pH5_c234_62_ts.fasta")$id, pH5.cladesamplingls.234[[6]] )] ),
             file = "cladeSampling/cs.geo.234", sep = "\t", quote = FALSE, row.names = FALSE )



## 232 ----------------

# rmd ml
ac_all_c232 <- str_match( ann_6407.tag$id[ which( ann_6407.tag$tag == "c232" ) ], "[A-Z0-9]+" )[,1]
ac_h5_c232  <- pTable$ac.ha[ na.omit( match( ac_all_c232, pTable$ac.ha ) ) ]
ac_n1_c232  <- pTable$ac.na[ intersect( na.omit( match( ac_all_c232, pTable$ac.ha ) ), 
                                        which( pTable$sero == "H5N1") ) ]

subfastaSeq( AC = TRUE, filedir = pH5_seq, AC_list = ac_h5_c232 )
subfastaSeq( AC = TRUE, filedir = pN1_seq, AC_list = ac_n1_c232 )

rmDup( fasfile = "./c232/raw/pH5_c232_1296.fasta", rmdup = TRUE )
rmDup( fasfile = "./c232/raw/pN1_c232_1283.fasta", rmdup = TRUE )

rmdup_plus( "./c232/rmd/pH5_c232_1007.fasta" )
rmdup_plus( "./c232/rmd/pN1_c232_939.fasta" )


# cn
rmDup( fasfile = "./c232/raw/pH5_c232_1296.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )
rmDup( fasfile = "./c232/raw/pN1_c232_1283.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )

reI.c232.ha <- tagExtra( "./filtering/cn/c232/pH5_c232_263.fasta.tre" )$id[ !is.na( tagExtra( "./filtering/cn/c232/pH5_c232_263.fasta.tre" )$tag ) ]
reI.c232.na <- tagExtra( "./filtering/cn/c232/pN1_c232_263.fasta.tre" )$id[ !is.na( tagExtra( "./filtering/cn/c232/pN1_c232_263.fasta.tre" )$tag ) ]

# reassort (most process below)

leafEx( "./c232/raw/pH5_c232_1296.fasta", setdiff( fastaEx( "./filtering/reassort/pH5_c232_rmR.fasta" )$id, reI.c232.ha ) )
leafEx( "./c232/raw/pN1_c232_1283.fasta", setdiff( fastaEx( "./filtering/reassort/pN1_c232_rmR.fasta" )$id, reI.c232.na ) )

# pre_tree

rmDup( "./Post_filter/c232/pH5_c232_235.fasta", rmdup = TRUE )
rmDup( "./Post_filter/c232/pN1_c232_235.fasta", rmdup = TRUE )

rmdup_plus( "./Post_filter/c232/rmd/pH5_c232_213.fasta" )
rmdup_plus( "./Post_filter/c232/rmd/pN1_c232_200.fasta" )

# select (keep unique)

leafEx( "./Post_filter/c232/pH5_c232_235.fasta", tagExtra( "./Post_filter/c232/adv_rmd/ml/pH5_c232_208.phy_phyml.tre")$id[ is.na( tagExtra("./Post_filter/c232/adv_rmd/ml/pH5_c232_208.phy_phyml.tre")$tag ) ] )


# check clockness and prep faslist ---

pH5.trelist.232      <- taxaInfo( "./Post_filter/c232/adv_rmd/select/pH5_c232_207_e1204.tre", useTree = TRUE, root2tip = TRUE)
pN5.r2t.232          <- pH5.trelist.232[[ 8 ]]
pH5.trelist.232[[9]] <- geoID( strings = pH5.trelist.232[[6]], host = TRUE )
pH5.trelist.232[[9]][ which( pH5.trelist.232[[9]] == "Unknown") ] <- 
  c( rep("nonML", 4), rep("ML", 14) )
  
pN1.trelist.232      <- taxaInfo( "./Post_filter/c232/adv_rmd/select/pN1_c232_193_e1204.tre", useTree = TRUE, root2tip = TRUE)
pN1.r2t.232          <- pN1.trelist.232[[ 8 ]]
pN1.trelist.232[[9]] <- geoID( strings = pN1.trelist.232[[6]], host = TRUE )
pN1.trelist.232[[9]][ which( pN1.trelist.232[[9]] == "Unknown") ] <- 
  c( rep("nonML", 4), rep("ML", 14) )

# prep for states (eco) - Baye_state.R ---

# ha
pH5.trelist.232[[10]] = pH5.trelist.232[[9]]
pH5.trelist.232[[10]] <- 
  read.csv("./c232/eco_232.csv", stringsAsFactors = FALSE)$states[ 
    match( pH5.trelist.232[[6]], read.csv("./c232/eco_232.csv", stringsAsFactors = FALSE)$name ) ]

pH5.trelist.232[[10]][ is.na(pH5.trelist.232[[10]]) ] <- "W_a"
pH5.trelist.232[[10]] <- ifelse( startsWith(pH5.trelist.232[[10]], "D"), "D", "W")


# na
pN1.trelist.232[[10]] = pN1.trelist.232[[9]]
pN1.trelist.232[[10]] <- 
  read.csv("./c232/eco_232.csv", stringsAsFactors = FALSE)$states[ 
    match( gsub( pattern = "^[A-Z0-9]+", "", pN1.trelist.232[[6]] ), 
           gsub( pattern = "^[A-Z0-9]+", "", read.csv("./c232/eco_232.csv", stringsAsFactors = FALSE)$name ) ) ]

pN1.trelist.232[[10]][ is.na(pN1.trelist.232[[10]]) ] <- c("D_a", "D_a", "D_a", "D_a", "D_a",
                                                           "D_a", "D_a", "W_a", "D_e", "W_m",
                                                           "D_a", "W_a", "D_a", "W_a")
pN1.trelist.232[[10]] <- ifelse( startsWith(pN1.trelist.232[[10]], "D"), "D", "W")



# export (for BEAST)
states.232 <- data.frame( id = pH5.trelist.232[[6]], states = pH5.trelist.232[[10]], stringsAsFactors = FALSE)
write.table( x = states.232, file = "eco.232", sep = "\t", quote = FALSE, row.names = FALSE)

# export .fasta
ac.c234_ha_12 <- acSearch( pH5.trelist.232, range = c(2002, 2011) )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  ac.c234_ha_12, no = "ha" )
ac.c234_na_12 <- acSearch( pN1.trelist.232, range = c(2002, 2011) )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pN1_c232_193.fasta", AC_list =  ac.c234_na_12, no = "na" )

ac.c232_D <- acSearch( pH5.trelist.232, keyword = "D", keyword.dir = 10 )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  ac.c232_D, no = "D" )
ac.c232_D2 <- acSearch( pH5.trelist.232, keyword = "D", keyword.dir = 10, range = c(2002, 2011) )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  ac.c232_D2, no = "D" )

ac.c232_W <- acSearch( pH5.trelist.232, keyword = "W", keyword.dir = 10 )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  ac.c232_W, no = "W" )
ac.c232_W2 <- acSearch( pH5.trelist.232, keyword = "W", keyword.dir = 10, range = c(2002, 2011) )
subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  ac.c232_W2, no = "W" )



# p2
for( i in 1:6 )
{
  tem.y <- seq( 2004, 2015, 1 ) 
  
  y = acSearch( faslist = pH5.trelist.232, range = c( tem.y[ (i*2-1) ], tem.y[ i*2 ] ) )
  subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  y, no = tem.y[ i*2 ] )
  
  tem.m <- match( y, str_match( read.table( "./BEAST/eco.232", stringsAsFactors = FALSE, header = TRUE)$id, "[A-Z0-9]+") )
  tem.s <- read.table( "./BEAST/eco.232", stringsAsFactors = FALSE, header = TRUE)$states[ tem.m ]
  
  tem.table <- data.frame( id     = read.table( "./BEAST/eco.232", stringsAsFactors = FALSE, header = TRUE)$id[ tem.m ],
                           states = tem.s, stringsAsFactors = FALSE )
  
  write.table( x = tem.table, file = paste0( "eco.232.", tem.y[ i*2 ]) , sep = "\t", quote = FALSE, row.names = FALSE)
}

# p4
for( i in 1:3 )
{
  tem.y <- c(2002, 2007, 2012)
  
  y = acSearch( faslist = pH5.trelist.232, range = c( tem.y[ i ], (tem.y[ i ]+4) ) )
  subfastaSeq( AC = TRUE, filedir = "./BEAST/pH5_c232_207.fasta", AC_list =  y, no = tem.y[ i ] )
  
  tem.m <- match( y, str_match( read.table( "./BEAST/eco.232", stringsAsFactors = FALSE, header = TRUE)$id, "[A-Z0-9]+") )
  tem.s <- read.table( "./BEAST/eco.232", stringsAsFactors = FALSE, header = TRUE)$states[ tem.m ]
  
  tem.table <- data.frame( id     = read.table( "./BEAST/eco.232", stringsAsFactors = FALSE, header = TRUE)$id[ tem.m ],
                           states = tem.s, stringsAsFactors = FALSE )
  
  write.table( x = tem.table, file = paste0( "eco.232.", tem.y[ i ]) , sep = "\t", quote = FALSE, row.names = FALSE)
}


# prep for states (geo) - Baye_state.R ---

pH5.trelist.232[[11]] <- 
  pH5_pri_cn_232[[9]]$geo[ match( pH5.trelist.232[[6]], pH5_pri_cn_232[[6]] ) ]

pH5.trelist.232[[11]][ is.na( pH5.trelist.232[[11]] ) ] <- c( "cnN", rep(NA, 3) )

write.table( x    = data.frame( id     = pH5.trelist.232[[6]][which( floor(pH5.trelist.232[[4]]) < 2012 )], 
                                states = pH5.trelist.232[[11]][which( floor(pH5.trelist.232[[4]]) < 2012 )] ),
             file = "./geo.232", sep = "\t", quote = FALSE, row.names = FALSE )

# clade sampling for geo ---

pH5.cladesamplingls.232 <- taxaInfo( file    = "cladeSampling/232_h5_1205_E.tre", 
                                     useTree = TRUE, root2tip = TRUE)

pH5.cladesamplingls.232[[9]] <- pH5.trelist.232[[11]][
  match( pH5.cladesamplingls.232[[6]], pH5.trelist.232[[6]] )]

cladeSampling( trefile   = "./cladeSampling/232_h5_1205_E.tre", 
               fasfile   = "./cladeSampling/input/pH5_c232_207.fasta", 
               suppList  = TRUE, 
               listinput = pH5.cladesamplingls.232,
               grid      = 1.1, 
               list.x    = c(6, 4, 9), 
               saveFasta = TRUE, showTree = TRUE)

temSample( fasfile    = "./cladeSampling/input/pH5_c232_207_cs.fasta", 
           faslist    = pH5.cladesamplingls.232, 
           list.x     = c(6, 4, 9), 
           samplelist = list( cnE  = c(2009, 4, 2014, 5, 2015, 6), 
                              cnSW = c(2015, 6),
                              cnS  = c(2007, 3) ) )



write.table( x = data.frame( id     = fastaEx("./cladeSampling/pH5_c232_96_ts.fasta")$id,
                             states = pH5.cladesamplingls.232[[9]][ 
                               match( fastaEx("./cladeSampling/pH5_c232_96_ts.fasta")$id, pH5.cladesamplingls.232[[6]] )] ),
             file = "cladeSampling/cs.geo.232", sep = "\t", quote = FALSE, row.names = FALSE )








### lable seq. length on the tree --------------------------------

annlth <- function( seqfile, trefile )
{
  seq.s <- getSequence( read.fasta( seqfile ) )
  seq.n <- attributes( read.fasta( seqfile ) )$names
  tre   <- read.nexus( trefile )

  seq.s.l <- sapply( seq.s , 
                     function(x)
                     {
                       y = grep( pattern = "a|t|c|g", x = x, value = TRUE )
                       
                       return( length( y ) ) 
                     } )

  tre$tip.label <- paste0( tre$tip.label, "___", seq.s.l[ match( gsub( "'", "", tre$tip.label ), seq.n ) ] )
  
  write.tree( tre, file = gsub( ".tre", ".lth.tre", trefile ) )
} 



## prepare cleantre for examination ----------------

tre.pH5_7326_c23x  <- read.nexus( "./tree/pH5_7326_c23x.tre" )
tre.pN1_4468       <- read.nexus( "./tree/pN1_4468_e1203.tre" )

cleantre_pH5_7326_c23x           = tre.pH5_7326_c23x
cleantre_pH5_7326_c23x$tip.label = sub( pattern = "^'[A-Z0-9]+_", replacement = "'", x = cleantre_pH5_7326_c23x$tip.label )
write.tree( cleantre_pH5_7326_c23x, "cleantre_pH5_7326_c23x.tre" )

cleantre_pN1_4468           = tre.pN1_4468
cleantre_pN1_4468$tip.label = sub( "^'[A-Z0-9]+_", "'", cleantre_pN1_4468$tip.label )
write.tree( cleantre_pN1_4468, "cleantre_pN1_4468.tre" )


## examine reassortment by tangleplot ----------------

cleantxt_pH5       <- fortify( cleantre_pH5_7326_c23x )
cleantxt_pH5       <- cleantxt_pH5[ which( cleantxt_pH5$isTip ), ]
cleantxt_pH5$label <- gsub( "'", "", cleantxt_pH5$label )
cleantxt_pH5       <- data.frame( cleantxt_pH5,  gene = "H5", stringsAsFactors = FALSE )

cleantxt_pH5$gene[ which( cleantxt_pH5$y >= 1319 ) ] = "h5_234"
cleantxt_pH5$gene[ which( cleantxt_pH5$y < 1319 & cleantxt_pH5$y >= 5 ) ] = "h5_232"

cleantxt_pN1       <- fortify( cleantre_pN1_4468 )
cleantxt_pN1       <- cleantxt_pN1[ which( cleantxt_pN1$isTip ), ]
cleantxt_pN1$label <- gsub( "'", "", cleantxt_pN1$label )
cleantxt_pN1       <- data.frame( cleantxt_pN1,  gene = "N1", stringsAsFactors = FALSE )


cleantxt <- rbind(cleantxt_pH5, cleantxt_pN1 )
cleantxt <- data.frame( cleantxt, 
                        geo  = str_match( cleantxt$label, "\\|([A-Za-z_]+)\\|" )[, 2],
                        stringsAsFactors = FALSE )

# all
cleantxt %>%
filter( gene == "h5_232" | gene == "N1" ) %>% 
#filter( geo == "China" | geo == "Hong_Kong" ) %>% 
ggplot( aes(x = gene, y = y, group = label) ) + 
  geom_point( size = 0.01 ) +
  geom_line( alpha = 0.5 ) + 
  theme_bw() + ylab("") + xlab("") + 
  geom_hline( yintercept = 3200, color = "blue" ) +
  geom_hline( yintercept = 2750, color = "blue" ) +
  geom_hline( yintercept = 2100, color = "blue" )

#234
cleantxt_pHA_234 <- cleantxt_pH5[ which( cleantxt_pH5$gene == "h5_234" ), ]
tem.match        <- match( cleantxt_pHA_234$label, cleantxt_pN1$label )
cleantxt_pNA_234 <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
cleantxt_pHA_234 <- cleantxt_pHA_234[ which( cleantxt_pHA_234$label %in% cleantxt_pNA_234$label ), ]
tem.list         <- 
intersect(  cleantxt_pHA_234[ which( cleantxt_pHA_234$y < 2000 ), ]$label, 
            cleantxt_pNA_234[ which( cleantxt_pNA_234$y < 3500 & 
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
     cleantxt_pNA_232[ which( cleantxt_pNA_232$y < 2750 & 
                              cleantxt_pNA_232$y > 2100 ), ]$label )
ls.232 <- 
  grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )


# extract the seq 

ha.id <- fastaEx( "./tree/trim_pH5_7326_3.fasta" )$id
na.id <- fastaEx( "./tree/trim_pN1_4468_2.fasta" )$id


sampled.ls <- c( "ls.234", "ls.232" )
for( k in 1: 2 )
{
  sampled       <- get( sampled.ls[k] )
  tem.strain.ha <- c()
  tem.strain.na <- c()
  
  for( j in 1: length(sampled) )
  {
    tem.strain.ha <- c( tem.strain.ha, grep( sampled[j], ha.id, value = TRUE, fixed = TRUE ) ) 
    tem.strain.na <- c( tem.strain.na, grep( sampled[j], na.id, value = TRUE, fixed = TRUE ) ) 
  }
  
  subfastaSeq( AC = TRUE, filedir = "./tree/trim_pH5_7326_3.fasta", 
               AC_list = str_match( tem.strain.ha, "^[A-Z0-9]+")[, 1], no = sampled.ls[k] )
  subfastaSeq( AC = TRUE, filedir = "./tree/trim_pN1_4468_2.fasta", 
               AC_list = str_match( tem.strain.na, "^[A-Z0-9]+")[, 1], no = sampled.ls[k] )
}






### hyphy input --------------------------------

# FORM: clade(+host)_geo_time-gene

## 1205 ----------------

t.hyphy1205 <- c( "c232a_CN_0306_H5", "c232a_CN_0711_H5", "c232a_CN_1216_H5",  
                  "c232a_CN_0306_N1", "c232a_CN_0711_N1", "c232a_CN_1216_N1",
                  "c234a_CN_0306_H5", "c234a_CN_0711_H5", 
                  "c234a_CN_0306_N1", "c234a_CN_0711_N1" )

ac.hyphy1205 <- list( t.hyphy1205 )

ac.hyphy1205[[2]]  <- acSearch( faslist = pH5.trelist.232, keyword = "nonML", keyword.dir = 9, range = c(2003, 2006) )
ac.hyphy1205[[3]]  <- acSearch( faslist = pH5.trelist.232, keyword = "nonML", keyword.dir = 9, range = c(2007, 2011) )
ac.hyphy1205[[4]]  <- acSearch( faslist = pH5.trelist.232, keyword = "nonML", keyword.dir = 9, range = c(2012, 2016) )
ac.hyphy1205[[5]]  <- acSearch( faslist = pN1.trelist.232, keyword = "nonML", keyword.dir = 9, range = c(2003, 2006) )
ac.hyphy1205[[6]]  <- acSearch( faslist = pN1.trelist.232, keyword = "nonML", keyword.dir = 9, range = c(2007, 2011) )
ac.hyphy1205[[7]]  <- acSearch( faslist = pN1.trelist.232, keyword = "nonML", keyword.dir = 9, range = c(2012, 2016) )

ac.hyphy1205[[8]]  <- acSearch( faslist = pH5.trelist.234, keyword = "nonML", keyword.dir = 9, range = c(2003, 2006) )
ac.hyphy1205[[9]]  <- acSearch( faslist = pH5.trelist.234, keyword = "nonML", keyword.dir = 9, range = c(2007, 2011) )
ac.hyphy1205[[10]] <- acSearch( faslist = pN1.trelist.234, keyword = "nonML", keyword.dir = 9, range = c(2003, 2006) )
ac.hyphy1205[[11]] <- acSearch( faslist = pN1.trelist.234, keyword = "nonML", keyword.dir = 9, range = c(2007, 2011) )


tem.dir = c( rep( "./c232/raw/pH5_c232_1296.fasta", 3 ), rep( "./c232/raw/pN1_c232_1283.fasta", 3 ),
             rep( "./c234/raw/pH5_c234_2429.fasta", 2 ), rep( "./c234/raw/pN1_c234_607.fasta", 2 ) )

for(j in 2:11 )
{
  subfastaSeq( AC = TRUE, AC_list = ac.hyphy1205[[j]], filedir = tem.dir[j-1], no = ac.hyphy1205[[1]][j-1] )
}

pat_sample.ls <- paste0( "hyphy/1205/", list.files( "hyphy/1205/" ))

for( k in 1: length(pat_sample.ls) )
{
  if( grepl( "H5", pat_sample.ls[k] ) )
  {
    ntpartition( position = c(1: 1017), filedir = pat_sample.ls[k], no = "-HA1.fasta")
    ntpartition( position = c(1018: 1683), filedir = pat_sample.ls[k], no = "-HA2.fasta") 
    
  }else
  {
    ntpartition( position = c(211: 1347), filedir = pat_sample.ls[k], no = "-NAh.fasta")
    ntpartition( position = c(1: 210), filedir = pat_sample.ls[k], no = "-NAs.fasta")   
    }
}




## 1107 ----------------

# t.hyphy1107 = c( "c232a_CN_0406_H5", "c232a_CN_0711_H5", "c232a_CN_0406_N1", "c232a_CN_0711_N1", "c234a_CN_0406_H5", "c234a_CN_0711_H5", "c234a_CN_0406_N1", "c234a_CN_0711_N1")

ac.hyphy1107 <- list( t.hyphy1107 )

ac.hyphy1107[[2]] <- faslist.h5.232[[1]][ intersect( which( faslist.h5.232[[8]] == "nonML" ), which( floor( faslist.h5.232[[4]] ) == 2004 | floor( faslist.h5.232[[4]] ) == 2005 | floor( faslist.h5.232[[4]] ) == 2006 )  ) ] 
ac.hyphy1107[[3]] <- faslist.h5.232[[1]][ intersect( which( faslist.h5.232[[8]] == "nonML" ), which( floor( faslist.h5.232[[4]] ) == 2007 | floor( faslist.h5.232[[4]] ) == 2008 | floor( faslist.h5.232[[4]] ) == 2009 | floor( faslist.h5.232[[4]] ) == 2010 | floor( faslist.h5.232[[4]] ) == 2011 )  ) ] 
ac.hyphy1107[[4]] <- faslist.n1.232[[1]][ intersect( which( faslist.n1.232[[8]] == "nonML" ), which( floor( faslist.n1.232[[4]] ) == 2004 | floor( faslist.n1.232[[4]] ) == 2005 | floor( faslist.n1.232[[4]] ) == 2006 )  ) ] 
ac.hyphy1107[[5]] <- faslist.n1.232[[1]][ intersect( which( faslist.n1.232[[8]] == "nonML" ), which( floor( faslist.n1.232[[4]] ) == 2007 | floor( faslist.n1.232[[4]] ) == 2008 | floor( faslist.n1.232[[4]] ) == 2009 | floor( faslist.n1.232[[4]] ) == 2010 | floor( faslist.n1.232[[4]] ) == 2011 )  ) ] 

ac.hyphy1107[[6]] <- faslist.h5.234[[1]][ intersect( which( faslist.h5.234[[8]] == "nonML" ), which( floor( faslist.h5.234[[4]] ) == 2004 | floor( faslist.h5.234[[4]] ) == 2005 | floor( faslist.h5.234[[4]] ) == 2006 )  ) ] 
ac.hyphy1107[[7]] <- faslist.h5.234[[1]][ intersect( which( faslist.h5.234[[8]] == "nonML" ), which( floor( faslist.h5.234[[4]] ) == 2007 | floor( faslist.h5.234[[4]] ) == 2008 | floor( faslist.h5.234[[4]] ) == 2009 | floor( faslist.h5.234[[4]] ) == 2010 | floor( faslist.h5.234[[4]] ) == 2011 )  ) ] 
ac.hyphy1107[[8]] <- faslist.n1.234[[1]][ intersect( which( faslist.n1.234[[8]] == "nonML" ), which( floor( faslist.n1.234[[4]] ) == 2004 | floor( faslist.n1.234[[4]] ) == 2005 | floor( faslist.n1.234[[4]] ) == 2006 )  ) ] 
ac.hyphy1107[[9]] <- faslist.n1.234[[1]][ intersect( which( faslist.n1.234[[8]] == "nonML" ), which( floor( faslist.n1.234[[4]] ) == 2007 | floor( faslist.n1.234[[4]] ) == 2008 | floor( faslist.n1.234[[4]] ) == 2009 | floor( faslist.n1.234[[4]] ) == 2010 | floor( faslist.n1.234[[4]] ) == 2011 )  ) ] 


subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[2]], filedir = "./Tree/GsGD_pH5_6084.fasta", no = ac.hyphy1107[[1]][1] )
subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[3]], filedir = "./Tree/GsGD_pH5_6084.fasta", no = ac.hyphy1107[[1]][2] )
subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[6]], filedir = "./Tree/GsGD_pH5_6084.fasta", no = ac.hyphy1107[[1]][5] )
subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[7]], filedir = "./Tree/GsGD_pH5_6084.fasta", no = ac.hyphy1107[[1]][6] )

subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[4]], filedir = "./aligntrim_pH5NA/trim_pN1_4428_lth1347.fasta", no = ac.hyphy1107[[1]][3] )
subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[5]], filedir = "./aligntrim_pH5NA/trim_pN1_4428_lth1347.fasta", no = ac.hyphy1107[[1]][4] )
subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[8]], filedir = "./aligntrim_pH5NA/trim_pN1_4428_lth1347.fasta", no = ac.hyphy1107[[1]][7] )
subfastaSeq( AC = TRUE, AC_list = ac.hyphy1107[[9]], filedir = "./aligntrim_pH5NA/trim_pN1_4428_lth1347.fasta", no = ac.hyphy1107[[1]][8] )


pat_sample1 <- paste0( "~/Desktop/data_souce/Clade_pH5NA/hyphy_1107/HA/", 
                       list.files( "~/Desktop/data_souce/Clade_pH5NA/hyphy_1107/HA/" ))
for(k in 1: length( pat_sample1 ) )
{
  ntpartition( position = c(1: 1017), filedir = pat_sample1[k], no = "-HA1.fasta")
  ntpartition( position = c(1018: 1683), filedir = pat_sample1[k], no = "-HA2.fasta")
}

pat_sample2 <- paste0( "~/Desktop/data_souce/Clade_pH5NA/hyphy_1107/NA/", 
                       list.files( "~/Desktop/data_souce/Clade_pH5NA/hyphy_1107/NA/" ))
for(k in 1: length( pat_sample2 ) )
{
  ntpartition( position = c(211: 1347), filedir = pat_sample2[k], no = "-NAh.fasta")
  ntpartition( position = c(1: 210), filedir = pat_sample2[k], no = "-NAs.fasta")
}


## 0105 ----------------

# default: CN

ac.hyphy0117 <- list()
for( k in 1:24)
{
  fls.t <- c( rep( "pH5.trelist.232", 6), rep( "pN1.trelist.232", 6), rep( "pH5.trelist.234", 6), rep( "pN1.trelist.234",6 ) )
  eco.t <- rep( c( rep("Da", 3),  rep("Wa", 3) ), 4)
  r.t   <- rep( c( 2002, 2007, 2012), 8)
  
  ac.hyphy0105[[k]] <- acSearch( faslist = get(fls.t[k]), keyword = c("nonML", eco.t[k] ), keyword.dir = c(9, 10), range = c(r.t[k], (r.t[k]+4) ) )
}

t.hyphy0105 <- c( "c232a_D_0306_H5", "c232a_D_0711_H5", "c232a_D_1216_H5",
                  "c232a_W_0306_H5", "c232a_W_0711_H5", "c232a_W_1216_H5",  
                  "c232a_D_0306_N1", "c232a_D_0711_N1", "c232a_D_1216_N1",
                  "c232a_W_0306_N1", "c232a_W_0711_N1", "c232a_W_1216_N1",
                  "c234a_D_0306_H5", "c234a_D_0711_H5", "c234a_D_1216_H5", 
                  "c234a_W_0306_H5", "c234a_W_0711_H5", "c234a_W_1216_H5", 
                  "c234a_D_0306_N1", "c234a_D_0711_N1", "c234a_D_1216_N1",
                  "c234a_W_0306_N1", "c234a_W_0711_N1", "c234a_W_1216_N1" )

dir.t <- c( rep( "./c232/raw/pH5_c232_1296.fasta", 6 ), rep( "./c232/raw/pN1_c232_1283.fasta", 6 ), 
            rep( "./c234/raw/pH5_c234_2429.fasta", 6 ), rep( "./c234/raw/pN1_c234_607.fasta", 6 ) )

for( j in 1:24 )
{
  subfastaSeq( AC = TRUE, AC_list = ac.hyphy0105[[j]], filedir = dir.t[j], no = t.hyphy0105[j] )
}


pat_sample.ls <- paste0( "hyphy/0105/", list.files( "hyphy/0105//" ))

for( k in 1: length(pat_sample.ls) )
{
  if( grepl( "H5", pat_sample.ls[k] ) )
  {
    ntpartition( position = c(1: 1017), filedir = pat_sample.ls[k], no = "-HA1.fasta")
    ntpartition( position = c(1018: 1683), filedir = pat_sample.ls[k], no = "-HA2.fasta") 
    
  }else
  {
    ntpartition( position = c(211: 1347), filedir = pat_sample.ls[k], no = "-NAh.fasta")
    ntpartition( position = c(1: 210), filedir = pat_sample.ls[k], no = "-NAs.fasta")   
  }
}

## 0117 ----------------

# default: CN

ac.hyphy0117 <- list()
for( k in 1:24)
{
  fls.t <- c( rep( "pH5.trelist.232", 6), rep( "pN1.trelist.232", 6), rep( "pH5.trelist.234", 6), rep( "pN1.trelist.234",6 ) )
  eco.t <- rep( c( rep("D_a", 3),  rep("W_a", 3) ), 4)
  r.t   <- rep( c( 2002, 2007, 2012), 8)
  
  ac.hyphy0117[[k]] <- acSearch( faslist = get(fls.t[k]), keyword = c("nonML", eco.t[k] ), keyword.dir = c(9, 10), range = c(r.t[k], (r.t[k]+4) ) )
}

t.hyphy0117 <- c( "c232a_Da_0306_H5", "c232a_Da_0711_H5", "c232a_Da_1216_H5",
                  "c232a_Wa_0306_H5", "c232a_Wa_0711_H5", "c232a_Wa_1216_H5",  
                  "c232a_Da_0306_N1", "c232a_Da_0711_N1", "c232a_Da_1216_N1",
                  "c232a_Wa_0306_N1", "c232a_Wa_0711_N1", "c232a_Wa_1216_N1",
                  "c234a_Da_0306_H5", "c234a_Da_0711_H5", "c234a_Da_1216_H5", 
                  "c234a_Wa_0306_H5", "c234a_Wa_0711_H5", "c234a_Wa_1216_H5", 
                  "c234a_Da_0306_N1", "c234a_Da_0711_N1", "c234a_Da_1216_N1",
                  "c234a_Wa_0306_N1", "c234a_Wa_0711_N1", "c234a_Wa_1216_N1" )

dir.t <- c( rep( "./c232/raw/pH5_c232_1296.fasta", 6 ), rep( "./c232/raw/pN1_c232_1283.fasta", 6 ), 
            rep( "./c234/raw/pH5_c234_2429.fasta", 6 ), rep( "./c234/raw/pN1_c234_607.fasta", 6 ) )

for( j in 1:24 )
{
  subfastaSeq( AC = TRUE, AC_list = ac.hyphy0117[[j]], filedir = dir.t[j], no = t.hyphy0117[j] )
}


pat_sample.ls <- paste0( "hyphy/0117/", list.files( "hyphy/0117/" ))

for( k in 1: length(pat_sample.ls) )
{
  if( grepl( "H5", pat_sample.ls[k] ) )
  {
    ntpartition( position = c(1: 1017), filedir = pat_sample.ls[k], no = "-HA1.fasta")
    ntpartition( position = c(1018: 1683), filedir = pat_sample.ls[k], no = "-HA2.fasta") 
    
  }else
  {
    ntpartition( position = c(211: 1347), filedir = pat_sample.ls[k], no = "-NAh.fasta")
    ntpartition( position = c(1: 210), filedir = pat_sample.ls[k], no = "-NAs.fasta")   
  }
}

