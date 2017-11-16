library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(tidyverse)
library(ggridges)


setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

# gisaid 

fas_allh5_g <- "./raw/H5_G_2092_20170905.fasta" 
csv_allh5_g <- "./raw/H5_G_2092_20170905.csv"

# ncbi 

fas_allh5_n <- "./raw/H5_N_7363_20170905.fasta"
csv_allh5_n <- "./raw/H5_N_7363_20170905.csv"


### read in & parsing ------------------------------

## GISAID ----------------
# n = 2092
# 2 sequence with the same accession number &
# only 2091 sequence info in csv file

allh5_g_seq      <- keepLongSeq( fastaEx( fas_allh5_g )$seq, fastaEx( fas_allh5_g )$id )$seq
allh5_g_id       <- keepLongSeq( fastaEx( fas_allh5_g )$seq, fastaEx( fas_allh5_g )$id )$id
infolist.allh5.g <- idInfo( rawid = allh5_g_id, datasource = "g", g.csv = csv_allh5_g)

# based on the Note column in csv 
infolist.allh5.g[[4]][ grep("178249", infolist.allh5.g[[1]] ) ] = "2014-12_(Day_unknown)"



## NCBI ----------------
# n = 7363

allh5_n_seq      <- fastaEx( fas_allh5_n )$seq
allh5_n_id       <- fastaEx( fas_allh5_n )$id
infolist.allh5.n <- idInfo( rawid = allh5_n_id )


## deal with duplicated strain name ----------------

allh5_seq      <- c( allh5_n_seq, allh5_g_seq )
infolist.allh5 <- list( c( infolist.allh5.n[[1]], infolist.allh5.g[[1]] ), 
                        c( infolist.allh5.n[[2]], infolist.allh5.g[[2]] ),
                        c( infolist.allh5.n[[3]], infolist.allh5.g[[3]] ),
                        c( infolist.allh5.n[[4]], infolist.allh5.g[[4]] ),
                        c( infolist.allh5.n[[5]], infolist.allh5.g[[5]] ),
                        allh5_seq )

# do the selection

s_infolist.allh5 <- strainSelect( infolist.allh5 )

# transform year format

s_infolist.allh5[[4]] <- seqDate( s_infolist.allh5[[4]] )


write.fasta( sequences = s_infolist.allh5[[6]], 
                 names = paste0(s_infolist.allh5[[1]], "_", 
                               s_infolist.allh5[[5]], "_|",
                               s_infolist.allh5[[3]], "|_",
                               s_infolist.allh5[[2]], "_",
                               s_infolist.allh5[[4]]
                               ), 
              file.out = "./allh5_9136.fasta")


## curation sequences ----------------

# remove sequence with > 1 % ambiguous nucleotides
# remove sequence < 1500 nt 
# remove duplicated sequences 

c_infolist.allh5 <- seqSelect( minlth = 1500, maxamb = 1, s_infolist.allh5, rmdup = TRUE)

write.fasta( sequences = c_infolist.allh5[[6]], 
             names = paste0(c_infolist.allh5[[1]], "_", 
                           c_infolist.allh5[[5]], "_|",
                           c_infolist.allh5[[3]], "|_",
                           c_infolist.allh5[[2]], "_",
                           c_infolist.allh5[[4]]
             ), 
             file.out = "./c_allh5_7368.fasta")


### tree and clade annotation --------------------------------

# mafft
# system("mafft --reorder c_allh5_7368.fasta > align_allh5_7368.fasta")

# manuall trim in BioEdit & Seqotron (final lth. = 1683)
# Fasttree
# system("~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./aligntrim_allh5/trim_allh5_7368_lth1683.fasta> ./tree/allh5_7368_lth1683.tre")


# manual extract GsGD lineages
# eliminate 2 sequence with apparent wrong phylogenetic inferences
# KT936697_common_teal_Shanghai_PD1108_8_2013_China_H5N1_2013.852
# KT936701_common_teal_Shanghai_PD1108_4_2013_China_H5N1_2013.852


subtreseq( seq_filedir  = "./aligntrim_allh5/trim_allh5_7368_lth1683.fasta", 
           list_filedir = "./Tree/allh5_GsGDlike_6322.txt" )

# re-Fasttree

# system("~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./tree/allh5_GsGDlike_6322.fasta> ./tree/allh5_GsGDlike_6322.tre")



### GsGD global circulation --------------------------------

trefile_GsGDlike_clade <- "./Tree/allh5_GsGDlike_6322_e0903.tre"

# readin

trelist_GsGD <- taxaInfo( trefile_GsGDlike_clade, useTree = TRUE )

col.code  <- unique( trelist_GsGD[[7]] )
col.clade <- c( "2.5", "1", "5-6", "2.2", "p",
                "2.3.4", "2.3.2", "7", "2.1", "0",
                "3", "2.4", "4", "8-9", "2.3.1", 
                "2", "2.3", "2.3.3")

trelist_GsGD[[7]] <- col.clade[ match(trelist_GsGD[[7]], col.code) ]

trelist_GsGD[[8]] <- trelist_GsGD[[7]]

for(k in 1: length( trelist_GsGD[[8]] ))
{
  if( trelist_GsGD[[8]][k] == "2.3.3" ){ trelist_GsGD[[8]][k] = "2.3.4"}
  if( trelist_GsGD[[8]][k] == "2.3.1" ){ trelist_GsGD[[8]][k] = "2.3.2"}
  if( trelist_GsGD[[8]][k] == "2.5" ){ trelist_GsGD[[8]][k] = "2.2" }
  
}

### clade-wise  --------------------------------

# extract clade 
clade <- c("2.3.4", "2.3.2", "2.2", "7")
for(i in 1: 4)
{
 subfastaSeq( AC      = TRUE, 
              no      = clade[i], 
              AC_list = trelist_GsGD[[1]][ which( trelist_GsGD[[8]] == clade[i] ) ],
              filedir = "./Tree/allh5_GsGDlike_6322.fasta")
  
}

# system("for f in $(ls ~/Desktop/data_souce/Clade_allh5/*.fasta); do ~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <$f> $f.tre ; done")



### 234 --------------------------------

rmDup( fasfile = "./Clade_allh5/234/raw/c234_1930.fasta", 
       sero    = "H5N1", 
       year    = c(1000, 2012) )

faslist.234 <- taxaInfo( file = "./Clade_allh5/234/rmdup/c234_508.fasta" )

## primary data ----------------
# geo

faslist.234[[8]] <- geoID( strings = faslist.234[[6]] )
faslist.234[[8]][ which( faslist.234[[8]] == "Unknown" ) ] = 
  c( "cnE", "cnE", "cnE", "cnE", "cnN", 
     "cnN", "cnE", "cnE", "cnC", "cnE", 
     "cnE", "cnC", "Unknown", "cnS", "cnS")
  
faslist.234[[2]][232] = "China"

# table( floor(faslist.234[[4]]), faslist.234[[8]] )
# epi.
# df.geo <- data.frame( y = faslist.234[[4]], g = faslist.234[[8]] )
#
# ggplot(data = df.geo, aes( x = y ) ) + 
#   geom_density( stat = 'bin', binwidth = 0.5) + 
#   facet_wrap(~g, ncol = 1) + xlab("")


# host
faslist.234[[9]] = geoID( faslist.234[[6]], host = TRUE)
faslist.234[[9]][ which(faslist.234[[9]] == "Unknown") ] = "ML"


# modification
# table( faslist.234[[8]] )
# cnBH     cnC     cnE     cnN    cnNE    cnNW     cnS 
# 11      49       8       2       1       3      63 
# cnSW    cnYZ     SEA Unknown 
# 74      18     278       1 
# 
# remove : Unknown ( EF624256_China_2006_|China|_H5N1_2006.496 ) 
# merge  : cnE / cnBH / cnYZ
# merge  : cnN / cnNW / cnNE

faslist.234[[8]][ which( faslist.234[[8]] == "Unknown" ) ] = NA

faslist.234[[8]][ which( faslist.234[[8]] == "cnBH" ) ]    = "cnE"
faslist.234[[8]][ which( faslist.234[[8]] == "cnYZ" ) ]    = "cnE"
faslist.234[[8]][ which( faslist.234[[8]] == "cnNW" ) ]    = "cnN"
faslist.234[[8]][ which( faslist.234[[8]] == "cnNE" ) ]    = "cnN"

# cnC  cnE  cnN  cnS cnSW  SEA 
# 49   37    6   63   74  278 

# remove duplicated sequence (seq. lth)
annlth( seqfile = "./Clade_allh5/234/rmdup/c234_508.fasta", 
        trefile = "./Clade_allh5/234/rmdup/c234_508.phy_phyml_e.tre" )

AC_234_rmLeaf <-         
  str_match( tagExtra( "./Clade_allh5/234/rmdup/c234_508.phy_phyml_e.lth.tre" )[,1][ grep( "ff0000", tagExtra( "./Clade_allh5/234/rmdup/c234_508.phy_phyml_e.lth.tre" )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]

subfastaSeq( AC = TRUE, AC_list = AC_234_rmLeaf, filedir = "./Clade_allh5/234/raw/c234_1930.fasta" )

# cladesampling ( 471 )
# remove : cnN 
faslist.234[[8]][ which( faslist.234[[8]] == "cnN" ) ] = NA

cladeSampling( trefile   = "Clade_allh5/234/rmdup/c234_471.phy_phyml_e.tre", 
               fasfile   = "Clade_allh5/234/rmdup/c234_508.fasta", 
               suppList  = TRUE, 
               listinput = faslist.234,
               grid      = 1, 
               list.x    = c(6, 4, 8), 
               saveFasta = TRUE, 
               showTree  = TRUE )

#       cnC cnE cnS cnSW SEA
# 2004   1   0   0    0   0
# 2005  10   3   7   11   4
# 2006  12   2   8   25   2
# 2007   7   3   6    0  34
# 2008   1   4   3    3  21
# 2009   3   4   0    3  12
# 2010   0   9   0    0  10
# 2011   1   3   0    0   1

# temSample

temSample( fasfile    = "./Clade_allh5/234/clade_tem_sample/c234_213_cs.fasta", 
           faslist    = faslist.234, 
           list.x     = c(6,4,8), 
           samplelist = list( cnC  = c(2005, 5, 2006, 6), 
                              cnE  = c(2010, 5),
                              cnSW = c(2005, 6, 2006, 13),
                              SEA  = c(2007, 12, 2008, 10, 2009, 6, 2010, 5) ) )

# cnC  cnE  cnS cnSW  SEA 
# 24   24   24   25   40 

list_c234_ts    <- taxaInfo( file = "./Clade_allh5/234/clade_tem_sample/c234_137.ts.phy_phyml_e.tre", useTree = TRUE, root2tip = TRUE)
ts.trait.234    <- data.frame( id  = list_c234_ts[[6]], 
                               geo = faslist.234[[8]][ match( list_c234_ts[[6]], faslist.234[[6]] ) ],
                               stringsAsFactors = FALSE )
write.table( x = ts.trait.234, file = "ts.trait.234", sep = "\t", quote = FALSE, row.names = FALSE)




# simple curate for phylogeo.
trelist.234_508     <- taxaInfo( file = "./Clade_allh5/234/rmdup/c234_508_clade.tre", useTree = TRUE, root2tip = TRUE)
trelist.234_508.r2t <- trelist.234_508[[8]]

# View( trelist.234_508.r2t$outlier )
# to remove: 
# EF624256_China_2006_ ( unknown isolation place )
# KC709802_chicken_Shandong_K0701_2010_ (304, divergence)
# KF169906_chicken_Hong_Kong_8825_2_2008_ (305, divergence)

ac.505 <- trelist.234_508[[1]][-which( trelist.234_508[[1]] %in% c("EF624256", "KC709802", "KF169906") ) ]
subfastaSeq( AC = TRUE, AC_list = ac.505, filedir = "./Clade_allh5/234/rmdup/c234_508.fasta")

# prepare table for beast


b.trait.234_505  <- data.frame( id  = taxaInfo( file = "./Clade_allh5/234/rmdup/c234_505.fasta")[[6]], 
                                geo = faslist.234[[8]][ match(taxaInfo( file = "./Clade_allh5/234/rmdup/c234_505.fasta")[[6]], faslist.234[[6]]) ], 
                                stringsAsFactors = FALSE )

write.table( x = b.trait.234_505, file = "b.trait.234_505", sep = "\t", quote = FALSE, row.names = FALSE)

## secondary data ----------------
# downsampling

cladeSampling( trefile   = "Clade_allh5/234/rmdup/c234_508_clade.tre", 
               fasfile   = "Clade_allh5/234/rmdup/c234_508.fasta", 
               suppList  = TRUE, 
               listinput = faslist.234,
               grid      = 1, 
               list.x    = c(6, 4, 8), 
               saveFasta = TRUE)

# remove over-divergence and reassortant

trelist.234     <- taxaInfo( file = "./Clade_allh5/234/downsample_int/c234_208_e.tre", useTree = TRUE, root2tip = TRUE)
trelist.234.r2t <- trelist.234[[8]]

# View( trelist.234.r2t$outlier )
# to remove: 
# EF624256_China_2006_ ( unknown isolation place )
# KC709802_chicken_Shandong_K0701_2010_ (151, divergence)
# KP735806_goose_Shandong_k1201_2009_ (18, divergence, reassortant Nx)
# JX534565_wild_duck_Shandong_628_2011_ (19, divergence, reassortant Nx)
# KF169906_chicken_Hong_Kong_8825_2_2008_ (152, divergence)
# EPI80661_chicken_Vietnam_NCVD_008_2008_ (185, divergence)

# prepare table for BEAST

b.trait.234.id <- taxaInfo( file = "./Clade_allh5/234/c234_202_e.tre", useTree = TRUE)[[6]]
b.trait.234    <- data.frame( id  = b.trait.234.id, 
                              geo = faslist.234[[8]][ match( b.trait.234.id, faslist.234[[6]] ) ],
                              stringsAsFactors = FALSE )
write.table( x = b.trait.234, file = "b.trait.234", sep = "\t", quote = FALSE, row.names = FALSE)


### 232 --------------------------------

rmDup( fasfile = "./Clade_allh5/232/raw/c232_1187.fasta", 
       sero    = "H5N1", 
       year    = c(1000, 2012) )

faslist.232 <- taxaInfo( file = "./Clade_allh5/232/rmdup/c232_446.fasta" )


## primary data ----------------
# geo
faslist.232[[8]] <- geoID( strings = faslist.232[[6]] )
faslist.232[[8]][ which( faslist.232[[8]] == "Unknown" ) ] = 
  c("cnN", "cnNW", "SEA", "cnC", "cnC", "cnC", "cnS")

# epi.
# df.geo <- data.frame( y = faslist.232[[4]], g = faslist.232[[8]] )
# 
# ggplot(data = df.geo, aes( x = y ) ) + 
#  geom_density( stat = 'bin', binwidth = 0.5) + 
#  facet_wrap(~g, ncol = 1) + xlab("")


# host

faslist.232[[9]] = geoID( strings = faslist.232[[6]], host = TRUE)
faslist.232[[9]][ which(faslist.232[[9]] == "Unknown") ][9:12] = "nonML"
faslist.232[[9]][ which(faslist.232[[9]] == "Unknown") ]       = "ML"


# modification
# table( faslist.232[[8]] )
# 
# cnBH  cnC  cnN cnNE cnNW  cnS cnSW cnYZ    E   nA  SEA 
#    5   39    1    2   18   39   42    9    1  115  175 
#
# remove: E ( CY110854_common_buzzard_Bulgaria_38WB_2010_|Bulgaria|_H5N1_2010.200 )
# merge : cnNW / cnN / cnNE
# merge : cnBH / cnYZ

faslist.232[[8]][ which( faslist.232[[8]] == "E" ) ]    = NA
faslist.232[[8]][ which( faslist.232[[8]] == "cnNW" | faslist.232[[8]] == "cnN" | faslist.232[[8]] == "cnNE" ) ] = "cnN"
faslist.232[[8]][ which( faslist.232[[8]] == "cnYZ" | faslist.232[[8]] == "cnBH" ) ] = "cnE"

# cnC  cnE  cnN  cnS cnSW   nA  SEA 
#  39   14   21   39   42  115  175 


# remove duplicated sequence (seq. lth)
annlth( seqfile = "./Clade_allh5/232/rmdup/c232_446.fasta", 
        trefile = "./Clade_allh5/232/rmdup/c232_446.phy_phyml_e.tre" )

AC_232_rmLeaf <-         
  str_match( tagExtra( "./Clade_allh5/232/rmdup/c232_446.phy_phyml_e.lth.tre" )[,1][ grep( "ff0000", tagExtra( "./Clade_allh5/232/rmdup/c232_446.phy_phyml_e.lth.tre" )[,2], invert = TRUE) ], "^[A-Za-z0-9]+" )[,1]

subfastaSeq( AC = TRUE, AC_list = AC_232_rmLeaf, filedir = "./Clade_allh5/232/raw/c232_1187.fasta" )

# cladesampling ( 160 )

cladeSampling( trefile   = "Clade_allh5/232/rmdup/c232_428.phy_phyml_e.tre", 
               fasfile   = "Clade_allh5/232/rmdup/c232_446.fasta", 
               suppList  = TRUE, 
               listinput = faslist.232,
               grid      = 1, 
               list.x    = c(6, 4, 8), 
               saveFasta = TRUE, showTree = TRUE)


#       cnC cnE cnN cnS cnSW nA SEA
# 2003   1   0   0   0    0  0   1
# 2004   2   0   0   4    3  0   0
# 2005   2   1   0   1    9  0   4
# 2006   0   0   0   2    0  0   1
# 2007   3   0   2   6    0  0   0
# 2008   4   0   0   6    0  1   2
# 2009   7   5   2   4    2  9   2
# 2010   2   1   2   3    0  7  17
# 2011   5   3   2   2    1  1  27

temSample( fasfile    = "./Clade_allh5/232/clade_tem_sample/c232_159.cs.fasta", 
           faslist    = faslist.232, 
           list.x     = c(6,4,8), 
           samplelist = list( cnC  = c(2008, 2, 2009, 3, 2011, 2), 
                              cnS  = c(2004, 2, 2007, 3, 2008, 3),
                              SEA  = c(2005, 2, 2010, 8, 2011, 10 ) ) )

# cnC  cnE  cnN  cnS cnSW   nA  SEA 
#  17   10    8   20   15   18   26 

list_c232_ts    <- taxaInfo( file = "./Clade_allh5/232/clade_tem_sample/c232_114.ts.phy_phyml_e.tre", useTree = TRUE, root2tip = TRUE)
ts.trait.232    <- data.frame( id  = list_c232_ts[[6]], 
                               geo = faslist.232[[8]][ match( list_c232_ts[[6]], faslist.232[[6]] ) ],
                               stringsAsFactors = FALSE )
write.table( x = ts.trait.232, file = "ts.trait.232", sep = "\t", quote = FALSE, row.names = FALSE)








# simple curate for phylogeo.
trelist.232_446     <- taxaInfo( file = "./Clade_allh5/232/rmdup/c232_446_clade.tre", useTree = TRUE, root2tip = TRUE)
trelist.232_446.r2t <- trelist.232_446[[8]]

# View( trelist.232_446.r2t$outlier )
# to remove: 
# CY110854_common_buzzard_Bulgaria_38WB_2010_ (Europe)
# JN807989_eurasian_eagle_owl_Korea_Q133_2011_ (431, divergence)
# JN807991_eurasian_sparrowhawk_Korea_Q94_2011_ (432, divergence)

ac.443 <- trelist.232_446[[1]][-which( trelist.232_446[[1]] %in% c("CY110854", "JN807989", "JN807991") ) ]
subfastaSeq( AC = TRUE, AC_list = ac.443, filedir = "./Clade_allh5/232/rmdup/c232_446.fasta")

# prepare table for beast


b.trait.232_443  <- data.frame( id  = taxaInfo( file = "./Clade_allh5/232/rmdup/c232_443.fasta")[[6]], 
                                geo = faslist.232[[8]][ match(taxaInfo( file = "./Clade_allh5/232/rmdup/c232_443.fasta")[[6]], faslist.232[[6]]) ], 
                                stringsAsFactors = FALSE )
write.table( x = b.trait.232_443, file = "b.trait.232_443", sep = "\t", quote = FALSE, row.names = FALSE)




## secondary data ----------------
# downsampling

cladeSampling( trefile   = "Clade_allh5/232/rmdup/c232_446_clade.tre", 
               fasfile   = "Clade_allh5/232/rmdup/c232_446.fasta", 
               suppList  = TRUE, 
               listinput = faslist.232,
               grid      = 1, 
               list.x    = c(6, 4, 8), 
               saveFasta = TRUE)

# remove over-divergence 

trelist.232     <- taxaInfo( file = "./Clade_allh5/232/downsample_int/c232_152_e.tre", useTree = TRUE, root2tip = TRUE)
trelist.232.r2t <- trelist.232[[8]]

# View( trelist.232.r2t$outlier )
# to remove:
# CY110854_common_buzzard_Bulgaria_38WB_2010_ (Europe)
# EPI169419_Chicken_Shandong_J06_2009_ (126, divergence)
# EPI169421_Duck_Shandong_Y02_2009_ (128, divergence)

b.trait.232.id <- taxaInfo( file = "./Clade_allh5/232/c232_149_e.tre", useTree = TRUE)[[6]]
b.trait.232    <- data.frame( id  = b.trait.232.id, 
                              geo = faslist.232[[8]][ match( b.trait.232.id, faslist.232[[6]] ) ],
                              stringsAsFactors = FALSE )

write.table( x = b.trait.232, file = "b.trait.232", sep = "\t", quote = FALSE, row.names = FALSE)



## joy plot  ----------------

# prepare data 

table.f <- data.frame( s = trelist_GsGD[[3]], g = trelist_GsGD[[2]], 
                       y = trelist_GsGD[[4]], cl = trelist_GsGD[[7]], 
                       n = trelist_GsGD[[6]], 
                       stringsAsFactors = FALSE )

table.f1 = 
  table.f               %>% 
  filter( cl != "p" )   %>%
  filter( s != "H5N0" ) %>%
  select( s, g, y, cl)  %>%
  mutate( China = ifelse( g == "China" | g == "Hong_Kong", "China", "Other countries") )  %>%
  mutate( Serotype = ifelse( s == "H5N1", "N1", "Nx") )                   
  
  Clade = table.f1$cl

  for(k in 1: length(Clade))
  {
    if( Clade[k] == "2.4" ){ Clade[k] = "2"}
    if( Clade[k] == "2.5" ){ Clade[k] = "2"}
    if( Clade[k] == "2.3" ){ Clade[k] = "2"}
    
    if( Clade[k] == "2.3.3" ){ Clade[k] = "2.3.4" }
    if( Clade[k] == "2.3.1" ){ Clade[k] = "2.3.2" }
    
    if( Clade[k] == "3" ){ Clade[k] = "Old" }
    if( Clade[k] == "4" ){ Clade[k] = "Old" }
    if( Clade[k] == "5-6" ){ Clade[k] = "Old" }
    if( Clade[k] == "8-9" ){ Clade[k] = "Old"}
  }
    
table.f1 <- data.frame( Serotype = table.f1$Serotype, 
                        Geo      = table.f1$China, 
                        Clade    = Clade, 
                        Year     = table.f1$y, stringsAsFactors = FALSE)


table.f1$Year  <- str_match( table.f1$Year, pattern = "[0-9]{4}" )[,1]

table.f1       <- data.frame( table( table.f1 ), stringsAsFactors = FALSE)
table.f1$Year  <- as.numeric( as.character( table.f1$Year ) )

table.f1$Clade <- factor( table.f1$Clade, 
                          levels = c("0", "Old", "1", "2", "2.1", "2.2", "7","2.3.2", "2.3.4") )

table.f1$Clade <- gsub( pattern = "Old", replacement = "3", x = table.f1$Clade)
table.f1$Clade <- factor(table.f1$Clade, 
                         levels = c("0", "3", "1", "2", "2.1", "2.2", "7","2.3.2", "2.3.4") )


# ggjoy

ggplot( table.f1, aes(x = Year-1996, y = Clade, height = Freq, fill = Serotype) ) + 
geom_density_ridges( stat = "identity", scale = 3, alpha = 0.8 ) + 
facet_wrap(~ Geo) +
scale_x_continuous( breaks = seq(0, 21, by = 2), labels = seq(1996, 2017, by = 2), 
                    limit = c(0, 21)) + 
  
scale_y_discrete( expand = c(0.001, 0)) +
scale_fill_cyclical( values = c("#2ca02c", "#d62728"), guide = "legend") + 
theme_minimal() + 
theme( panel.grid.major.y =  element_blank(), 
       panel.grid.minor.x =  element_blank(), 
       strip.background = element_rect(fill="white", color = "white"), 
       text = element_text(family = "Helvetica"),
       strip.text.x     = element_text(size = 20), 
       axis.text.y      = element_text(size = 20), 
       legend.position = c(0,1) ) +
labs(x = "", y = "")

# ggsave("a.pdf", width = 8, height = 6)


## Root-to-tip plot ----------------

trefile_GsGDlike_nwk <- "./Tree/allh5_GsGDlike_6322_nwk"
  
GsGDlike_tree     <- read.tree( trefile_GsGDlike_nwk )
GsGDlike_treedata <- fortify( GsGDlike_tree )
rootnode          <- length( GsGDlike_tree$tip.label ) + 1

GsGDlike_dis      <- dist.nodes( GsGDlike_tree )[rootnode, 1: (rootnode - 1)]
GsGDlike_tree_id  <- gsub("'", "", GsGDlike_treedata$label[1: rootnode - 1 ])

GsGDlike_table    <- cbind( GsGDlike_table, 
                            Divergence = GsGDlike_dis[ match(GsGDlike_table$GsGDlike_name, GsGDlike_tree_id) ] )


# plot 

ggplot() + 
  theme_bw() + 
  xlab("") + ylab("Root-to-tip divergence (sub./site)") + 
  
  scale_x_continuous( limits = c(2000, 2016),
                      breaks = seq(2000, 2016, by = 2),
                      labels = seq(2000, 2016, by = 2) ) +
  # all
  geom_point( data = GsGDlike_table[ which(GsGDlike_table$GsGDlike_name_geo == "China" | 
                                             GsGDlike_table$GsGDlike_name_geo == "Hong_Kong"), ],
              
              aes(x = GsGDlike_name_year, y = Divergence), color = "gray", size = 3, alpha = 0.2, stroke = 0) + 
  
  # 2.3.4
  geom_point( data = GsGDlike_table[ which( ( GsGDlike_table$GsGDlike_name_geo   == "China" | 
                                              GsGDlike_table$GsGDlike_name_geo   == "Hong_Kong") &
                                              GsGDlike_table$GsGDlike_clade      == "2.3.4" & 
                                              GsGDlike_table$GsGDlike_name_sero  == "H5N1" ), ],
              
              aes(x = GsGDlike_name_year, y = Divergence), color = "#d62728", size = 3, alpha = 0.8, stroke = 0) +
  
  # 2.3.2
  geom_point( data = GsGDlike_table[ which( ( GsGDlike_table$GsGDlike_name_geo   == "China" | 
                                              GsGDlike_table$GsGDlike_name_geo   == "Hong_Kong") &
                                              GsGDlike_table$GsGDlike_clade      == "2.3.2" ), ],
              
              aes(x = GsGDlike_name_year, y = Divergence), color = "#2ca02c", size = 3, alpha = 0.4, stroke = 0) + 
  
  # 2.2
  geom_point( data = GsGDlike_table[ which( ( GsGDlike_table$GsGDlike_name_geo   == "China" | 
                                              GsGDlike_table$GsGDlike_name_geo   == "Hong_Kong") &
                                              GsGDlike_table$GsGDlike_clade      == "2.2" ), ],
              
              aes(x = GsGDlike_name_year, y = Divergence), color = "#ff7f0e", size = 3, alpha = 0.4, stroke = 0) 



### 2344 H5Nx -------------------------------- 

subtreseq( list_filedir = "./Nx/taxalist_2344.txt", seq_filedir = "./Tree/allh5_GsGDlike_6322.fasta")

rmDup( fasfile = "./Nx/c2344_1365.fasta", 
       year    = c(1000, 2012) )

# system("~/Raxmldata/raxml_AVX2 -f a -p 666 -s /Users/yaosmacbook/Desktop/data_souce/Nx/c2344_21.fasta -x 616 -#autoMRE -m GTRGAMMA --HKY85 -n h5nx")
# check in TempEst & taxaInfo(useTree = TRUE, file = "./Nx/c2344_21_e.tre", root2tip = TRUE)
# I remove : KP735813_goose_Yangzhou_ZG60_2009_|China|_H5N5_2009.929
# II keep 8 as previous 





### hyphy output  --------------------------------

# FORM: clade(+host)_geo_time_gene


## 1007 ( NEED re-do ) ----------------

# t.hyphy1007a = c( "c234a_SW_T", "c234a_nSW_T", "c234a_SEA_T", "c232a_SW_T", "c232a_nSW_T", "c232a_SEA_T" ) 
ac.hyphy1007a <- list(t.hyphy1007a)

c234.nonML <- which( faslist.234[[9]] == "nonML" )
c234.nSW   <- setdiff( which( faslist.234[[2]] == "China"| faslist.234[[2]] == "Hong_Kong" ), which( faslist.234[[8]] == "cnSW") )
ac.hyphy1007a[[2]] <- faslist.234[[1]][ intersect( c234.nonML, which( faslist.234[[8]] == "cnSW") ) ] 
ac.hyphy1007a[[3]] <- faslist.234[[1]][ intersect( c234.nonML, c234.nSW ) ] 
ac.hyphy1007a[[4]] <- faslist.234[[1]][ intersect( c234.nonML, which( faslist.234[[8]] == "SEA") ) ] 

c232.nonML <- which( faslist.232[[9]] == "nonML" )
c232.nSW   <- setdiff( which( faslist.232[[2]] == "China"| faslist.232[[2]] == "Hong_Kong" ), which( faslist.232[[8]] == "cnSW") )
ac.hyphy1007a[[5]] <- faslist.232[[1]][ intersect( c232.nonML, which( faslist.232[[8]] == "cnSW") ) ] 
ac.hyphy1007a[[6]] <- faslist.232[[1]][ intersect( c232.nonML, c232.nSW ) ] 
ac.hyphy1007a[[7]] <- faslist.232[[1]][ intersect( c232.nonML, which( faslist.232[[8]] == "SEA") ) ] 


# t.hyphy1007b = c( "c234a_CN_0405", "c234a_CN_0607", "c234a_CN_0809", "c234a_CN_1011", "c232a_CN_0405", "c232a_CN_0607", "c232a_CN_0809", "c232a_CN_1011" )
ac.hyphy1007b <- list(t.hyphy1007b)

c234.nonML_CN <- intersect( which( faslist.234[[9]] == "nonML" ), which( faslist.234[[2]] == "China"| faslist.234[[2]] == "Hong_Kong" ) ) 
ac.hyphy1007b[[2]] <- faslist.234[[1]][ intersect( c234.nonML, which( floor( faslist.234[[4]] ) == 2004 | floor( faslist.234[[4]] ) == 2005 )  ) ] 
ac.hyphy1007b[[3]] <- faslist.234[[1]][ intersect( c234.nonML, which( floor( faslist.234[[4]] ) == 2006 | floor( faslist.234[[4]] ) == 2007 )  ) ] 
ac.hyphy1007b[[4]] <- faslist.234[[1]][ intersect( c234.nonML, which( floor( faslist.234[[4]] ) == 2008 | floor( faslist.234[[4]] ) == 2009 )  )  ] 
ac.hyphy1007b[[5]] <- faslist.234[[1]][ intersect( c234.nonML, which( floor( faslist.234[[4]] ) == 2010 | floor( faslist.234[[4]] ) == 2011 )  )  ] 

c232.nonML_CN <- intersect( which( faslist.232[[9]] == "nonML" ), which( faslist.232[[2]] == "China"| faslist.232[[2]] == "Hong_Kong" ) ) 
ac.hyphy1007b[[6]] <- faslist.232[[1]][ intersect( c232.nonML, which( floor( faslist.232[[4]] ) == 2004 | floor( faslist.232[[4]] ) == 2005 )  ) ] 
ac.hyphy1007b[[7]] <- faslist.232[[1]][ intersect( c232.nonML, which( floor( faslist.232[[4]] ) == 2006 | floor( faslist.232[[4]] ) == 2007 )  ) ] 
ac.hyphy1007b[[8]] <- faslist.232[[1]][ intersect( c232.nonML, which( floor( faslist.232[[4]] ) == 2008 | floor( faslist.232[[4]] ) == 2009 )  )  ] 
ac.hyphy1007b[[9]] <- faslist.232[[1]][ intersect( c232.nonML, which( floor( faslist.232[[4]] ) == 2010 | floor( faslist.232[[4]] ) == 2011 )  )  ] 


# t.hyphy1007c = c( "c234a_SW_0405", "c234a_SW_06", "c234a_SW_07", "c234a_SW_0809", "c234a_nSW_0405", "c234a_nSW_06", "c234a_nSW_07", "c234a_nSW_0809" )

ac.hyphy1007c <- list(t.hyphy1007c)

c234.nonML_SW  <- intersect( which( faslist.234[[9]] == "nonML" ),  which( faslist.234[[8]] == "cnSW") ) 
c234.nonML_nSW <- intersect( which( faslist.234[[9]] == "nonML" ), setdiff( which( faslist.234[[2]] == "China"| faslist.234[[2]] == "Hong_Kong" ), which( faslist.234[[8]] == "cnSW" ) )  ) 

ac.hyphy1007c[[2]] <- faslist.234[[1]][ intersect( c234.nonML_SW, which( floor( faslist.234[[4]] ) == 2004 | floor( faslist.234[[4]] ) == 2005 )  ) ] 
ac.hyphy1007c[[3]] <- faslist.234[[1]][ intersect( c234.nonML_SW, which( floor( faslist.234[[4]] ) == 2006 )  ) ] 
ac.hyphy1007c[[4]] <- faslist.234[[1]][ intersect( c234.nonML_SW, which( floor( faslist.234[[4]] ) == 2007 )  )  ] 
ac.hyphy1007c[[5]] <- faslist.234[[1]][ intersect( c234.nonML_SW, which( floor( faslist.234[[4]] ) == 2008 | floor( faslist.234[[4]] ) == 2009 )  )  ] 

ac.hyphy1007c[[6]] <- faslist.234[[1]][ intersect( c234.nonML_nSW, which( floor( faslist.234[[4]] ) == 2004 | floor( faslist.234[[4]] ) == 2005 )  ) ] 
ac.hyphy1007c[[7]] <- faslist.234[[1]][ intersect( c234.nonML_nSW, which( floor( faslist.234[[4]] ) == 2006 )  ) ] 
ac.hyphy1007c[[8]] <- faslist.234[[1]][ intersect( c234.nonML_nSW, which( floor( faslist.234[[4]] ) == 2007 )  )  ] 
ac.hyphy1007c[[9]] <- faslist.234[[1]][ intersect( c234.nonML_nSW, which( floor( faslist.234[[4]] ) == 2008 | floor( faslist.234[[4]] ) == 2009 )  )  ] 



for(i in 1: 3)
{
  aclist <- c("ac.hyphy1007a", "ac.hyphy1007b", "ac.hyphy1007c")
  
  aclist.i <- get( aclist[i] )
  
  for(k in 2: length(aclist.i) )
  {
    subfastaSeq( AC = TRUE, AC_list = aclist.i[[k]], filedir = "Tree/allh5_GsGDlike_6322.fasta", no = aclist.i[[1]][k-1] )
    
  }
}


## 1019 ----------------


# t.hyphy1019 = c( "c232a_CN_0406", "c232a_CN_0711", "c234a_nSW_0406", "c234a_nSW_0711", "c234a_SW_0406", "c234a_SW_0711" ) 

ac.hyphy1019 <- list(t.hyphy1019)

c232.nonML_CN <- intersect( which( faslist.232[[9]] == "nonML" ), which( faslist.232[[2]] == "China"| faslist.232[[2]] == "Hong_Kong" ) ) 
ac.hyphy1019[[2]] <- faslist.232[[1]][ intersect( c232.nonML_CN, which( floor( faslist.232[[4]] ) == 2004 | floor( faslist.232[[4]] ) == 2005 | floor( faslist.232[[4]] ) == 2006 )  ) ] 
ac.hyphy1019[[3]] <- faslist.232[[1]][ intersect( c232.nonML_CN, which( floor( faslist.232[[4]] ) == 2007 | floor( faslist.232[[4]] ) == 2008 | floor( faslist.232[[4]] ) == 2009 | floor( faslist.232[[4]] ) == 2010 | floor( faslist.232[[4]] ) == 2011 )  ) ] 

c234.nonML_nSW <- intersect( which( faslist.234[[9]] == "nonML" ), setdiff( which( faslist.234[[2]] == "China"| faslist.234[[2]] == "Hong_Kong" ), which( faslist.234[[8]] == "cnSW" ) )  ) 
c234.nonML_SW  <- intersect( which( faslist.234[[9]] == "nonML" ),  which( faslist.234[[8]] == "cnSW") ) 

ac.hyphy1019[[4]] <- faslist.234[[1]][ intersect( c234.nonML_nSW, which( floor( faslist.234[[4]] ) == 2004 | floor( faslist.234[[4]] ) == 2005 | floor( faslist.234[[4]] ) == 2006 )  ) ] 
ac.hyphy1019[[5]] <- faslist.234[[1]][ intersect( c234.nonML_nSW, which( floor( faslist.234[[4]] ) == 2007 | floor( faslist.234[[4]] ) == 2008 | floor( faslist.234[[4]] ) == 2009 | floor( faslist.234[[4]] ) == 2010 | floor( faslist.234[[4]] ) == 2011 )  ) ] 

ac.hyphy1019[[6]] <- faslist.234[[1]][ intersect( c234.nonML_SW, which( floor( faslist.234[[4]] ) == 2004 | floor( faslist.234[[4]] ) == 2005 | floor( faslist.234[[4]] ) == 2006 )  ) ] 
ac.hyphy1019[[7]] <- faslist.234[[1]][ intersect( c234.nonML_SW, which( floor( faslist.234[[4]] ) == 2007 | floor( faslist.234[[4]] ) == 2008 | floor( faslist.234[[4]] ) == 2009 | floor( faslist.234[[4]] ) == 2010 | floor( faslist.234[[4]] ) == 2011 )  ) ] 



for(k in 2: length(ac.hyphy1019) )
{
  subfastaSeq( AC = TRUE, AC_list = ac.hyphy1019[[k]], filedir = "Tree/allh5_GsGDlike_6322.fasta", no = ac.hyphy1019[[1]][k-1] )
  
}

pat_sample <- paste0( "~/Desktop/data_souce/Clade_allh5/hyphy1019/HA/", 
                      list.files( "~/Desktop/data_souce/Clade_allh5/hyphy1019/HA/" ))

for(k in 1: length(pat_sample) )
{
  ntpartition( position = c(1: 1017), filedir = pat_sample[k], no = "-HA1.fasta")
  ntpartition( position = c(1018: 1683), filedir = pat_sample[k], no = "-HA2.fasta")
}







