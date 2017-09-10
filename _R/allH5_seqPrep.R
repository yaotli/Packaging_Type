library(seqinr)
library(stringr)
library(ggtree)
library(ape)


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

allh5_g_seq <- keepLongSeq( fastaEx( fas_allh5_g )$seq, fastaEx( fas_allh5_g )$id )$seq
allh5_g_id  <- keepLongSeq( fastaEx( fas_allh5_g )$seq, fastaEx( fas_allh5_g )$id )$id

a.string.g  <- "EPI_ISL_([0-9]+)"
s.string.g  <- "_A_/_(H5N[0-9a-zA-Z]{1,2})_"
y.string.g  <- "_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)" 
n.Nstring.g <- "_EPI_ISL_([0-9]+)|_A_/_H5N[0-9a-zA-Z]{1,2}|_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)"


# accession number 
allh5_g_id_ac   <- gsub("_ISL_", "", 
                        str_match(allh5_g_id, a.string.g)[, 1] )

# serotype
allh5_g_id_sero <- str_match( pattern = s.string.g, 
                              string = allh5_g_id)[,2] 

# time
allh5_g_id_year <- str_match( pattern = y.string.g, 
                              string = allh5_g_id) 

allh5_g_id_year <- gsub( pattern = "^_", replacement = "", x = allh5_g_id_year)


# based on the 'Note' column in csv 
allh5_g_id_year[ grep("178249", allh5_g_id_ac) ] = gsub( pattern     = "_\\(Month_and_day_unknown\\)", 
                                                         replacement = "-12_(Day_unknown)",
                                                         x           = allh5_g_id_year[ grep("178249", allh5_g_id_ac) ] )

# strain name
allh5_g_id_name <- gsub( pattern = n.Nstring.g, 
                         x = allh5_g_id, replacement = "")

allh5_g_id_name[ which( startsWith(allh5_g_id_name, "A/") == FALSE) ] <- 
  gsub(pattern     = "_A/", 
       replacement = "A/", 
       x           = allh5_g_id_name[ which( startsWith(allh5_g_id_name, "A/") == FALSE) ])

allh5_g_id_name_c  <- gsub("\\(|\\)|\\[|\\]|\\.|:|-|/", "_", allh5_g_id_name)
allh5_g_id_name_c  <- gsub("__", "_", allh5_g_id_name_c)
allh5_g_id_name_c  <- gsub("\\'|\\?|>", "", allh5_g_id_name_c)
allh5_g_id_name_c  <- gsub("A_", "", allh5_g_id_name_c)
allh5_g_id_name_c  <- gsub("_$", "", allh5_g_id_name_c)


# geo 
g_geo <- gsub( pattern = " ", replacement = "_", 
               read.csv( csv_allh5_g, header = TRUE, stringsAsFactors = FALSE)$Location )

g_geo <- gsub("_$", "", 
              str_match( g_geo, "([A-Za-z_]+)_/_([A-Za-z_]+)" )[,3] ) 

g_geo[ which( is.na(g_geo) == TRUE ) ] = "Unknown"
g_geo[ which(g_geo == "Russian_Federation") ] = "Russia"
g_geo[ which(g_geo == "United_States") ] = "USA"
g_geo[ which(g_geo == "Korea") ] = "South_Korea"

allh5_g_geo <- 
g_geo[
match( allh5_g_id_ac, 
       gsub("_ISL_", "", read.csv( csv_allh5_g, header = TRUE, stringsAsFactors = FALSE)$Isolate_Id ) ) ]


# TRUE %in% is.na( c(allh5_g_id_ac, allh5_g_id_sero, allh5_g_id_year, allh5_g_id_name, g_geo) )


## NCBI ----------------
# n = 7363

allh5_n_seq <- fastaEx( fas_allh5_n )$seq
allh5_n_id  <- fastaEx( fas_allh5_n )$id

a.string.n   <- "[A-Z]{1,2}[0-9]{5,6}"
s.g.string.n <- "_(H5[N0-9]{0,2})_\\|([a-zA-Z_\\']+)\\|"
y.string.n   <- "_[0-9]{4}-[0-9]{2}-[0-9]{2}|_[0-9]{4}-[0-9]{2}-|_[0-9]{4}--|_--"
n.Nstring.g  <- "[A-Z]{1,2}[0-9]{5,6}_|_(H5[N0-9]{0,2})_\\|([a-zA-Z_\\']+)\\||_([0-9]{4}[-0-9]{2,6})"

# accession number
allh5_n_id_ac   <- str_match( string = allh5_n_id, 
                              pattern = a.string.n)[,1]

# serotype
allh5_n_id_sero <- str_match( string = allh5_n_id, 
                              pattern = s.g.string.n)[,2]
allh5_n_id_sero[ which(allh5_n_id_sero == "H5")  ] = "H5N0"

# geo
allh5_n_id_geo  <- str_match( string = allh5_n_id, 
                              pattern = s.g.string.n)[,3]

allh5_n_id_geo[ which(allh5_n_id_geo == "Viet_Nam") ] = "Vietnam"
allh5_n_id_geo[ which(allh5_n_id_geo == "Cote_d'Ivoire") ] = "Cote_dIvoire"


# time 
allh5_n_id_year <- str_match( string = allh5_n_id, 
                              pattern = y.string.n)

allh5_n_id_year <- gsub( pattern = "_--", 
                         replacement = "1900-01-01", allh5_n_id_year)

allh5_n_id_year <- gsub( pattern = "^_", 
                         replacement = "", allh5_n_id_year)
                         

# strain name 
allh5_n_id_name <- gsub( pattern = n.Nstring.g, x = allh5_n_id, 
                         replacement = "")

allh5_n_id_name[ which( startsWith(allh5_n_id_name, "A/") == FALSE) ] <- 
  paste0("A/", allh5_n_id_name[ which( startsWith(allh5_n_id_name, "A/") == FALSE) ])

allh5_n_id_name_c  <- gsub("\\(|\\)|\\[|\\]|\\.|:|-|/|__", "_", allh5_n_id_name)
allh5_n_id_name_c  <- gsub("__", "_", allh5_n_id_name_c)
allh5_n_id_name_c  <- gsub("\\'|\\?|>", "", allh5_n_id_name_c)
allh5_n_id_name_c  <- gsub("A_", "", allh5_n_id_name_c)
allh5_n_id_name_c  <- gsub("_$", "", allh5_n_id_name_c)


# TRUE %in% is.na( c( allh5_n_id_ac, allh5_n_id_sero, allh5_n_id_geo, allh5_n_id_year, allh5_n_id_name_c ) )


## deal with duplicated strain name ----------------

allh5_seq  <- c( allh5_n_seq, allh5_g_seq )

allh5_name <- c( allh5_n_id_name_c, allh5_g_id_name_c )
allh5_ac   <- c( allh5_n_id_ac, allh5_g_id_ac )
allh5_sero <- c( allh5_n_id_sero, allh5_g_id_sero )
allh5_year <- c( allh5_n_id_year, allh5_g_id_year )
allh5_geo  <- c( allh5_n_id_geo, allh5_g_geo )

# inspect the dupplicated strain name
d.n = allh5_name[ duplicated( allh5_name ) | duplicated( allh5_name, fromLast = TRUE ) ] 
d.y = allh5_year[ duplicated( allh5_name ) | duplicated( allh5_name, fromLast = TRUE ) ] 
d.a = allh5_ac[ duplicated( allh5_name ) | duplicated( allh5_name, fromLast = TRUE ) ] 
d.l = sapply( allh5_seq[ duplicated( allh5_name ) | duplicated( allh5_name, fromLast = TRUE ) ],
              function(x)
              {
                y = c2s(x)
                z = gsub("-|~", "", y)
                z = grep( pattern = "a|t|c|g", x = y, ignore.case = TRUE, value = TRUE)
                
                l = length( s2c(z) )
                return(l)
                
              } )
  
du = data.frame(d.n, d.y, d.a, d.l)

# do the selection

toberemove <- c()
dup        <- which( duplicated(allh5_name) )

for(i in 1: length( dup ) ) 
{
  id_dup_ii <- which( allh5_name %in% allh5_name[ dup[i] ] == TRUE )
  
  
  lth_ii    <- sapply( allh5_seq[id_dup_ii],
                       
                       function(x)
                        {
                          y = c2s(x)
                          z = gsub("-|~", "", y)
                          z = grep( pattern = "a|t|c|g", 
                                    x = y, 
                                    ignore.case = TRUE, value = TRUE )
                          
                          l = length( s2c(z) )
                          
                          return(l)
                          
                        } )
  
  SeqL      <- which.max( lth_ii )   
  
  
  if ( length(  which( lth_ii == max(lth_ii) )  ) > 1 )
  {
    
    id_dup_jj  <- id_dup_ii[ which( lth_ii == max(lth_ii) ) ] 
    
    nchar_jj   <- nchar( gsub( pattern     = "[-\\(\\)A-Za-z]+", 
                               replacement = "",
                               x           = allh5_year[id_dup_jj] ) )

    id_dup_j   <- which.max( nchar_jj )
    SeqL       <- which( lth_ii == max(lth_ii) )[id_dup_j]
    
    
    if (  length( which( nchar_jj == max(nchar_jj) ) ) > 1  )
      {
      
        id_dup_kk <- id_dup_jj[ which(nchar_jj == max(nchar_jj) ) ]
        
        ac        <- allh5_ac[id_dup_kk]
        ac.a      <- nchar( gsub( pattern = "[0-9]+", replacement = "", x = ac) )
        ac.d      <- as.numeric( gsub( pattern = "[a-zA-Z]+", replacement = "", x = ac ) )
        
        ac.df     <- data.frame( id_dup_kk, ac.a, ac.d )
        
        SeqL      <- which( id_dup_ii == ac.df[order( ac.df[,2], ac.df[,3] ),][1,1] )
      
      }
  }
  
  toberemove = c(toberemove, id_dup_ii[-SeqL] )
  
  print( allh5_ac[ id_dup_ii[SeqL] ] )
}

remain <- seq( 1, length(allh5_seq) )[- toberemove]


c_allh5_seq  <- allh5_seq[remain]
c_allh5_name <- allh5_name[remain]
c_allh5_ac   <- allh5_ac[remain]
c_allh5_sero <- allh5_sero[remain]
c_allh5_year <- allh5_year[remain]
c_allh5_geo  <- allh5_geo[remain]
 
c_allh5_info <- ifelse( grepl( pattern = "--$|Month", x = c_allh5_year ), 1, 0)


# transform year format

c_allh5_year <- 
gsub( pattern     = "_\\(Day_unknown\\)", 
      replacement = "-15", 
      x           = gsub( pattern     = "_\\(Month_and_day_unknown\\)", 
                          replacement = "-07-01", 
                          x           = c_allh5_year ) )
c_allh5_year <-
gsub( pattern     = "-$", 
      replacement = "-15", 
      x           = gsub( pattern = "--$", 
                      replacement = "-07-01", 
                      x           = c_allh5_year)   ) 


          
d_allh5_year <- phylo_date(c_allh5_year) 

# remove seq without time info.

c_allh5      <- data.frame( c_allh5_ac, c_allh5_name, c_allh5_geo, c_allh5_sero, d_allh5_year,
                            c_allh5_info, 
                            stringsAsFactors = FALSE ) 

d_allh5      <- c_allh5[-1, ]
d_allh5_seq  <- c_allh5_seq[ seq(1, length(c_allh5_seq) )[-1] ]

o_allh5      <- cbind( d_allh5, no = seq(1, length(d_allh5_seq) ) )
o_allh5      <- o_allh5[ order( o_allh5[,5], o_allh5[,1], o_allh5[,4]), ]


write.fasta( sequences = d_allh5_seq[o_allh5$no], 
                 names = paste(o_allh5$c_allh5_ac, 
                               o_allh5$c_allh5_name,
                               o_allh5$c_allh5_geo,
                               o_allh5$c_allh5_sero,
                               o_allh5$d_allh5_year, sep = "_"), 
             
              file.out = "./o_allh5_9135.fasta")


o_allh5_seq  <- d_allh5_seq[o_allh5$no]
o_allh5_id   <- paste(o_allh5$c_allh5_ac, 
                      o_allh5$c_allh5_name,
                      o_allh5$c_allh5_geo,
                      o_allh5$c_allh5_sero,
                      o_allh5$d_allh5_year, sep = "_")


## curation sequences ----------------

# remove sequence with > 1 % ambiguous nucleotides
# remove sequence < 1500 nt 
# remove duplicated sequences 

ATCG = c("a", "t", "c", "g")

i_allh5_id  <- o_allh5_id[ sort(o_allh5$c_allh5_info, index.return = TRUE )$ix ]
i_allh5_seq <- o_allh5_seq[ sort(o_allh5$c_allh5_info, index.return = TRUE )$ix ]

# n = 558 
lth_amb     <- which( sapply( i_allh5_seq, 
                              function(x)
                              {
                                x.s   <- gsub( "-|~", "", c2s( x ) )
                                x.l   <- length( s2c(x.s) )
                                
                                x.c.s <- grep( "a|t|c|g", s2c(x.s) )[1]
                                x.c.e <- grep( "a|t|c|g", s2c(x.s) )[ length( grep( "a|t|c|g", s2c(x.s) ) ) ]
                                
                                x.a   <- length( which(! s2c(x.s)[x.c.s: x.c.e] %in% ATCG ) )
                       
                                return( x.l < 1500 | x.a > 0.01*x.l )
                                
                              })  
                      == TRUE ) 


# n = 1273
dup          <- which( duplicated( sapply( i_allh5_seq,  
                                           function(x)
                                           {
                                            x.s <- gsub( "~|-", "", c2s(x) )
                                             
                                            return(x.s) 
                                           }
                                           ) ) == TRUE  )

# remain = 7368

remain_s    <- seq(1, length(i_allh5_id) )[ - unique( sort( c(lth_amb, dup) ) ) ]

u_allh5_id  <- i_allh5_id[remain_s]
u_allh5_seq <- i_allh5_seq[remain_s]

write.fasta( sequences = u_allh5_seq, 
             names     = u_allh5_id,
             file.out  = "./curated_allh5_7368.fasta")



### tree and clade annotation --------------------------------

system("mafft curated_allh5_7368.fasta > align_allh5_7368.fasta")

system("mkdir aligntrim_allh5")
system("mv align_allh5_7368.fasta aligntrim_allh5/align_allh5_7368.fasta")

# manuall trim in BioEdit (lth = 1677)
# Fasttree

system("~/./FastTree -nt -nni 10 -spr 4 -gtr -cat 20 -gamma -notop <./aligntrim_allh5/trim_allh5_7368_lth1677.fasta> ./tree/allh5_7368.tre")

# manual extract GsGD 








