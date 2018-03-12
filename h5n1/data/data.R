source("functions.R")

require(seqinr)
require(stringr)
require(lubridate)



# read-in
fas_ha_g <- "data/raw/pH5_G_2136_20171124.fasta" 
fas_na_g <- "data/raw/pNA_G_2135_20171124.fasta"
csv_g    <- "data/raw/G_2135_20171124.csv"

fas_ha_n <- "data/raw/pH5_N_5753_20171124.fasta" 
fas_na_n <- "data/raw/pNA_N_5721_20171124.fasta"


# remove erroneous sequence (GISAID)
seq.ha.g <- keepLongSeq( fastaEx(fas_ha_g)$seq, fastaEx(fas_ha_g)$id )$seq
id.ha.g  <- keepLongSeq( fastaEx(fas_ha_g)$seq, fastaEx(fas_ha_g)$id )$id

seq.na.g <- fastaEx(fas_na_g)$seq
id.na.g  <- fastaEx(fas_na_g)$id

# parse the info. in seqence names
infolist.ha.g <- idInfo( rawid = id.ha.g, datasource = "g", g.csv = csv_g)
infolist.na.g <- idInfo( rawid = id.na.g, datasource = "g", g.csv = csv_g)

# according to the note column 
infolist.ha.g[[4]][grep("178249", infolist.ha.g[[1]])] = "2014-12_(Day_unknown)"
infolist.na.g[[4]][grep("178249", infolist.na.g[[1]])] = "2014-12_(Day_unknown)"


# discard strain without NA (NCBI)
id.ha.n <- fastaEx( fas_ha_n )$id
id.na.n <- fastaEx( fas_na_n )$id

remain = seq(1, length(id.ha.n) )[- grep( "yellow_billed_teal_Chile_14_2014", idInfo(id.ha.n)[[5]] )]

id.ha.n  <- fastaEx( fas_ha_n )$id[ remain ]
seq.ha.n <- fastaEx( fas_ha_n )$seq[ remain ]
seq.na.n <- fastaEx( fas_na_n )$seq

infolist.ha.n <- idInfo( id.ha.n )
infolist.na.n <- idInfo( id.na.n )


# combine 
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


# remove duplicated *strain*
s_infolist.ha <- strainSelect( infolist.ha )
s_infolist.na <- strainSelect( infolist.na )

s_infolist.ha[[4]] <- seqDate( s_infolist.ha[[4]] )
s_infolist.na[[4]] <- seqDate( s_infolist.na[[4]] )


# curate sequences
c_infolist.ha <- seqSelect( minlth = 1500, maxamb = 1, s_infolist.ha, rmdup = FALSE )
c_infolist.na <- seqSelect( minlth = 1200, maxamb = 1, s_infolist.na, rmdup = FALSE )


# merge HA-NA info.
int_hana <- intersect( c_infolist.ha[[5]], c_infolist.na[[5]] )

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

c_infolist.ha[[3]][ match(int_hana, c_infolist.ha[[5]]) ][conf][3] = "China"
int_idx.ha <- match( int_hana, c_infolist.ha[[5]] )
int_idx.na <- match( int_hana, c_infolist.na[[5]] )  


paired_infolist_ha <- list()
for( l in 1: length(c_infolist.ha) )
{
  paired_infolist_ha[[ l ]] = c_infolist.ha[[l]][int_idx.ha]
}

paired_infolist_na <- list()
for( l in 1: length(c_infolist.na) )
{
  paired_infolist_na[[ l ]] = c_infolist.na[[l]][int_idx.na]
}

paired_infolist_na[[3]] = paired_infolist_ha[[3]]
paired_infolist_na[[4]] = paired_infolist_ha[[4]]


# export
write.csv( file = "data/pH5NA.csv", x =
             data.frame( ac.ha = paired_infolist_ha[[1]], 
                         name  = paired_infolist_ha[[5]],
                         geo   = paired_infolist_ha[[3]],
                         sero  = paired_infolist_ha[[2]],
                         year  = paired_infolist_ha[[4]],
                         ac.na = paired_infolist_na[[1]], 
                         stringsAsFactors = FALSE), row.names = FALSE)

write.fasta( sequences = paired_infolist_na[[6]], 
             names     = paste0( paired_infolist_na[[1]], "_",
                                 paired_infolist_na[[5]], "_|",
                                 paired_infolist_na[[3]], "|_",
                                 paired_infolist_na[[2]], "_",
                                 paired_infolist_na[[4]]
             ), file.out  = paste0( "data/pNA_", length( paired_infolist_na[[1]] ) ,".fasta") )

write.fasta( sequences = c_infolist.ha[[6]], 
             names     = paste0( c_infolist.ha[[1]], "_",
                                 c_infolist.ha[[5]], "_|",
                                 c_infolist.ha[[3]], "|_",
                                 c_infolist.ha[[2]], "_",
                                 c_infolist.ha[[4]]
             ), file.out  = paste0( "data//pH5_", length( c_infolist.ha[[1]] ) ,".fasta")  )


# extract N1 
rmDup( fasfile = "data/pNA_7138.fasta", sero = "H5N1", rmdup = FALSE)

# trim
trimtool( propblank = 0.9, filedir = "data/processed/align_pN1_4468.fasta" )





