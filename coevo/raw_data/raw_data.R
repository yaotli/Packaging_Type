source( "./function.coevo.R" )


# database download ------------
#
# GISAID -> 2689 isolates 
# 1. type 'A'; 2. H '5' 
# 3. host undo 'Lab derived, Unknown, NA'
# 4. location 'all'
# 5. required seqments 'HA, NA'
# 6. mim length '1000'; 7. only GISAID 
# 
# NCBI -> 6452 / 6387 
# 1. type 'A'; 2. host undo 'Unknown' 
# 3. country 'any'; 4. segment 'HA'
# 5. H '5'; N 'any' 
# 6. required segments 'HA, NA' 
# 7. 'exclude' pH1N1, 'exclude' lab strains
#    'include' lineage defining strains 
#    'include' FLU project, 
#    'include' vaccine strains, 'exclude' mixed 
#

fas_ha_g <- "./raw_data/sources/pH5_G_2690_20180912.fasta"
fas_na_g <- "./raw_data/sources/pNA_G_2689_20180912.fasta"
csv_g    <- "./raw_data/sources/G_20180912.csv"
fas_ha_n <- "./raw_data/sources/pH5_N_6452_20180912.fasta"
fas_na_n <- "./raw_data/sources/pNA_N_6387_20180912.fasta"


# data cleaning  ------------

# GISAID 
# to remove one isolate which appears twice in the 'seq_ha_g' data
#
seq_ha_g <- keepLongSeq( fastaEx( fas_ha_g )$seq, fastaEx( fas_ha_g )$id )$seq 
id_ha_g  <- keepLongSeq( fastaEx( fas_ha_g )$seq, fastaEx( fas_ha_g )$id )$id

seq_na_g <- fastaEx( fas_na_g )$seq
id_na_g  <- fastaEx( fas_na_g )$id

infols_ha_g <- idInfo( rawid = id_ha_g, datasource = "g", g.csv = csv_g )
infols_na_g <- idInfo( rawid = id_na_g, datasource = "g", g.csv = csv_g )

# check notes by 
# read.csv( csv_g, stringsAsFactors = FALSE)$Note[1529]
# &
# read.csv( csv_g, stringsAsFactors = FALSE)$Isolate_Id[1529]

infols_ha_g[[4]][ grep( "178249", infols_ha_g[[1]] ) ] = "2014-12_(Day_unknown)"
infols_na_g[[4]][ grep( "178249", infols_na_g[[1]] ) ] = "2014-12_(Day_unknown)"


# NCBI

# setdiff( unique( idInfo( fastaEx(fas_ha_n)$id )[[5]] ),  unique( idInfo( fastaEx(fas_na_n)$id )[[5]] ) )
# n = 32
# setdiff( unique( idInfo( fastaEx(fas_na_n)$id )[[5]] ),  unique( idInfo( fastaEx(fas_ha_n)$id )[[5]] ) )
# n = 0

rm <- paste0( setdiff( unique( idInfo( fastaEx(fas_ha_n)$id )[[5]] ),  
                       unique( idInfo( fastaEx(fas_na_n)$id )[[5]] ) ), 
              collapse = "|" )

# remove isolates do not contain na sequences
rm.i     <- grep( rm, idInfo( fastaEx(fas_ha_n)$id )[[5]] )
id_ha_n  <- fastaEx( fas_ha_n )$id[ -rm.i ]
seq_ha_n <- fastaEx( fas_ha_n )$seq[ -rm.i ]
seq_na_n <- fastaEx( fas_na_n )$seq

infols_ha_n = idInfo( id_ha_n )
infols_na_n = idInfo( fastaEx(fas_na_n)$id )


# combine and remove the duplicate isolate ------------

seq_ha <- c( seq_ha_n, seq_ha_g )
seq_na <- c( seq_na_n, seq_na_g )

infols_ha <- lapply( as.list( seq( 1, length(infols_ha_n) ) ), function(x){ c( infols_ha_n[[x]], infols_ha_g[[x]])  }  )
infols_na <- lapply( as.list( seq( 1, length(infols_na_n) ) ), function(x){ c( infols_na_n[[x]], infols_na_g[[x]])  }  )

infols_ha[[ length(infols_ha) + 1 ]] = seq_ha
infols_na[[ length(infols_na) + 1 ]] = seq_na

s_infols_ha <- strainSelect( infols_ha )
s_infols_na <- strainSelect( infols_na )

s_infols_ha[[4]] <- seqDate( s_infols_ha[[4]] )
s_infols_na[[4]] <- seqDate( s_infols_na[[4]] )

c_infols_ha <- seqSelect( minlth = 1500, maxamb = 1, s_infols_ha, rmdup = FALSE ) # n = 8524
c_infols_na <- seqSelect( minlth = 1200, maxamb = 1, s_infols_na, rmdup = FALSE ) # n = 8480


# merge HA-NA information ------------

int_hana <- intersect( c_infols_ha[[5]], c_infols_na[[5]] )

conf = c()
ha_n = c_infols_ha[[5]][ match( int_hana, c_infols_ha[[5]] ) ] #name
na_n = c_infols_na[[5]][ match( int_hana, c_infols_na[[5]] ) ]
ha_g = c_infols_ha[[3]][ match( int_hana, c_infols_ha[[5]] ) ] #geo
na_g = c_infols_na[[3]][ match( int_hana, c_infols_na[[5]] ) ]
ha_y = as.numeric( c_infols_ha[[4]][ match( int_hana, c_infols_ha[[5]] ) ] ) #time
na_y = as.numeric( c_infols_na[[4]][ match( int_hana, c_infols_na[[5]] ) ] )

for(i in 1: length(int_hana) )
{
  if( ha_n[i] != na_n[i] ){ conf = c(conf, i) } 
  if( ha_g[i] != na_g[i] ){ conf = c(conf, i) } 
  if( abs(ha_y[i] - na_y[i]) > 0.5 ){ conf = c(conf, i) } 
}

# c_infols_ha[[3]][ grep( paste0( int_hana[ conf ], collapse = "|" ), c_infols_ha[[5]] ) ]
# c_infols_na[[3]][ grep( paste0( int_hana[ conf ], collapse = "|" ), c_infols_na[[5]] ) ]
# c_infols_ha[[4]][ grep( paste0( int_hana[ conf ], collapse = "|" ), c_infols_ha[[5]] ) ]
# c_infols_na[[4]][ grep( paste0( int_hana[ conf ], collapse = "|" ), c_infols_na[[5]] ) ]


c_infols_ha[[3]][ match(int_hana, c_infols_ha[[5]]) ][conf][3] = 
  c_infols_na[[3]][ match(int_hana, c_infols_na[[5]]) ][conf][3]
  
int_idx.ha <- match( int_hana, c_infols_ha[[5]] )
int_idx.na <- match( int_hana, c_infols_na[[5]] )  

paired_infols_ha <- lapply( c_infols_ha, function(x) x[int_idx.ha] )
paired_infols_na <- lapply( c_infols_na, function(x) x[int_idx.na] )

paired_infols_na[[3]] = paired_infols_ha[[3]]
paired_infols_na[[4]] = paired_infols_ha[[4]]


# export ------------

write.csv( file = "raw_data/processed/pH5NA.csv", 
           x    = data.frame( ac.ha = paired_infols_ha[[1]], 
                              name  = paired_infols_ha[[5]],
                              geo   = paired_infols_ha[[3]],
                              sero  = paired_infols_ha[[2]],
                              year  = paired_infols_ha[[4]],
                              ac.na = paired_infols_na[[1]], 
                              stringsAsFactors = FALSE), row.names = FALSE )  
             
write.fasta( sequences = paired_infols_ha[[6]], 
             names     = paste0( paired_infols_ha[[1]], "_",
                                 paired_infols_ha[[5]], "_|",
                                 paired_infols_ha[[3]], "|_",
                                 paired_infols_ha[[2]], "_",
                                 paired_infols_ha[[4]]
             ), file.out  = paste0( "raw_data/processed/pH5_", length( paired_infols_ha[[1]] ) ,".fasta") )


write.fasta( sequences = paired_infols_na[[6]], 
             names     = paste0( paired_infols_na[[1]], "_",
                                 paired_infols_na[[5]], "_|",
                                 paired_infols_na[[3]], "|_",
                                 paired_infols_na[[2]], "_",
                                 paired_infols_na[[4]]
             ), file.out  = paste0( "raw_data/processed/pNA_", length( paired_infols_na[[1]] ) ,".fasta") )


# seq manipulation ------------

gapFill( "./raw_data/processed/pH5_8334_trim2.1.fasta", s.start = 1012, s.end = 1088 )


# NA tree ------------

rmDup( "./raw_data/processed/pNA_8334.fasta", sero = "H5N1", rmdup = "FALSE" )
trimtool( propblank = 0.9, filedir = "./raw_data/processed/pN1_4696_align.fasta" )












