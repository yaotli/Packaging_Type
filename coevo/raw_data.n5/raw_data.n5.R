source( "./function.coevo.R" )


# database download ------------
#
# GISAID -> 35 isolates ( 35 in the file)
# 1. type 'A'; 2. N '5' 
# 3. host undo 'Lab derived, Unknown, NA'
# 4. location 'all'
# 5. required seqments 'HA, NA'
# 6. mim length '1000'; 7. only GISAID 
# format 'Isolate name Type Collection date Isolate ID'
# 
# NCBI -> 589 isolates
# 1. type 'A'; 2. host undo 'Unknown' 
# 3. country 'any'; 4. segment 'HA'
# 5. H 'any'; N '6' 
# 6. required segments 'HA, NA' 
# 7. 'exclude' pH1N1, 'exclude' lab strains
#    'include' lineage defining strains 
#    'include' FLU project, 
#    'include' vaccine strains, 'exclude' mixed 
# format '>{accession}_{strain}_{serotype}_|{country}|_{year}-{month}-{day}'
# 


fas_n5_g <- "./raw_data.n5/sources/pN5_G_35_20190115.fasta"
csv_g    <- "./raw_data.n5/sources/G_20190115.csv"
fas_n5_n <- "./raw_data.n5/sources/pN5_N_589_20190115.fasta"

fas_curated_h5 <- "./raw_data/processed/pH5_8334.fasta"


# data cleaning  ------------

# GISAID 
# to remove one isolate which appears twice 
#

seq_n5_g    <- fastaEx( fas_n5_g )$seq
id_n5_g     <- fastaEx( fas_n5_g )$id

infols_n5_g <- idInfo.na( rawid = id_n5_g, datasource = "g", g.csv = csv_g, na_subtype = 5 )

# NCBI
# 

seq_n5_n    <- fastaEx( fas_n5_n )$seq
id_n5_n     <- fastaEx( fas_n5_n )$id
infols_n5_n <- idInfo.na( rawid = id_n5_n, datasource = "n", na_subtype = 6 )


# combine and remove the duplicate isolate ------------

seq_n5    <- c( seq_n5_g, seq_n5_n )
infols_n5 <- lapply( as.list( seq( 1, length( infols_n5_n ) ) ), function(x){ c( infols_n5_g[[x]], infols_n5_n[[x]])  }  )

infols_n5[[ length(infols_n5) + 1 ]] = seq_n5

s_infols_n5      <- strainSelect( infols_n5 )
s_infols_n5[[4]] <- seqDate( s_infols_n5[[4]] )

c_infols_n5 <- seqSelect( minlth = 1200, maxamb = 1, s_infols_n5, rmdup = FALSE ) # n = 567


# examine HA-NA information ------------

info_h5 <- taxaInfo( fas_curated_h5 )

h5.i <- grep( "H5", c_infols_n5[[2]] )
m.i  <- match( c_infols_n5[[5]][h5.i], info_h5[[5]] )


dismatch = c()
for( i in 1: length( m.i ) )
{
  if( !is.na(m.i[i]) )
  {
    if( c_infols_n5[[3]][h5.i][i] !=  info_h5[[2]][m.i][i] ){ dismatch = c(dismatch, i) }
    if( abs( as.numeric( c_infols_n5[[4]][h5.i][i] ) - info_h5[[4]][m.i][i] ) > 0.1 ){ dismatch = c(dismatch, i) }
    
  }
}

# export ------------

write.fasta( sequences = c_infols_n5[[6]], 
             names     = paste0( c_infols_n5[[1]], "_",
                                 c_infols_n5[[5]], "_|",
                                 c_infols_n5[[3]], "|_",
                                 c_infols_n5[[2]], "_",
                                 c_infols_n5[[4]]
             ), file.out  = paste0( "raw_data.n5/processed/pn5_", length( c_infols_n5[[1]] ) ,".fasta") )


# seq manipulation ------------

trimtool( propblank = 0.9, filedir = "./raw_data.n5/processed/pN5_567_align.fasta" )

# remove outliers ------------

# two outliers found by TempEst
# EPI24750_Duck_Hunnan_70_2004_|China|_H5N5_2004.497
# MH597119_ruddy_turnstone_Delaware_1016391_2003_|USA|_H9N5_2003.384

n567.seqname = fastaEx( "./raw_data.n5/processed/pN5_567_trim2.fasta" )$id
n567.seq     = fastaEx( "./raw_data.n5/processed/pN5_567_trim2.fasta" )$seq

rm <- grep( "EPI24750|MH597119", n567.seqname )
write.fasta( n567.seq[-rm], n567.seqname[-rm], file.out = "./raw_data.n5/processed/pN5_565.fasta" )






