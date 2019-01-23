source( "./function.coevo.R" )


# database download ------------
#
# GISAID -> 752 isolates ( 752 in the file)*
# 1. type 'A'; 2. N '8' 
# 3. host undo 'Lab derived, Unknown, NA'
# 4. location 'all'
# 5. required seqments 'HA, NA'
# 6. mim length '1000'; 7. only GISAID 
# format 'Isolate name Type Collection date Isolate ID'
# 
# NCBI -> 2853 isolates
# 1. type 'A'; 2. host undo 'Unknown' 
# 3. country 'any'; 4. segment 'HA'
# 5. H 'any'; N '2' 
# 6. required segments 'HA, NA' 
# 7. 'exclude' pH1N1, 'exclude' lab strains
#    'include' lineage defining strains 
#    'include' FLU project, 
#    'include' vaccine strains, 'exclude' mixed 
# format '>{accession}_{strain}_{serotype}_|{country}|_{year}-{month}-{day}'
# 


fas_n8_g <- "./raw_data.n8/sources/pN8_G_752_20190105.fasta"
csv_g    <- "./raw_data.n8/sources/G_20190105.csv"
fas_n8_n <- "./raw_data.n8/sources/pN8_N_2853_20190105.fasta"

fas_curated_h5 <- "./raw_data/processed/pH5_8334.fasta"


# data cleaning  ------------

# GISAID 
# to remove one isolate which appears twice 
#

seq_n8_g    <- fastaEx( fas_n8_g )$seq
id_n8_g     <- fastaEx( fas_n8_g )$id
id_n8_g     = gsub( "\\&|#", "", id_n8_g )
infols_n8_g <- idInfo.na( rawid = id_n8_g, datasource = "g", g.csv = csv_g, na_subtype = 8 )

# NCBI
# 

seq_n8_n    <- fastaEx( fas_n8_n )$seq
id_n8_n     <- fastaEx( fas_n8_n )$id
infols_n8_n <- idInfo.na( rawid = id_n8_n, datasource = "n", na_subtype = 8 )

# remove no 'year' or 'country' information

rm.i = which( is.na( infols_n8_n[[2]] ) )

seq_n8_n    <- seq_n8_n[ -rm.i ]
id_n8_n     <- id_n8_n[ -rm.i ]
infols_n8_n <- idInfo.na( rawid = id_n8_n, datasource = "n", na_subtype = 2 )


# combine and remove the duplicate isolate ------------

seq_n8    <- c( seq_n8_g, seq_n8_n )
infols_n8 <- lapply( as.list( seq( 1, length( infols_n8_n ) ) ), function(x){ c( infols_n8_g[[x]], infols_n8_n[[x]])  }  )

infols_n8[[ length(infols_n8) + 1 ]] = seq_n8

s_infols_n8      <- strainSelect( infols_n8 )
s_infols_n8[[4]] <- seqDate( s_infols_n8[[4]] )

c_infols_n8 <- seqSelect( minlth = 1200, maxamb = 1, s_infols_n8, rmdup = FALSE ) # n = 3471


# examine HA-NA information ------------

info_h5 <- taxaInfo( fas_curated_h5 )

h5.i <- grep( "H5", c_infols_n8[[2]] )
m.i  <- match( c_infols_n8[[5]][h5.i], info_h5[[5]] )


dismatch = c()
for( i in 1: length( m.i ) )
{
  if( !is.na(m.i[i]) )
  {
    if( c_infols_n8[[3]][h5.i][i] !=  info_h5[[2]][m.i][i] ){ dismatch = c(dismatch, i) }
    if( abs( as.numeric( c_infols_n8[[4]][h5.i][i] ) - info_h5[[4]][m.i][i] ) > 0.1 ){ dismatch = c(dismatch, i) }
    
  }
}

c_infols_n8[[4]][h5.i][ dismatch ] = info_h5[[4]][ m.i[ dismatch ] ]

    
# export ------------

write.fasta( sequences = c_infols_n8[[6]], 
             names     = paste0( c_infols_n8[[1]], "_",
                                 c_infols_n8[[5]], "_|",
                                 c_infols_n8[[3]], "|_",
                                 c_infols_n8[[2]], "_",
                                 c_infols_n8[[4]]
             ), file.out  = paste0( "raw_data.n8/processed/pN8_", length( c_infols_n8[[1]] ) ,".fasta") )


# seq manipulation ------------

trimtool( propblank = 0.9, filedir = "./raw_data.n8/processed/pN8_3471_align.fasta" )


# remove outliers ------------

# 4 outliers found by TempEst
# MH597205_ruddy_turnstone_Delaware_1016417_2003_|USA|_H9N8_2003.384
# MH597724_ruddy_turnstone_New_Jersey_1321395_2005_|USA|_H3N8_2005.386
# MH597774_ruddy_turnstone_New_Jersey_1321408_2005_|USA|_H3N8_2005.397
# MH597785_ruddy_turnstone_New_Jersey_1321408B_2005_|USA|_H3N8_2005.397

n3471.seqname = fastaEx( "./raw_data.n8/processed/pN8_3471_trim2.fasta" )$id
n3471.seq     = fastaEx( "./raw_data.n8/processed/pN8_3471_trim2.fasta" )$seq

rm <- grep( "MH597205|MH597724|MH597774|MH597785", n3471.seqname )
write.fasta( n3471.seq[-rm], n3471.seqname[-rm], file.out = "./raw_data.n8/processed/pN8_3467.fasta" )






