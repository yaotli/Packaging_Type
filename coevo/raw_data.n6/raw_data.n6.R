source( "./function.coevo.R" )


# database download ------------
#
# GISAID -> 1499 isolates (1499 in the file)
# 1. type 'A'; 2. N '6' 
# 3. host undo 'Lab derived, Unknown, NA'
# 4. location 'all'
# 5. required seqments 'HA, NA'
# 6. mim length '1000'; 7. only GISAID 
# format 'Isolate name Type Collection date Isolate ID'
# 
# NCBI -> 3039 isolates
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


fas_n6_g <- "./raw_data.n6/sources/pN6_G_1499_20190105.fasta"
csv_g    <- "./raw_data.n6/sources/G_20190105.csv"
fas_n6_n <- "./raw_data.n6/sources/pN6_N_3039_20190105.fasta"

fas_curated_h5 <- "./raw_data/processed/pH5_8334.fasta"


# data cleaning  ------------

# GISAID 
# to remove one isolate which appears twice 
#

seq_n6_g    <- fastaEx( fas_n6_g )$seq
id_n6_g     <- fastaEx( fas_n6_g )$id

id_n6_g[4]  = gsub( "'", "", id_n6_g[4] )
infols_n6_g <- idInfo.na( rawid = id_n6_g, datasource = "g", g.csv = csv_g, na_subtype = 6 )

# NCBI
# 

seq_n6_n    <- fastaEx( fas_n6_n )$seq
id_n6_n     <- fastaEx( fas_n6_n )$id
infols_n6_n <- idInfo.na( rawid = id_n6_n, datasource = "n", na_subtype = 6 )


# combine and remove the duplicate isolate ------------

seq_n6    <- c( seq_n6_g, seq_n6_n )
infols_n6 <- lapply( as.list( seq( 1, length( infols_n6_n ) ) ), function(x){ c( infols_n6_g[[x]], infols_n6_n[[x]])  }  )

infols_n6[[ length(infols_n6) + 1 ]] = seq_n6

s_infols_n6      <- strainSelect( infols_n6 )
s_infols_n6[[4]] <- seqDate( s_infols_n6[[4]] )

c_infols_n6 <- seqSelect( minlth = 1200, maxamb = 1, s_infols_n6, rmdup = FALSE ) # n = 4286


# examine HA-NA information ------------

info_h5 <- taxaInfo( fas_curated_h5 )

h5.i <- grep( "H5", c_infols_n6[[2]] )
m.i  <- match( c_infols_n6[[5]][h5.i], info_h5[[5]] )


dismatch = c()
for( i in 1: length( m.i ) )
{
  if( !is.na(m.i[i]) )
  {
    if( c_infols_n6[[3]][h5.i][i] !=  info_h5[[2]][m.i][i] ){ dismatch = c(dismatch, i) }
    if( abs( as.numeric( c_infols_n6[[4]][h5.i][i] ) - info_h5[[4]][m.i][i] ) > 0.1 ){ dismatch = c(dismatch, i) }
    
  }
}

c_infols_n6[[4]][h5.i][ dismatch ] = info_h5[[4]][ m.i[ dismatch ] ]

    
# export ------------

write.fasta( sequences = c_infols_n6[[6]], 
             names     = paste0( c_infols_n6[[1]], "_",
                                 c_infols_n6[[5]], "_|",
                                 c_infols_n6[[3]], "|_",
                                 c_infols_n6[[2]], "_",
                                 c_infols_n6[[4]]
             ), file.out  = paste0( "raw_data.n6/processed/pN6_", length( c_infols_n6[[1]] ) ,".fasta") )


# seq manipulation ------------

trimtool( propblank = 0.9, filedir = "./raw_data.n6/processed/pN6_4286_align.fasta" )


# # scale down the tree ------------
# 
# rmDup( "./raw_data.n6/processed/pN6_4286_trim2.fasta", rmdup = TRUE )
# rmdup_plus( "./raw_data.n6/processed/pN6_4286_rmd.fasta" ) #n = 2552
# 
# 
# # extract EA gene pool and rebuild a tree
# 
# n6_EA_pool <- tagExtra( "./raw_data.n6/processed/tree/fasttree_pN6_2552_r.tre" )
# 
# leafEx( filedir = "./raw_data.n6/processed/pN6_4286_trim2.fasta", n6_EA_pool$id[ !is.na( n6_EA_pool$tag ) ] )







