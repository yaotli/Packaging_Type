source( "./function.coevo.R" )


# database download ------------
#
# GISAID -> 1107 isolates (1108 in the file)
# 1. type 'A'; 2. N '2' 
# 3. host undo 'Lab derived, Unknown, NA, Human*'
# 4. location 'all'
# 5. required seqments 'HA, NA'
# 6. mim length '1000'; 7. only GISAID 
# format 'Isolate name Type Collection date Isolate ID'
# 
# NCBI -> 31535 isolates
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


fas_n2_g <- "./raw_data.n2/sources/pN2_G_1108_20181205.fasta"
csv_g    <- "./raw_data.n2/sources/G_20181205.csv"
fas_n2_n <- "./raw_data.n2/sources/pN2_N_31535_20181205.fasta"

fas_curated_h5 <- "./raw_data/processed/pH5_8334.fasta"


# data cleaning  ------------

# GISAID 
# to remove one isolate which appears twice 
#

seq_n2_g    <- keepLongSeq( fastaEx( fas_n2_g )$seq, fastaEx( fas_n2_g )$id )$seq 
id_n2_g     <- keepLongSeq( fastaEx( fas_n2_g )$seq, fastaEx( fas_n2_g )$id )$id
infols_n2_g <- idInfo.na( rawid = id_n2_g, datasource = "g", g.csv = csv_g, na_subtype = 2 )

# NCBI
# 

seq_n2_n    <- fastaEx( fas_n2_n )$seq
id_n2_n     <- fastaEx( fas_n2_n )$id
infols_n2_n <- idInfo.na( rawid = id_n2_n, datasource = "n", na_subtype = 2 )

# remove no 'year' information

rm.i = which( is.na(infols_n2_n[[4]]) )

seq_n2_n    <- seq_n2_n[ -rm.i ]
id_n2_n     <- id_n2_n[ -rm.i ]
infols_n2_n <- idInfo.na( rawid = id_n2_n, datasource = "n", na_subtype = 2 )


# combine and remove the duplicate isolate ------------

seq_n2    <- c( seq_n2_g, seq_n2_n )
infols_n2 <- lapply( as.list( seq( 1, length( infols_n2_n ) ) ), function(x){ c( infols_n2_g[[x]], infols_n2_n[[x]])  }  )

infols_n2[[ length(infols_n2) + 1 ]] = seq_n2

s_infols_n2      <- strainSelect( infols_n2 )
s_infols_n2[[4]] <- seqDate( s_infols_n2[[4]] )

s_infols_n2[[4]][ which( s_infols_n2[[4]] < 1900 ) ] = "1900.00"

c_infols_n2 <- seqSelect( minlth = 1200, maxamb = 1, s_infols_n2, rmdup = FALSE ) # n = 29606


# examine HA-NA information ------------

info_h5 <- taxaInfo( fas_curated_h5 )

h5.i <- grep( "H5", c_infols_n2[[2]] )
m.i  <- match( c_infols_n2[[5]][h5.i], info_h5[[5]] )


dismatch = c()
for( i in 1: length( m.i ) ) 
{
  if( !is.na(m.i[i]) )
  {
    if( c_infols_n2[[3]][h5.i][i] !=  info_h5[[2]][m.i][i] ){ dismatch = c(dismatch, i) }
    if( abs( as.numeric( c_infols_n2[[4]][h5.i][i] ) - info_h5[[4]][m.i][i] ) > 0.1 ){ dismatch = c(dismatch, i) }
    
  }
}

# c_infols_n2[[4]][h5.i][ dismatch ]
# info_h5[[4]][ m.i ][ dismatch ]
# temporally ignore at this instance

    
# export ------------

write.fasta( sequences = c_infols_n2[[6]], 
             names     = paste0( c_infols_n2[[1]], "_",
                                 c_infols_n2[[5]], "_|",
                                 c_infols_n2[[3]], "|_",
                                 c_infols_n2[[2]], "_",
                                 c_infols_n2[[4]]
             ), file.out  = paste0( "raw_data.n2/processed/pN2_", length( c_infols_n2[[1]] ) ,".fasta") )


# seq manipulation ------------

trimtool( propblank = 0.9, filedir = "./raw_data.n2/processed/pN2_29606_align.fasta" )


# exclude seasonal flu ------------

n2_treetag <- tagExtra( "./raw_data.n2/processed/tree/fasttree_pN2_29606.tre" )

id_list = c( n2_treetag$id[ is.na(n2_treetag$tag) ], 
             n2_treetag$id[ !is.na(n2_treetag$tag) ][ which.min( str_match( n2_treetag$id[ !is.na(n2_treetag$tag) ], "[0-9.]+$" ) ) ] )


leafEx( "./raw_data.n2/processed/pN2_29606_trim2.fasta", leaflist = id_list )


# remove outlier ------------

# I
# outliers found by TempEst
# GU247937_chicken_Korea_580_2006_|South_Korea|_H9N2_2006.318
# MF919794_swine_Minnesota_MT_13_01_S79_2013_|USA|_H1N2_2013.277
# GU247942_chicken_Korea_97_2007_|South_Korea|_H9N2_2007.041
# GU247936_chicken_Korea_196_2006_|South_Korea|_H9N2_2006.129

n7166.seqname = fastaEx( "./raw_data.n2/processed/pN2_7166.fasta" )$id
n7166.seq     = fastaEx( "./raw_data.n2/processed/pN2_7166.fasta" )$seq
rm <- grep( "GU247937|MF919794|GU247942|GU247936", n7166.seqname )
write.fasta( n7166.seq[-rm], n7166.seqname[-rm], file.out = "./raw_data.n2/processed/pN2_7162.fasta" )


# II
# EPI208666_Chicken_Cambodia_Pror_Yuk_030501_2011_H6N2_|Cambodia|_H6N2_1900.00
# GU122057_Pekin_duck_Korea_LBM186_2007_|South_Korea|_H9N2_2007.123
# MH271630_chicken_China_233_2017_|China|_H9N2_2017.315
# CY116845_mallard_MT_Y61_|Russia|_H2N2_1900.000

n7162.seqname = fastaEx( "./raw_data.n2/processed/pN2_7162.fasta" )$id
n7162.seq     = fastaEx( "./raw_data.n2/processed/pN2_7162.fasta" )$seq
rm <- grep( "EPI208666|GU122057|MH271630|CY116845", n7162.seqname )
write.fasta( n7162.seq[-rm], n7162.seqname[-rm], file.out = "./raw_data.n2/processed/pN2_7158.fasta" ) # resulting = 7158






# # locate gsgd h5n2 viruses ------------
# 
# h5_raw_tre           <- treeio::read.nexus( "./gsgd/processed/tree/fasttree_pH5_gsgd_7252_e0921.tre" )
# h5_raw_tre$tip.label <- gsub( "^[A-Za-z0-9']+", "", h5_raw_tre$tip.label )
# treeio::write.tree( h5_raw_tre, "./gsgd/processed/tree/fasttree_pH5_gsgd_7252_.tre" )
# 
# n2_raw_tre           <- treeio::read.nexus( "./raw_data.n2/processed/tree/fasttree_pN2_7166_.tre" )
# n2_raw_tre$tip.label <- gsub( "^[A-Za-z0-9']+", "", n2_raw_tre$tip.label )
# treeio::write.tree( n2_raw_tre, "./raw_data.n2/processed/tree/fasttree_pN2_7166_.nwk" )
# 
# 
# # scale down the tree ------------
# 
# rmDup( "./raw_data.n2/processed/pN2_7166.fasta", rmdup = TRUE )
# rmdup_plus( "./raw_data.n2/processed/pN2_7166_rmd.fasta" ) #n = 5583
# 
# # remove sequences which lack time information 
# 
# tem_id <- fastaEx( "./raw_data.n2/processed/pN2_7166_rmd_rmd2.fasta" )$id
# rm.j   <- which( as.numeric( str_match( tem_id, "[0-9.]+$" ) ) == 1900 ) 
# 
# write.fasta( sequences = fastaEx( "./raw_data.n2/processed/pN2_7166_rmd_rmd2.fasta" )$seq[-rm.j],
#              names     = fastaEx( "./raw_data.n2/processed/pN2_7166_rmd_rmd2.fasta" )$id[-rm.j], 
#              file.out  = "./raw_data.n2/processed/pN2_5581.fasta")
# 
# # I
# # find  outliers in TempEst 
# # MF919794_swine_Minnesota_MT_13_01_S79_2013_|USA|_H1N2_2013.277
# # GU247937_chicken_Korea_580_2006_|South_Korea|_H9N2_2006.318
# # GU247942_chicken_Korea_97_2007_|South_Korea|_H9N2_2007.041
# # GU247936_chicken_Korea_196_2006_|South_Korea|_H9N2_2006.129
# 
# n5581.seqname = fastaEx( "./raw_data.n2/processed/pN2_5581.fasta" )$id
# n5581.seq     = fastaEx( "./raw_data.n2/processed/pN2_5581.fasta" )$seq
# 
# rm.k <- grep( "MF919794|GU247937|GU247942|GU247936", n5581.seqname )
# write.fasta( n5581.seq[-rm.k], n5581.seqname[-rm.k], file.out = "./raw_data.n2/processed/pN2_5577.fasta" )
# 
# 
# # II
# # HQ143715_silky_fowl_Korea_LBM446_2007_|South_Korea|_H9N2_2007.200
# # MH271630_chicken_China_233_2017_|China|_H9N2_2017.315
# 
# n5577.seqname = fastaEx( "./raw_data.n2/processed/pN2_5577.fasta" )$id
# n5577.seq     = fastaEx( "./raw_data.n2/processed/pN2_5577.fasta" )$seq
# 
# rm.l <- grep( "HQ143715|MH271630", n5577.seqname )
# write.fasta( n5577.seq[-rm.l], n5577.seqname[-rm.l], file.out = "./raw_data.n2/processed/pN2_5575.fasta" )
# 
# 
# # finalize the working alignment 
# 
# N2subtreeNorAm <- treeio::read.nexus( "./raw_data.n2/processed/tree/pN2_NorthAm.tre" )
# leafEx( "./raw_data.n2/processed/pN2_29606_trim2.fasta", leaflist = gsub( "'", "", N2subtreeNorAm$tip.label) )
# 
# # fix EPI237995_turkey_Minnesota_9845_4_2015_|USA|_H5N2_2005.000 as h5 alignemnt 
# 
# n542.seqname = fastaEx( "./raw_data.n2/processed/pN2_542.fasta" )$id
# n542.seq     = fastaEx( "./raw_data.n2/processed/pN2_542.fasta" )$seq
# 
# n542.seqname[ grep( "EPI237995", n542.seqname ) ] = gsub( "2005.000$", "2015.000", n542.seqname[ grep( "EPI237995", n542.seqname ) ] )
# write.fasta( n542.seq,  n542.seqname, file.out = "./raw_data.n2/processed/pN2_542.fasta")







