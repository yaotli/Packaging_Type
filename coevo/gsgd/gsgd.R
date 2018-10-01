source( "./function.coevo.R" )

# extract gsgd-like and rebuild a tree 

ls.gsgd <- tagExtra( "./raw_data/processed/tree/fasttree_pH5_8334_e0921.tre" )
ls.gsgd <- ls.gsgd$id[ !is.na( ls.gsgd$tag ) ]
leafEx( "./raw_data/processed/pH5_8334_trim2.2.fasta", ls.gsgd )

# extract c2344 and remove duplicate 

ls.c2344 <- tagExtra( "./gsgd/processed/tree/fasttree_pH5_gsgd_7252_e0921.tre" )
ls.c2344 <- ls.c2344$id[ !is.na( ls.c2344$tag ) ]
leafEx( "./raw_data/processed/pH5_8334_trim2.2.fasta", ls.c2344 )

rmDup( "./gsgd/processed/pH5_8334_trim2.2_2659.fasta", rmdup = TRUE )
rmdup_plus( "./gsgd/processed/pH5_8334_trim2.2_2659_rmd.fasta" )