source("./function.coevo.R")
require( ape )
require( tidyverse )
require( ggtree )
require( ggpubr )
require( seqinr )


# extract sequence and annotated the aa change 

beast.c2344anc_139 <- "./anc_seq_reconstruct/beast/c2344anc_139/20181020_pH5_c2344anc_139-anno.tre"

# print node number
# tre.c2344anc_139   <- read.beast( beast.c2344anc_139 )
# ggtree( tre.c2344anc_139, mrsd = "2014-02-06", right = TRUE ) + geom_text2( aes( label = node ), size = 2 )

branchAA( trefile = beast.c2344anc_139, saveFasta = T, writeTre = T )

# examine the findings 
# pos1 - branch h5n1-h5nx
pos1 = c( 98, 139, 172, 199, 234,239, 256, 279 )

for( p in 1: length(pos1) )
{
  t = 
  viewResi( trefile = "./gsgd/processed/tree/raxml_gsgd_5180/raxml_pH5_gsgd_5180_e1026.tre",
            fasfile = "./gsgd/processed/pH5_8334_trim2.2_7245_rmd_rmd2.fasta",
            pos     =  pos1[p]  )
  ggsave( plot = t, filename = paste0( "resi_", pos1[p], ".pdf"), device = "pdf", width = 4, height = 4, units = "in")
}

# pos2 - branch between first and the tmrca of 4 clades
pos2 = unique( c( 88, 69, 171, 178, 208, 214, 289, 541, 61, 111, 131, 143, 145, 149, 185, 205 ) )
for( p in 1: length(pos2) )
{
  t = 
    viewResi( trefile = "./gsgd/processed/tree/raxml_gsgd_5180/raxml_pH5_gsgd_5180_e1026.tre",
              fasfile = "./gsgd/processed/pH5_8334_trim2.2_7245_rmd_rmd2.fasta",
              pos     =  pos2[p]  )
  ggsave( plot = t, filename = paste0( "resi_", pos2[p], ".pdf"), device = "pdf", width = 4, height = 4, units = "in")
}

# pos3 - sites found by MEME and 4 high freq changing sites found by branchAA 
pos3 = unique( c( 2, 3, 14, 21, 27, 31, 52, 102, 154, 211, 405, 455, 10, 14, 139, 140 ) )
for( p in 1: length(pos3) )
{
  t = 
    viewResi( trefile = "./gsgd/processed/tree/raxml_gsgd_5180/raxml_pH5_gsgd_5180_e1026.tre",
              fasfile = "./gsgd/processed/pH5_8334_trim2.2_7245_rmd_rmd2.fasta",
              pos     =  pos3[p]  )
  ggsave( plot = t, filename = paste0( "resi_", pos3[p], ".pdf"), device = "pdf", width = 4, height = 4, units = "in")
}


# pos4 - sites found by MEME in the 4 clades (most are clade D)
pos4 = unique( c( 552, 14, 10, 131, 142, 154, 239, 291, 292) )
for( p in 1: length(pos4) )
{
  t = 
    viewResi( trefile = "./gsgd/processed/tree/raxml_gsgd_5180/raxml_pH5_gsgd_5180_e1026.tre",
              fasfile = "./gsgd/processed/pH5_8334_trim2.2_7245_rmd_rmd2.fasta",
              pos     =  pos4[p]  )
  ggsave( plot = t, filename = paste0( "resi_", pos4[p], ".pdf"), device = "pdf", width = 4, height = 4, units = "in")
}




pos = 291
viewResi( trefile = "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1110_.tre",
          fasfile = "./gsgd/processed/pH5_c2344_1598.fasta",
          pos     =  pos  )
ha_num( ref_fas = "./aa_numbering/align_aa_ref.fasta", ref_csv = "./aa_numbering/ref_numbering.csv", data_pos = 140)