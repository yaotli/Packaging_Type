source("./function.coevo.R")
require( ape )
require( tidyverse )
require( ggtree )
require( ggpubr )
require( seqinr )


# extract sequence and annotated the aa change 

beast.c2344anc_139 <- "./anc_seq_reconstruct/beast/c2344anc_139/20181020_pH5_c2344anc_139-anno.tre"
tre.c2344anc_139   <- read.beast( beast.c2344anc_139 )

# print node number
# ggtree( tre.c2344anc_139, mrsd = "2014-02-06", right = TRUE ) + geom_text2( aes( label = node ), size = 2 )

branchAA( trefile = beast.c2344anc_139, saveFasta = T, writeTre = T )

