source( "./function.coevo.R" )
#source( "./ha_classi.R" )
source( "./f.aa_recon.R")
require(ggtree)
require(ape)
require(tidyverse)
require(stringr)


H5_treefile = "./gsgd/processed/tree/raxml_c2344_2657/raxml_pH5_2657.tre"
N6_trefile  = "./raw_data.n6/processed/tree/fasttree_pN6_4286_.tre"

H5_seq = "./gsgd/processed/pH5_c2344_2657.fasta"
N6_seq = "./raw_data.n6/processed/pN6_4286_trim2.fasta"


# pairing  H5-N6 -----------------------------------------------------

n6_tre   <- read.nexus( N6_trefile )
n6_root  <- length( n6_tre$tip.label ) + 1
n6_dismx <- dist.nodes( n6_tre )
n6_table <- fortify( n6_tre )

n6_i     <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.4$id ), collapse = "|" ), gsub( "\\|", "",  n6_table$label ) )
# length( ha_mdr.4$id) ==  length(n6_i) 

# fix discrepancy in id 
n6_id  <- gsub( "'", "", n6_table$label[ n6_i ] )
ha_na  <- match(  gsub(  "^[0-9A-Za-z]+|\\|", "", n6_id ), gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.4$id ) )

dm = c()
for( i in 1:length(ha_mdr.4$id ) ){ if( ! i %in% ha_na ){ dm = c(dm, i) }  }
dm.r = grep( paste0( gsub( "[0-9.]+$", "", ha_mdr.4$id[ dm ] ), collapse = "|" ), gsub( "\\|", "",  n6_table$label ) )
n6_table$label[ dm.r ][1] <- sub(  "[0-9.]+\\'$", str_match( ha_mdr.4$id[ dm ], "[0-9.]+$" )[1], n6_table$label[ dm.r ][1]  )
n6_table$label[ dm.r ][2] <- sub(  "[0-9.]+\\'$", str_match( ha_mdr.4$id[ dm ], "[0-9.]+$" )[2], n6_table$label[ dm.r ][2]  )


n6_i   <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.4$id ), collapse = "|" ), gsub( "\\|", "",  n6_table$label ) )
# length( ha_mdr.4$id) ==  length(n6_i) 

n6_id  <- gsub( "'", "", n6_table$label[ n6_i ] )
n6_mdr <- treeMDR( n6_i, n6_dismx )

ha_na <- match(  gsub(  "^[0-9A-Za-z]+|\\|", "", n6_id ), gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.4$id ) )

n6_mdr$ix    = n6_i
n6_mdr$group = ha_mdr.4$group[ ha_na ]
n6_mdr$id    = gsub(  "^[0-9A-Za-z]+|\\|", "", n6_id )
n6_mdr$type  = "N"
n6_mdr$sero  = "H5N6"

ha_mdr.4$type = "H"


# # V1
# ggplot( n6_mdr, aes( x = Dim_1, y = Dim_2, label = id ) )  + geom_point( aes(color = group), alpha = 0.5, size = 5) 
# 
# # V2
# ggplot( rbind( ha_mdr.4, n6_mdr ), aes( x = Dim_1, y = Dim_2, label = id ) )  + 
#   geom_point( aes(color = group, alpha = type ), size = 5) +
#   geom_line( aes(group = id), size = 0.1) + 
#   geom_rect( aes( xmin = h5n6_g1[1], xmax = h5n6_g1[2], ymin = h5n6_g1[3], ymax = h5n6_g1[4] ), inherit.aes = FALSE, color = "red", fill = NA) + 
#   geom_rect( aes( xmin = h5n6_g2[1], xmax = h5n6_g2[2], ymin = h5n6_g2[3], ymax = h5n6_g2[4] ), inherit.aes = FALSE, color = "red", fill = NA) +
#   geom_rect( aes( xmin = h5n6_g3[1], xmax = h5n6_g3[2], ymin = h5n6_g3[3], ymax = h5n6_g3[4] ), inherit.aes = FALSE, color = "red", fill = NA)  
# # coord_cartesian( xlim = c(0, 0.05), ylim = c(0, 0.02) ) +
# # geom_text( aes(alpha = type), size = 2, vjust = 1) + 
# # scale_y_continuous( limits = c( -0.1,  0.2) ) 
# 
# # V3
# N6_trein = treeio::read.nexus( N6_trefile )
# N6_tredf = fortify( N6_trein )
# N6_tredf$shape = NA
# N6_tredf$group = NA
# 
# N6_tredf$shape[ n6_mdr$ix ] = 1
# ggtree( N6_trein, right = TRUE ) %<+% N6_tredf + geom_tippoint( aes( shape = I(shape) ), color = "red", size = 5, alpha = 0.5 )
# 
# # V4
# N6_tredf$group[ n6_mdr$ix ] = n6_mdr$group
# ggtree( N6_trein, right = TRUE ) %<+% N6_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5, alpha = 0.5 )
# 
# # V5
# N6_tredf$shape[ g3_out$ix ] = 19
# ggtree( N6_trein, right = TRUE ) %<+% N6_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5 )


# grouping   
h5n6_g1 <- c( -0.2, -0.175, -0.125, -0.1 )
h5n6_g2 <- c( -0.075, -0.05, 0.075, 0.12 )
h5n6_g3 <- c( -0.025, 0.05, -0.025, 0.025 )


# extract 1
g1_out = 
  n6_mdr %>% 
  filter( Dim_1 > h5n6_g1[1] & Dim_1 < h5n6_g1[2] ) %>%
  filter( Dim_2 > h5n6_g1[3] & Dim_2 < h5n6_g1[4] ) %>%
  filter( group == 2  ) %>%
  select( ix )

# extract 2
g2_out = 
  n6_mdr %>% 
  filter( Dim_1 > h5n6_g2[1] & Dim_1 < h5n6_g2[2] ) %>%
  filter( Dim_2 > h5n6_g2[3] & Dim_2 < h5n6_g2[4] ) %>%
  filter( group == 1  ) %>%
  select( ix )

# extract 3
g3_out = 
  n6_mdr %>% 
  filter( Dim_1 > h5n6_g3[1] & Dim_1 < h5n6_g3[2] ) %>%
  filter( Dim_2 > h5n6_g3[3] & Dim_2 < h5n6_g3[4] ) %>%
  filter( group == 3  ) %>%
  select( ix )


# output 

#
g1_na <- gsub( "'", "", n6_table$label )[ g1_out$ix ]
g1_ha <- gsub( "'", "", ha_table$label[ ha_mdr.4$ix[ ha_na[ match( g1_out$ix, n6_mdr$ix ) ] ] ] )
leafEx( H5_seq, g1_ha, seq.out = "./h5_n6/pHA_h5n6_g1.fasta")
m.k = match( str_match( g1_na, "^[A-Z0-9a-z]+" )[,1], str_match( fastaEx( N6_seq )$id, "^[A-Z0-9a-z]+" )[,1] )
write.fasta( sequences = fastaEx( N6_seq )$seq[ m.k ], names = fastaEx( N6_seq )$id[ m.k ], 
             file.out  = "./h5_n6/pNA_h5n6_g1.fasta" ) # manually fix the year 

#
g2_na <- gsub( "'", "", n6_table$label )[ g2_out$ix ]
g2_ha <- gsub( "'", "", ha_table$label[ ha_mdr.4$ix[ ha_na[ match( g2_out$ix, n6_mdr$ix ) ] ] ] )
leafEx( H5_seq, g2_ha, seq.out = "./h5_n6/pHA_h5n6_g2.fasta")
leafEx( N6_seq, g2_na, seq.out = "./h5_n6/pNA_h5n6_g2.fasta")  

#
g3_na <- gsub( "'", "", n6_table$label )[ g3_out$ix ]
g3_ha <- gsub( "'", "", ha_table$label[ ha_mdr.4$ix[ ha_na[ match( g3_out$ix, n6_mdr$ix ) ] ] ] )
leafEx( H5_seq, g3_ha, seq.out = "./h5_n6/pHA_h5n6_g3.fasta")
m.l = match( str_match( g3_na, "^[A-Z0-9a-z]+" )[,1], str_match( fastaEx( N6_seq )$id, "^[A-Z0-9a-z]+" )[,1] )
write.fasta( sequences = fastaEx( N6_seq )$seq[ m.l ], names = fastaEx( N6_seq )$id[ m.l ], 
             file.out  = "./h5_n6/pNA_h5n6_g3.fasta" ) # manually fix the year 


# aa reconstruction  -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n6/dS/h5n6_g1_h5/" )
.aa_recon( folderdir = "./h5_n6/dS/h5n6_g1_n6/" )

# sam2
.aa_recon( folderdir = "./h5_n6/dS/h5n6_g2_h5/" )
.aa_recon( folderdir = "./h5_n6/dS/h5n6_g2_n6/" )

# sam3
.aa_recon( folderdir = "./h5_n6/dS/h5n6_g3_h5/" )
.aa_recon( folderdir = "./h5_n6/dS/h5n6_g3_n6/" )



# 
# 2nd samples  -----------------------------------------------------

# g1
.root_seq(  seqfile = "./h5_n6/pHA_h5n6_g1.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n6/pNA_h5n6_g1.fasta", treefile = N6_trefile, originfas = N6_seq )

# g2
.root_seq(  seqfile = "./h5_n6/pHA_h5n6_g2.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n6/pNA_h5n6_g2.fasta", treefile = N6_trefile, originfas = N6_seq )

# g3
.root_seq(  seqfile = "./h5_n6/pHA_h5n6_g3.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n6/pNA_h5n6_g3.fasta", treefile = N6_trefile, originfas = N6_seq )


# aa reconstruction  -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n6/dS_r/h5n6_g1_h5/" )
.aa_recon( folderdir = "./h5_n6/dS_r/h5n6_g1_n6/" )

# sam2
.aa_recon( folderdir = "./h5_n6/dS_r/h5n6_g2_h5/" )
.aa_recon( folderdir = "./h5_n6/dS_r/h5n6_g2_n6/" )

# sam3
.aa_recon( folderdir = "./h5_n6/dS_r/h5n6_g3_h5/" )
.aa_recon( folderdir = "./h5_n6/dS_r/h5n6_g3_n6/" )


