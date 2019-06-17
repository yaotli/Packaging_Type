source( "./function.coevo.R" )
#source( "./ha_classi.R" )
source( "./f.aa_recon.R")
require(ggtree)
require(ape)
require(tidyverse)
require(stringr)


H5_treefile = "./gsgd/processed/tree/raxml_c2344_2657/raxml_pH5_2657.tre"
N1_trefile  = "./raw_data/processed/tree/fasttree_pN1_4696_.tre"

H5_seq = "./gsgd/processed/pH5_c2344_2657.fasta"
N1_seq = "./raw_data/processed/pN1_4696_trim2.fasta"


# pairing  H5-N1 -----------------------------------------------------

n1_tre   <- read.nexus( N1_trefile )
n1_root  <- length( n1_tre$tip.label ) + 1
n1_dismx <- dist.nodes( n1_tre )
n1_table <- fortify( n1_tre )

n1_i     <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.1$id ), collapse = "|" ), gsub( "\\|", "",  n1_table$label ) )
# length( ha_mdr.1$id) ==  length(n1_i) 

n1_id  <- gsub( "'", "", n1_table$label[ n1_i ] )
n1_mdr <- treeMDR( n1_i, n1_dismx )

ha_na <- match(  gsub(  "^[0-9A-Za-z]+|\\|", "", n1_id ), gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.1$id ) )

n1_mdr$ix    = n1_i
n1_mdr$group = ha_mdr.1$group[ ha_na ]
n1_mdr$id    = gsub(  "^[0-9A-Za-z]+|\\|", "", n1_id )
n1_mdr$type  = "N"
n1_mdr$sero  = "H5N1"

ha_mdr.1$type = "H"



# # V1
# ggplot( n1_mdr, aes( x = Dim_1, y = Dim_2, label = id ) )  + geom_point( aes(color = group), alpha = 0.5, size = 5)
# 
# # V2
# ggplot( rbind( ha_mdr.1, n1_mdr ), aes( x = Dim_1, y = Dim_2, label = id ) )  +
#   geom_point( aes(color = group, alpha = type ), size = 5) +
#   geom_line( aes(group = id), size = 0.1) +
#   geom_rect( aes( xmin = h5n1_g1[1], xmax = h5n1_g1[2], ymin = h5n1_g1[3], ymax = h5n1_g1[4] ), inherit.aes = FALSE, color = "red", fill = NA) +
#   geom_rect( aes( xmin = h5n1_g2[1], xmax = h5n1_g2[2], ymin = h5n1_g2[3], ymax = h5n1_g2[4] ), inherit.aes = FALSE, color = "red", fill = NA)
# # coord_cartesian( xlim = c(0, 0.05), ylim = c(0, 0.02) ) +
# # geom_text( aes(alpha = type), size = 2, vjust = 1) +
# # scale_y_continuous( limits = c( -0.1,  0.2) )
# 
# # V3
# N1_trein = treeio::read.nexus( N1_trefile )
# N1_tredf = fortify( N1_trein )
# N1_tredf$shape = NA
# N1_tredf$group = NA
# 
# N1_tredf$shape[ n1_mdr$ix ] = 1
# ggtree( N1_trein, right = TRUE ) %<+% N1_tredf + geom_tippoint( aes( shape = I(shape) ), color = "red", size = 5, alpha = 0.5 )
# 
# # V4
# N1_tredf$group[ n1_mdr$ix ] = n1_mdr$group
# ggtree( N1_trein, right = TRUE ) %<+% N1_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5, alpha = 0.5 )
# 
# # V5
# N1_tredf$shape[ g1_out$ix ] = 19
# ggtree( N1_trein, right = TRUE ) %<+% N1_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5 )


# grouping   
h5n1_g1 <- c( 0, 0.01, 0.0075, 0.0175 )
h5n1_g2 <- c( -0.075, -0.025, -0.03, -0.01 )

# extract 1
g1_out = 
  n1_mdr %>% 
  filter( Dim_1 > h5n1_g1[1] & Dim_1 < h5n1_g1[2] ) %>%
  filter( Dim_2 > h5n1_g1[3] & Dim_2 < h5n1_g1[4] ) %>%
  filter( group == 1  ) %>%
  select( ix )

# extract 2
g2_out = 
  n1_mdr %>% 
  filter( Dim_1 > h5n1_g2[1] & Dim_1 < h5n1_g2[2] ) %>%
  filter( Dim_2 > h5n1_g2[3] & Dim_2 < h5n1_g2[4] ) %>%
  filter( group == 2  ) %>%
  select( ix )


# output 

#
g1_na <- gsub( "'", "", n1_table$label )[ g1_out$ix ]
g1_ha <- gsub( "'", "", ha_table$label[ ha_mdr.1$ix[ ha_na[ match( g1_out$ix, n1_mdr$ix ) ] ] ] )
leafEx( H5_seq, g1_ha, seq.out = "./h5_n1/pHA_h5n1_g1.fasta")
leafEx( N1_seq, g1_na, seq.out = "./h5_n1/pNA_h5n1_g1.fasta" )  

#
g2_na <- gsub( "'", "", n1_table$label )[ g2_out$ix ]
g2_ha <- gsub( "'", "", ha_table$label[ ha_mdr.1$ix[ ha_na[ match( g2_out$ix, n1_mdr$ix ) ] ] ] )
leafEx( H5_seq, g2_ha, seq.out = "./h5_n1/pHA_h5n1_g2.fasta")
leafEx( N1_seq, g2_na, seq.out = "./h5_n1/pNA_h5n1_g2.fasta")  




# aa reconstruction  -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n1/dS/h5n1_g1_h5/" )
.aa_recon( folderdir = "./h5_n1/dS/h5n1_g1_n1/" )

# sam2
.aa_recon( folderdir = "./h5_n1/dS/h5n1_g2_h5/" )
.aa_recon( folderdir = "./h5_n1/dS/h5n1_g2_n1/" )



# 
# 2nd samples  -----------------------------------------------------

.root_seq(  seqfile = "./h5_n1/pHA_h5n1_g2.fasta", H5_treefile, H5_seq )
.root_seq(  seqfile = "./h5_n1/pNA_h5n1_g2.fasta", N1_trefile, N1_seq )


# aa reconstruction  -----------------------------------------------------

# sam2
.aa_recon( folderdir = "./h5_n1/dS_r//h5n1_g2_h5/" )
.aa_recon( folderdir = "./h5_n1/dS_r/h5n1_g2_n1/" )























