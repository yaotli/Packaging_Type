source( "./function.coevo.R" )
#source( "./ha_classi.R" )
source( "./f.aa_recon.R")
require(ggtree)
require(ape)
require(tidyverse)
require(stringr)


H5_treefile = "./gsgd/processed/tree/raxml_c2344_2657/raxml_pH5_2657.tre"
N8_trefile  = "./raw_data.n8/processed/tree/fasttree_pN8_3467_sub1237.tre"  

H5_seq = "./gsgd/processed/pH5_c2344_2657.fasta"
N8_seq = "./raw_data.n8/processed/pN8_3467.fasta"


# pairing  H5-N8 -----------------------------------------------------

n8_tre   <- read.nexus( N8_trefile )
n8_root  <- length( n8_tre$tip.label ) + 1
n8_dismx <- dist.nodes( n8_tre )
n8_table <- fortify( n8_tre )

n8_i     <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.5$id ), collapse = "|" ), gsub( "\\|", "",  n8_table$label ) )
# length( ha_mdr.5$id) ==  length(n8_i) 

n8_id  <- gsub( "'", "", n8_table$label[ n8_i ] )
n8_mdr <- treeMDR( n8_i, n8_dismx )

ha_na <- match(  gsub(  "^[0-9A-Za-z]+|\\|", "", n8_id ), gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.5$id ) )

n8_mdr$ix    = n8_i
n8_mdr$group = ha_mdr.5$group[ ha_na ]
n8_mdr$id    = gsub(  "^[0-9A-Za-z]+|\\|", "", n8_id )
n8_mdr$type  = "N"
n8_mdr$sero  = "H5N8"

ha_mdr.5$type = "H"


# # V1
# ggplot( n8_mdr, aes( x = Dim_1, y = Dim_2, label = id ) )  + geom_point( aes(color = group), alpha = 0.5, size = 5) 
# + geom_text( aes(alpha = type), size = 2, vjust = 1)  
# 
# # V2
# ggplot( rbind( ha_mdr.5, n8_mdr ), aes( x = Dim_1, y = Dim_2, label = id ) )  + 
#   geom_point( aes(color = group, alpha = type ), size = 5) +
#   geom_line( aes(group = id), size = 0.1) + 
#   geom_rect( aes( xmin = h5n8_g1[1], xmax = h5n8_g1[2], ymin = h5n8_g1[3], ymax = h5n8_g1[4] ), inherit.aes = FALSE, color = "red", fill = NA) + 
#   geom_rect( aes( xmin = h5n8_g2[1], xmax = h5n8_g2[2], ymin = h5n8_g2[3], ymax = h5n8_g2[4] ), inherit.aes = FALSE, color = "red", fill = NA) 
# # coord_cartesian( xlim = c(0, 0.05), ylim = c(0, 0.02) ) +
# # geom_text( aes(alpha = type), size = 2, vjust = 1) + 
# # scale_y_continuous( limits = c( -0.1,  0.2) ) 
# 
# # V3
# N8_trein = treeio::read.nexus( N8_trefile )
# N8_tredf = fortify( N8_trein )
# N8_tredf$shape = NA
# N8_tredf$group = NA
# 
# N8_tredf$shape[ n8_mdr$ix ] = 1
# ggtree( N8_trein, right = TRUE ) %<+% N8_tredf + geom_tippoint( aes( shape = I(shape) ), color = "red", size = 5, alpha = 0.5 )
# 
# # V4
# N8_tredf$group[ n8_mdr$ix ] = n8_mdr$group
# ggtree( N8_trein, right = TRUE ) %<+% N8_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5, alpha = 0.5 )
# 
# # V5
# N8_tredf$shape[ g2_out$ix ] = 19
# ggtree( N8_trein, right = TRUE ) %<+% N8_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5 )


# grouping   
h5n8_g1 <- c( -0.06, -0.025, -0.01, 0.02 )
h5n8_g2 <- c( -0.025, 0.04, -0.04, 0.01 )

# extract 1
g1_out = 
  n8_mdr %>% 
  filter( Dim_1 > h5n8_g1[1] & Dim_1 < h5n8_g1[2] ) %>%
  filter( Dim_2 > h5n8_g1[3] & Dim_2 < h5n8_g1[4] ) %>%
  filter( group == 1  ) %>%
  select( ix )

# extract 2
g2_out = 
  n8_mdr %>% 
  filter( Dim_1 > h5n8_g2[1] & Dim_1 < h5n8_g2[2] ) %>%
  filter( Dim_2 > h5n8_g2[3] & Dim_2 < h5n8_g2[4] ) %>%
  filter( group == 2  ) %>%
  select( ix )


# output 

#
g1_na <- gsub( "'", "", n8_table$label )[ g1_out$ix ]
g1_ha <- gsub( "'", "", ha_table$label[ ha_mdr.5$ix[ ha_na[ match( g1_out$ix, n8_mdr$ix ) ] ] ] )
leafEx( H5_seq, g1_ha, seq.out = "./h5_n8/pHA_h5n8_g1.fasta")
leafEx( N8_seq, g1_na, seq.out = "./h5_n8/pNA_h5n8_g1.fasta")  

#
g2_na <- gsub( "'", "", n8_table$label )[ g2_out$ix ]
g2_ha <- gsub( "'", "", ha_table$label[ ha_mdr.5$ix[ ha_na[ match( g2_out$ix, n8_mdr$ix ) ] ] ] )
leafEx( H5_seq, g2_ha, seq.out = "./h5_n8/pHA_h5n8_g2.fasta")
leafEx( N8_seq, g2_na, seq.out = "./h5_n8/pNA_h5n8_g2.fasta")  




# aa reconstruction I -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n8/dS/h5n8_g1_h5/" )
.aa_recon( folderdir = "./h5_n8/dS/h5n8_g1_n8/" )

# sam2
.aa_recon( folderdir = "./h5_n8/dS/h5n8_g2_h5/" )
.aa_recon( folderdir = "./h5_n8/dS/h5n8_g2_n8/" )


# 
# 2nd samples  -----------------------------------------------------

# g1
.root_seq(  seqfile = "./h5_n8/pHA_h5n8_g1.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n8/pNA_h5n8_g1.fasta", treefile = N8_trefile, originfas = N8_seq )

# g2
.root_seq(  seqfile = "./h5_n8/pHA_h5n8_g2.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n8/pNA_h5n8_g2.fasta", treefile = N8_trefile, originfas = N8_seq )


# aa reconstruction II -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n8/dS_r/h5n8_g1_h5/" )
.aa_recon( folderdir = "./h5_n8/dS_r/h5n8_g1_n8/" )

# sam2
.aa_recon( folderdir = "./h5_n8/dS_r/h5n8_g2_h5/" )
.aa_recon( folderdir = "./h5_n8/dS_r/h5n8_g2_n8/" )


# aa reconstruction III -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n8/dS_com/h5n8_g0_h5/" )
.aa_recon( folderdir = "./h5_n8/dS_com/h5n8_g0_n8/" )




