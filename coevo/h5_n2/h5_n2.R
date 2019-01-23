source( "./function.coevo.R" )
source( "./ha_classi.R" )
source( "./f.aa_recon.R")
require(ggtree)
require(ape)
require(tidyverse)
require(stringr)


H5_treefile = "./gsgd/processed/tree/raxml_c2344_2657/raxml_pH5_2657.tre"
N2_trefile  = "./raw_data.n2/processed/tree/fasttree_pN2_7158_.tre"

H5_seq = "./gsgd/processed/pH5_c2344_2657.fasta"
N2_seq = "./raw_data.n2/processed/pN2_7158.fasta"


# pairing  H5-N2 -----------------------------------------------------

n2_tre   <- read.nexus( N2_trefile )
n2_root  <- length( n2_tre$tip.label ) + 1
n2_dismx <- dist.nodes( n2_tre )
n2_table <- fortify( n2_tre )

n2_i     <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.2$id ), collapse = "|" ), gsub( "\\|", "",  n2_table$label ) )
# length( ha_mdr.2$id) ==  length(n2_i) 

dismatch <- grep( "_turkey_Minnesota_9845_4_", n2_table$label )
n2_table$label[ dismatch ] =  gsub( "_2005.000", "_2015.000", n2_table$label[ dismatch ] )

n2_i   <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.2$id ), collapse = "|" ), gsub( "\\|", "",  n2_table$label ) )

n2_id  <- gsub( "'", "", n2_table$label[ n2_i ] )
n2_mdr <- treeMDR( n2_i, n2_dismx )

ha_na <- match(  gsub(  "^[0-9A-Za-z]+|\\|", "", n2_id ), gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.2$id ) )

n2_mdr$ix    = n2_i
n2_mdr$group = ha_mdr.2$group[ ha_na ]
n2_mdr$id    = gsub(  "^[0-9A-Za-z]+|\\|", "", n2_id )
n2_mdr$type  = "N"
n2_mdr$sero  = "H5N2"

ha_mdr.2$type = "H"



# # V1
# ggplot( n2_mdr, aes( x = Dim_1, y = Dim_2, label = id ) )  + geom_point( aes(color = group), alpha = 0.5, size = 5) 
# 
# # V2
# ggplot( rbind( ha_mdr.2, n2_mdr ), aes( x = Dim_1, y = Dim_2, label = id ) )  + 
#   geom_point( aes(color = group, alpha = type ), size = 5) +
#   geom_line( aes(group = id), size = 0.1) + 
#   geom_rect( aes( xmin = h5n2_g1[1], xmax = h5n2_g1[2], ymin = h5n2_g1[3], ymax = h5n2_g1[4] ), inherit.aes = FALSE, color = "red", fill = NA) + 
#   geom_rect( aes( xmin = h5n2_g2[1], xmax = h5n2_g2[2], ymin = h5n2_g2[3], ymax = h5n2_g2[4] ), inherit.aes = FALSE, color = "red", fill = NA) +
#   geom_rect( aes( xmin = h5n2_g3[1], xmax = h5n2_g3[2], ymin = h5n2_g3[3], ymax = h5n2_g3[4] ), inherit.aes = FALSE, color = "red", fill = NA) + 
#   geom_rect( aes( xmin = h5n2_g4[1], xmax = h5n2_g4[2], ymin = h5n2_g4[3], ymax = h5n2_g4[4] ), inherit.aes = FALSE, color = "red", fill = NA) 
# 
# # coord_cartesian( xlim = c(0.15, 0.25), ylim = c(0.025, 0.1) ) +
# # geom_text( aes(alpha = type), size = 2, vjust = 1) + 
# # scale_y_continuous( limits = c( 0.1,  0.2) ) 
# 
# # V3
# N2_trein = treeio::read.nexus( N2_trefile )
# N2_tredf = fortify( N2_trein )
# N2_tredf$shape = NA
# N2_tredf$group = NA
# 
# N2_tredf$shape[ n2_mdr$ix ] = 1
# ggtree( N2_trein, right = TRUE ) %<+% N2_tredf + geom_tippoint( aes( shape = I(shape) ), color = "red", size = 5, alpha = 0.5 )
# 
# # V4
# N2_tredf$group[ n2_mdr$ix ] = n2_mdr$group
# ggtree( N2_trein, right = TRUE ) %<+% N2_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5, alpha = 0.5 )
# 
# # V5
# N2_tredf$shape[ g4_out$ix ] = 19
# ggtree( N2_trein, right = TRUE ) %<+% N2_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5 )


# grouping   
h5n2_g1 <- c( 0.15, 0.2, -0.01, 0.025 )
h5n2_g2 <- c( -0.04, -0.01, -0.005, 0.01 )
h5n2_g3 <- c( 0.12, 0.13, -0.15, -0.11 )
h5n2_g4 <- c( 0.175, 0.25, 0.04, 0.09 )


# extract 1
g1_out = 
  n2_mdr %>% 
  filter( Dim_1 > h5n2_g1[1] & Dim_1 < h5n2_g1[2] ) %>%
  filter( Dim_2 > h5n2_g1[3] & Dim_2 < h5n2_g1[4] ) %>%
  filter( group == 1  ) %>%
  select( ix )

# extract 2
g2_out = 
  n2_mdr %>% 
  filter( Dim_1 > h5n2_g2[1] & Dim_1 < h5n2_g2[2] ) %>%
  filter( Dim_2 > h5n2_g2[3] & Dim_2 < h5n2_g2[4] ) %>%
  filter( group == 4  ) %>%
  select( ix )

# extract 3
g3_out = 
  n2_mdr %>% 
  filter( Dim_1 > h5n2_g3[1] & Dim_1 < h5n2_g3[2] ) %>%
  filter( Dim_2 > h5n2_g3[3] & Dim_2 < h5n2_g3[4] ) %>%
  filter( group == 4  ) %>%
  select( ix )

# extract 4
g4_out = 
  n2_mdr %>% 
  filter( Dim_1 > h5n2_g4[1] & Dim_1 < h5n2_g4[2] ) %>%
  filter( Dim_2 > h5n2_g4[3] & Dim_2 < h5n2_g4[4] ) %>%
  filter( group == 3  ) %>%
  select( ix )



# output 

#
g1_na <- gsub( "'", "", n2_table$label )[ g1_out$ix ]
g1_ha <- gsub( "'", "", ha_table$label[ ha_mdr.2$ix[ ha_na[ match( g1_out$ix, n2_mdr$ix ) ] ] ] )
leafEx( H5_seq, g1_ha, seq.out = "./h5_n2/pHA_h5n2_g1.fasta")
leafEx( N2_seq, g1_na, seq.out = "./h5_n2/pNA_h5n2_g1.fasta" )  

#
g2_na <- gsub( "'", "", n2_table$label )[ g2_out$ix ]
g2_ha <- gsub( "'", "", ha_table$label[ ha_mdr.2$ix[ ha_na[ match( g2_out$ix, n2_mdr$ix ) ] ] ] )
leafEx( H5_seq, g2_ha, seq.out = "./h5_n2/pHA_h5n2_g2.fasta")
m.j = match( str_match( g2_na, "^[A-Z0-9a-z]+" )[,1], str_match( fastaEx( N2_seq )$id, "^[A-Z0-9a-z]+" )[,1] )
write.fasta( sequences = fastaEx( N2_seq )$seq[ m.j ], names = fastaEx( N2_seq )$id[ m.j ], 
             file.out  = "./h5_n2/pNA_h5n2_g2.fasta" ) # manually fix the year 

#
g3_na <- gsub( "'", "", n2_table$label )[ g3_out$ix ]
g3_ha <- gsub( "'", "", ha_table$label[ ha_mdr.2$ix[ ha_na[ match( g3_out$ix, n2_mdr$ix ) ] ] ] )
leafEx( H5_seq, g3_ha, seq.out = "./h5_n2/pHA_h5n2_g3.fasta")
leafEx( N2_seq, g3_na, seq.out = "./h5_n2/pNA_h5n2_g3.fasta" )  

#
g4_na <- gsub( "'", "", n2_table$label )[ g4_out$ix ]
g4_ha <- gsub( "'", "", ha_table$label[ ha_mdr.2$ix[ ha_na[ match( g4_out$ix, n2_mdr$ix ) ] ] ] )
leafEx( H5_seq, g4_ha, seq.out = "./h5_n2/pHA_h5n2_g4.fasta")
leafEx( N2_seq, g4_na, seq.out = "./h5_n2/pNA_h5n2_g4.fasta" )  




# aa reconstruction  -----------------------------------------------------

# sam1-ha
.aa_recon( folderdir = "./h5_n1/sam1/HA/" )
# sam1-na
.aa_recon( folderdir = "./h5_n1/sam1/NA/" )


# sam2-ha
.aa_recon( folderdir = "./h5_n1/sam2/HA/" )
# sam2-na
.aa_recon( folderdir = "./h5_n1/sam2/NA/" )




# dist. test  -----------------------------------------------------

