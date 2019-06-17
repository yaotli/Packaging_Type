source( "./function.coevo.R" )
source( "./ha_classi.R" )
source( "./f.aa_recon.R")
require(ggtree)
require(ape)
require(tidyverse)
require(stringr)


H5_treefile = "./gsgd/processed/tree/raxml_c2344_2657/raxml_pH5_2657.tre"
N5_trefile  = "./raw_data.n5/processed/tree/fasttree_pN5_565_.tre"

H5_seq = "./gsgd/processed/pH5_c2344_2657.fasta"
N5_seq = "./raw_data.n5/processed/pN5_565.fasta"



# pairing  H5-N5 -----------------------------------------------------

n5_tre   <- read.nexus( N5_trefile )
n5_root  <- length( n5_tre$tip.label ) + 1
n5_dismx <- dist.nodes( n5_tre )
n5_table <- fortify( n5_tre )

n5_i     <- grep( paste0( gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.3$id ), collapse = "|" ), gsub( "\\|", "",  n5_table$label ) )
# length( n5_i ) == length( ha_mdr.3$id )

n5_id  <- gsub( "'", "", n5_table$label[ n5_i ] )
n5_mdr <- treeMDR( n5_i, n5_dismx )

ha_na <- match(  gsub(  "^[0-9A-Za-z]+|\\|", "", n5_id ), gsub(  "^[0-9A-Za-z]+|\\|", "", ha_mdr.3$id ) )

n5_mdr$ix    = n5_i
n5_mdr$group = ha_mdr.3$group[ ha_na ]
n5_mdr$id    = gsub(  "^[0-9A-Za-z]+|\\|", "", n5_id )
n5_mdr$type  = "N"
n5_mdr$sero  = "H5N5"

ha_mdr.3$type = "H"



# V1
ggplot( n5_mdr, aes( x = Dim_1, y = Dim_2, label = id ) )  + geom_point( aes(color = group), alpha = 0.5, size = 5) 

# V2
ggplot( rbind( ha_mdr.3, n5_mdr ), aes( x = Dim_1, y = Dim_2, label = id ) )  + 
  geom_point( aes(color = group, alpha = type ), size = 5) +
  geom_line( aes(group = id), size = 0.1) + 
  geom_rect( aes( xmin = h5n5_g1[1], xmax = h5n5_g1[2], ymin = h5n5_g1[3], ymax = h5n5_g1[4] ), inherit.aes = FALSE, color = "red", fill = NA) + 
  geom_rect( aes( xmin = h5n5_g2[1], xmax = h5n5_g2[2], ymin = h5n5_g2[3], ymax = h5n5_g2[4] ), inherit.aes = FALSE, color = "red", fill = NA) 

# coord_cartesian( xlim = c(0.15, 0.25), ylim = c(0.025, 0.1) ) +
# geom_text( aes(alpha = type), size = 2, vjust = 1) + 
# scale_y_continuous( limits = c( 0.1,  0.2) ) 

# V3
N5_trein = treeio::read.nexus( N5_trefile )
N5_tredf = fortify( N5_trein )
N5_tredf$shape = NA
N5_tredf$group = NA

N5_tredf$shape[ n5_mdr$ix ] = 1
ggtree( N5_trein, right = TRUE ) %<+% N5_tredf + geom_tippoint( aes( shape = I(shape) ), color = "red", size = 5, alpha = 0.5 )

# V4
N5_tredf$group[ n5_mdr$ix ] = n5_mdr$group
ggtree( N5_trein, right = TRUE ) %<+% N5_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5, alpha = 0.5 )

# V5
N5_tredf$shape[ g2_out$ix ] = 19
ggtree( N5_trein, right = TRUE ) %<+% N5_tredf + geom_tippoint( aes( shape = I(shape), color = group ), size = 5 )


# grouping   
h5n5_g1 <- c( -0.075, -0.03, -0.04, 0 )
h5n5_g2 <- c( 0.025, 0.035, -0.02, 0.01 )


# extract 1
g1_out = 
  n5_mdr %>% 
  filter( Dim_1 > h5n5_g1[1] & Dim_1 < h5n5_g1[2] ) %>%
  filter( Dim_2 > h5n5_g1[3] & Dim_2 < h5n5_g1[4] ) %>%
  filter( group == 1  ) %>%
  select( ix )

# extract 2
g2_out = 
  n5_mdr %>% 
  filter( Dim_1 > h5n5_g2[1] & Dim_1 < h5n5_g2[2] ) %>%
  filter( Dim_2 > h5n5_g2[3] & Dim_2 < h5n5_g2[4] ) %>%
  filter( group == 2  ) %>%
  select( ix )


# output 

#
g1_na <- gsub( "'", "", n5_table$label )[ g1_out$ix ]
g1_ha <- gsub( "'", "", ha_table$label[ ha_mdr.3$ix[ ha_na[ match( g1_out$ix, n5_mdr$ix ) ] ] ] )
leafEx( H5_seq, g1_ha, seq.out = "./h5_n5/pHA_h5n5_g1.fasta")
leafEx( N5_seq, g1_na, seq.out = "./h5_n5/pNA_h5n5_g1.fasta" )  

#
g2_na <- gsub( "'", "", n5_table$label )[ g2_out$ix ]
g2_ha <- gsub( "'", "", ha_table$label[ ha_mdr.3$ix[ ha_na[ match( g2_out$ix, n5_mdr$ix ) ] ] ] )
leafEx( H5_seq, g2_ha, seq.out = "./h5_n5/pHA_h5n5_g2.fasta")
leafEx( N5_seq, g2_na, seq.out = "./h5_n5/pNA_h5n5_g2.fasta" )  





# aa reconstruction  -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n5/dS/h5n5_g1_h5/" )
.aa_recon( folderdir = "./h5_n5/dS/h5n5_g1_n5/" )

# sam2
.aa_recon( folderdir = "./h5_n5/dS/h5n5_g2_h5/" )
.aa_recon( folderdir = "./h5_n5/dS/h5n5_g2_n5/" )


# 
# 2nd samples  -----------------------------------------------------

# g1
.root_seq(  seqfile = "./h5_n5/pHA_h5n5_g1.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n5/pNA_h5n5_g1.fasta", treefile = N5_trefile, originfas = N5_seq )

# g2
.root_seq(  seqfile = "./h5_n5/pHA_h5n5_g2.fasta", treefile = H5_treefile, originfas = H5_seq )
.root_seq(  seqfile = "./h5_n5/pNA_h5n5_g2.fasta", treefile = N5_trefile, originfas = N5_seq )


# aa reconstruction  -----------------------------------------------------

# sam1
.aa_recon( folderdir = "./h5_n5/dS_r/h5n5_g1_h5/" )
.aa_recon( folderdir = "./h5_n5/dS_r/h5n5_g1_n5/" )

# sam2
.aa_recon( folderdir = "./h5_n5/dS_r/h5n5_g2_h5/" )
.aa_recon( folderdir = "./h5_n5/dS_r/h5n5_g2_n5/" )


