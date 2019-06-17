source( "./function.coevo.R" )
require(ggtree)
require(ape)
require(tidyverse)
require(stringr)


H5_treefile = "./gsgd/processed/tree/raxml_c2344_2657/raxml_pH5_2657.tre"
N1_trefile  = "./raw_data/processed/tree/fasttree_pN1_4696_.tre"
N2_trefile  = "./raw_data.n2/processed/tree/fasttree_pN2_7158_.tre"
N5_trefile  = "./raw_data.n5/processed/tree/fasttree_pN5_565_.tre"
N6_trefile  = "./raw_data.n6/processed/tree/fasttree_pN6_4286_.tre"
N8_trefile  = "./raw_data.n8/processed/tree/fasttree_pN8_3467_sub1237.tre"  

H5_seq = "./gsgd/processed/pH5_c2344_2657.fasta"
N1_seq = "./raw_data/processed/pN1_4696_trim2.fasta"
N2_seq = "./raw_data.n2/processed/pN2_7158.fasta"
N5_seq = "./raw_data.n5/processed/pN5_565.fasta"
N6_seq = "./raw_data.n6/processed/pN6_4286_trim2.fasta"
N8_seq = "./raw_data.n8/processed/pN8_3467.fasta"


# visualization  ----------------------------------------------------------

# H5_trein = treeio::read.nexus( H5_treefile )
# H5_tredf = fortify( H5_trein )
# H5_tredf$shape = NA
# H5_tredf$group = NA
# 
# # to map the targets on the tree
# H5_tredf$shape[ ha_mdr.2$ix ] = "20"
# ggtree( H5_trein, right = TRUE ) %<+% H5_tredf + geom_tippoint( aes( shape = shape ), color = "red", size = 5, alpha = 0.5 )
# 
# H5_tredf$group[ ha_mdr.2$ix ] = ha_mdr.2$group
# ggtree( H5_trein, right = TRUE ) %<+% H5_tredf + geom_tippoint( aes( shape = shape, color = group ), size = 5, alpha = 0.5 )
# 
# # for specific subtype with classification
# ggplot( ha_mdr.2, aes( x = Dim_1, y = Dim_2, label = id ) )  + geom_point( aes(color = group), alpha = 0.5, size = 4)
# 
# # for reduced dimension plot
# ggplot( rbind( ha_mdr.1, ha_mdr.2, ha_mdr.3, ha_mdr.4, ha_mdr.5 ), aes( x = Dim_1, y = Dim_2, label = id ) )  +
#   geom_point( aes(color = sero, alpha = sero ), size = 4) +
#   scale_alpha_manual( values = c( 0, 1, 0, 0, 0) ) +
#   # geom_text( aes( alpha = sero ), size = 2, vjust = 1)


# locate the Nx virus on HA tree  -----------------------------------------------------

# H5-Nx
h5_tre   <- read.nexus( H5_treefile )
h5_root  <- length( h5_tre$tip.label ) + 1
h5_dismx <- dist.nodes( h5_tre )
ha_table <- fortify( h5_tre )

sero = c( "H5N1", "H5N2", "H5N5", "H5N6", "H5N8" )
for( i in 1: 5 )
{

  ha_i   <- which( str_match( ha_table$label, "_(H5N[0-9]{1,2})_" )[,2] == sero[i] )
  ha_id  <- gsub( "'", "", ha_table$label[ ha_i ] )
  ha_mdr <- treeMDR( ha_i, h5_dismx )
  
  ha_mdr$ix    = ha_i
  ha_mdr$group = "0"
  ha_mdr$id    = gsub(  "^[0-9A-Za-z]+|\\|", "", ha_id )
  ha_mdr$sero  = sero[i]
  
  assign( paste0( "ha_mdr.", i ), ha_mdr )
}

# cc = rbind( ha_mdr.1, ha_mdr.2, ha_mdr.3, ha_mdr.4, ha_mdr.5 )
# max(cc$Dim_1)
# min(cc$Dim_1)
# max(cc$Dim_2)
# min(cc$Dim_2)



# H5 - n1  
h5n1_sam1 = c( -0.05, -0.03, -0.01, 0.01 )
h5n1_sam2 = c( 0.02, 0.04, -0.0125, 0 )
h5n1_sam3 = c( 0.02, 0.03, 0.03, 0.055 )

for( j in 1: 3 )
{
  sam = get( paste0( "h5n1_sam", j ) )
  y = 
  ha_mdr.1 %>% 
  filter( Dim_1 > sam[1] &  Dim_1 < sam[2] & Dim_2 > sam[3] & Dim_2 < sam[4] ) %>%
  select( ix )
  ha_mdr.1$group[ match( y$ix, ha_mdr.1$ix ) ] = as.character(j)
}

# H5 - n2
h5n2_sam1 = c( -0.06, -0.03, -0.0125, 0.01 )
h5n2_sam2 = c( -0.05, -0.04, -0.03, -0.0125 )
h5n2_sam3 = c( -0.02, 0, -0.0125, 0 )
h5n2_sam4 = c( 0, 0.015, -0.0125, 0.01 )

for( j in 1: 4 )
{
  sam = get( paste0( "h5n2_sam", j ) )
  y = 
    ha_mdr.2 %>% 
    filter( Dim_1 > sam[1] &  Dim_1 < sam[2] & Dim_2 > sam[3] & Dim_2 < sam[4] ) %>%
    select( ix )
  ha_mdr.2$group[ match( y$ix, ha_mdr.2$ix ) ] = as.character(j)
}



# H5 - n5
h5n5_sam1 = c( -0.06, -0.02, -0.0125, 0.0125 )
h5n5_sam2 = c( 0.01, 0.03, -0.0125, 0.0125 )

for( j in 1: 2 )
{
  sam = get( paste0( "h5n5_sam", j ) )
  y = 
    ha_mdr.3 %>% 
    filter( Dim_1 > sam[1] &  Dim_1 < sam[2] & Dim_2 > sam[3] & Dim_2 < sam[4] ) %>%
    select( ix )
  ha_mdr.3$group[ match( y$ix, ha_mdr.3$ix ) ] = as.character(j)
}



# H5 - n6
h5n6_sam1 = c( -0.065, -0.04, -0.03, 0 )
h5n6_sam2 = c( -0.04, -0.01, 0, 0.055 )
h5n6_sam3 = c( -0.01, 0.03, -0.0125, 0.0125 )

for( j in 1: 3 )
{
  sam = get( paste0( "h5n6_sam", j ) )
  y = 
    ha_mdr.4 %>% 
    filter( Dim_1 > sam[1] &  Dim_1 < sam[2] & Dim_2 > sam[3] & Dim_2 < sam[4] ) %>%
    select( ix )
  ha_mdr.4$group[ match( y$ix, ha_mdr.4$ix ) ] = as.character(j)
}


# H5 - n8
h5n8_sam1 = c( -0.06, -0.03, -0.0125, 0.02 )
h5n8_sam2 = c( -0.02, 0.03, -0.04, 0.0125 )

for( j in 1: 2 )
{
  sam = get( paste0( "h5n8_sam", j ) )
  y = 
    ha_mdr.5 %>% 
    filter( Dim_1 > sam[1] &  Dim_1 < sam[2] & Dim_2 > sam[3] & Dim_2 < sam[4] ) %>%
    select( ix )
  ha_mdr.5$group[ match( y$ix, ha_mdr.5$ix ) ] = as.character(j)
}



