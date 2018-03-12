require(ape)
require(ggtree)
require(stringr)

c234_h5_annotre <- "./dynamics/ride_all/result/ride_all_234_h5-anno.tre"
c234_n1_annotre <- "./dynamics/ride_all/result/ride_all_234_n1-anno.tre"
c232_h5_annotre <- "./dynamics/ride_all/result/ride_all_232_h5pika-anno.tre"
c232_n1_annotre <- "./dynamics/ride_all/result/ride_all_232_n1pika-anno.tre"
c234_h5_table <- "./eco/eco.234"
c234_n1_table <- "./eco/eco.234.n1"
c232_h5_table <- "./eco/eco.232"
c232_n1_table <- "./eco/eco.232.n1"

clade     <- c()
gene      <- c()
time      <- c()
eco       <- c()
tmrca_med <- c()
tmrca_u   <- c()
tmrca_l   <- c()

ls_annotre <- c( c234_h5_annotre, c234_n1_annotre, c232_h5_annotre, c232_n1_annotre )
ls_table   <- c( c234_h5_table, c234_n1_table, c232_h5_table, c232_n1_table)

for( i in 1:4)
{
  cl     <- str_match( ls_annotre[i], "ride_all_([0-9]+)_([a-z0-9]{2})" )[,2]
  g      <- str_match( ls_annotre[i], "ride_all_([0-9]+)_([a-z0-9]{2})" )[,3]
  tre    <- read.beast( ls_annotre[i] )
  tre.df <- fortify( tre )
  labels <- attributes(tre)$phylo$tip.label
  sp     <- read.table( ls_table[i], header = TRUE, stringsAsFactors = FALSE)$states[ match( labels, read.table( ls_table[i], header = TRUE, stringsAsFactors = FALSE)$id ) ]
  y      <- floor( as.numeric( str_match( attributes(tre)$phylo$tip.label, "_([0-9.]+)$" )[,2] ) )
  l.t    <- max( as.numeric( str_match( attributes(tre)$phylo$tip.label, "_([0-9.]+)$" )[,2] ) )
  
  tg <- c(2006)
  
  for(s in 1:2)
  {
    for( t in 1:length(tg) )
    {
     
      taxa = intersect( which( sp == c("D", "W")[s]), c( which( y == tg[t] ), which( y == tg[t]+1 ) ) )
      
      if( length(taxa) > 1 )
      {
        tmrca_med <- c(tmrca_med, l.t - as.numeric( tre.df[ MRCA( tre, labels[taxa] ), ]$height ) )
        tmrca_l   <- c(tmrca_l, l.t - as.numeric( str_match( tre.df[ MRCA( tre, labels[taxa] ), ]$height_0.95_HPD, "([0-9.]+),([0-9.]+)" )[,2] ) )
        tmrca_u   <- c(tmrca_u, l.t - as.numeric( str_match( tre.df[ MRCA( tre, labels[taxa] ), ]$height_0.95_HPD, "([0-9.]+),([0-9.]+)" )[,3] ) )
        
        clade <- c(clade, cl)
        gene  <- c(gene, g)
        time  <- c(time, tg[t])
        eco   <- c(eco, ifelse(s == 1, "D", "W") )
        
      }
    }
  }
}

df_mrca <- data.frame( clade, gene, time, eco, tmrca_med, tmrca_l, tmrca_u, stringsAsFactors = FALSE)


ggplot( data = df_mrca[ which( df_mrca$gene == "n1"), ] ) + 
  geom_point( aes( x = time, y = tmrca_med, color = eco, group = eco) ) + 
  geom_line( aes( x = time, y = tmrca_med, color = eco, group = eco) ) + 
  #scale_x_continuous( labels = seq(2005, 2015, by = 1), breaks = seq(2005, 2015, by = 1) ) +
  facet_wrap(~clade)


