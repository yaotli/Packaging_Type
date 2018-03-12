source("./functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(ggpubr)


raw_pH5_c234 <- "./dynamics/pH5_c234_159.fasta"
raw_pN1_c234 <- "./dynamics/pN1_c234_152.fasta"
raw_pH5_c232 <- "./dynamics/pH5_c232_204.fasta"
raw_pN1_c232 <- "./dynamics/pN1_c232_189.fasta"
csv_eco_c234 <- "./state/eco_234.csv"
csv_eco_c232 <- "./state/eco_232.csv"

#c234
sub_pH5_c234_159        <- taxaInfo( file = raw_pH5_c234, useTree = FALSE )
sub_pH5_c234_159[[ 8 ]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5_c234_159[[6]], read.csv( csv_eco_c234, stringsAsFactors = FALSE )$name ) ]

sub_pH5_c234_159[[ 8 ]][ is.na( sub_pH5_c234_159[[ 8 ]] ) ] <- "W_a"
sub_pH5_c234_159[[ 8 ]] <- ifelse( startsWith( sub_pH5_c234_159[[ 8 ]], "D" ), "D", "W" )

eco.234 <- data.frame( id = sub_pH5_c234_159[[ 6 ]], states = sub_pH5_c234_159[[ 8 ]], stringsAsFactors = FALSE )
write.table( x = eco.234, file = "./eco/eco.234", sep = "\t", quote = FALSE, row.names = FALSE )

#c234_n1
sub_pN1_c234_152      <- taxaInfo( file = raw_pN1_c234, useTree = FALSE )
sub_pN1_c234_152[[8]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE)$states[ 
  match( gsub( pattern = "^[A-Z0-9]+", "", sub_pN1_c234_152[[6]] ), 
         gsub( pattern = "^[A-Z0-9]+", "", read.csv( csv_eco_c234, stringsAsFactors = FALSE)$name ) ) ]

sub_pN1_c234_152[[8]][is.na(sub_pN1_c234_152[[8]])] <- c( "D_a", "W_a", "D_e", "D_e", "D_m",
                                                          "W_a", "D_a", "D_a", "D_a", "D_a",
                                                          "D_a", "D_a", "D_a", "W_a", "D_a" )

sub_pN1_c234_152[[8]] <- ifelse( startsWith( sub_pN1_c234_152[[8]], "D"), "D", "W")
eco.234.n1            <- data.frame( id = sub_pN1_c234_152[[6]], states = sub_pN1_c234_152[[8]], stringsAsFactors = FALSE )
write.table( x = eco.234.n1, file = "./eco/eco.234.n1", sep = "\t", quote = FALSE, row.names = FALSE )




#c232 
sub_pH5_c232_204        <- taxaInfo( file = raw_pH5_c232, useTree = FALSE )
sub_pH5_c232_204[[ 8 ]] <- read.csv( csv_eco_c232, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5_c232_204[[6]], read.csv( csv_eco_c232, stringsAsFactors = FALSE )$name ) ]

sub_pH5_c232_204[[ 8 ]][ is.na( sub_pH5_c232_204[[ 8 ]] ) ] <- "W_a"
sub_pH5_c232_204[[ 8 ]] <- ifelse( startsWith( sub_pH5_c232_204[[ 8 ]], "D" ), "D", "W" )

eco.232 <- data.frame( id = sub_pH5_c232_204[[ 6 ]], states = sub_pH5_c232_204[[ 8 ]], stringsAsFactors = FALSE )
write.table( x = eco.232, file = "./eco/eco.232", sep = "\t", quote = FALSE, row.names = FALSE )

eco.232_2012 <- eco.232[ which( floor( as.numeric( str_match( eco.232$id, "[0-9.]+$" ) ) ) < 2012 ), ]
write.table( x = eco.232_2012, file = "./eco/eco.232_2012", sep = "\t", quote = FALSE, row.names = FALSE )

#c232_n1
sub_pN1_c232_189      <- taxaInfo( file = raw_pN1_c232, useTree = FALSE )
sub_pN1_c232_189[[8]] <- read.csv( csv_eco_c232, stringsAsFactors = FALSE)$states[ 
  match( gsub( pattern = "^[A-Z0-9]+", "", sub_pN1_c232_189[[6]] ), 
         gsub( pattern = "^[A-Z0-9]+", "", read.csv( csv_eco_c232, stringsAsFactors = FALSE)$name ) ) ]

sub_pN1_c232_189[[ 8 ]][ is.na(sub_pN1_c232_189[[8]]) ] <- c( "D_a", "D_a", "D_a", "D_a", "W_a",
                                                              "D_a", "D_a", "D_a", "D_a", "W_a", 
                                                              "D_a", "D_e", "W_a" )

sub_pN1_c232_189[[8]] <- ifelse( startsWith( sub_pN1_c232_189[[8]], "D"), "D", "W")
eco.232.n1            <- data.frame( id = sub_pN1_c232_189[[6]], states = sub_pN1_c232_189[[8]], stringsAsFactors = FALSE )
write.table( x = eco.232.n1, file = "./eco/eco.232.n1", sep = "\t", quote = FALSE, row.names = FALSE )




# eco-stratified c234
ac_c234h5_D <- acSearch( faslist = sub_pH5_c234_159, keyword.dir = 8, keyword = "D")
ac_c234h5_W <- acSearch( faslist = sub_pH5_c234_159, keyword.dir = 8, keyword = "W")
subfastaSeq( AC = TRUE, AC_list = ac_c234h5_D, filedir = raw_pH5_c234 )
subfastaSeq( AC = TRUE, AC_list = ac_c234h5_W, filedir = raw_pH5_c234 )


# eco-stratified c232
rmDup( fasfile = raw_pH5_c232, year = c( 2000, 2012 ), rmdup = FALSE )
rmDup( fasfile = raw_pN1_c232, year = c( 2000, 2012 ), rmdup = FALSE )

ac_c232h5_D <- acSearch( faslist = sub_pH5_c232_204, keyword.dir = 8, keyword = "D", range = c( 2000, 2011 ), range.dir = 4 )
ac_c232h5_W <- acSearch( faslist = sub_pH5_c232_204, keyword.dir = 8, keyword = "W", range = c( 2000, 2011 ), range.dir = 4 )

subfastaSeq( AC = TRUE, AC_list = ac_c232h5_D, filedir = "./dynamics/ride_2012/pH5_c232_204_2012.fasta" )
subfastaSeq( AC = TRUE, AC_list = ac_c232h5_W, filedir = "./dynamics/ride_2012/pH5_c232_204_2012.fasta" )




## figures TREE ----------------

eco_tre.232  <- "./eco/eco_all/result/eco_all_232_h5pika-anno.tre"
eco_tre.234  <- "./eco/eco_all/result/eco_all_234_h5-anno.tre"
rawbeast.232 <- read.beast( eco_tre.232 )
rawbeast.234 <- read.beast( eco_tre.234 )


g1 <- 
  ggtree( rawbeast.232, right = TRUE, size = 0.6, mrsd = "2016-04-12") + 
  aes( color = states ) +
  scale_color_manual( values = c( "#8c564b", "#17becf" )  ) + 
  theme_tree2( axis.text.x = element_text( size = 12  ) ) + 
  scale_y_continuous( expand = c(0,1) ) + 
  # ecom_tippoint(aes(fill = eco), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2016.5, by = 2), 
                      labels = seq(2002, 2016, by = 2), 
                      limits = c(2001, 2017))  + 
  theme( axis.ticks  = element_blank() )

rectdf <- data.frame( xstart = seq( 2002, 2017, 2), 
                      xend   = seq( 2003, 2017, 2))

g11 <- g1 + geom_rect( data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                       fill = "gray", alpha = 0.2, inherit.aes=FALSE)


g2 <- 
  ggtree( rawbeast.234, right = TRUE, size = 0.6, mrsd = "2011-12-15") + 
  aes( color = states ) +
  scale_color_manual( values = c( "#8c564b", "#17becf" )  ) + 
  theme_tree2( axis.text.x = element_text( size = 12  ) ) + 
  scale_y_continuous( expand = c(0,1) ) + 
  # ecom_tippoint(aes(fill = eco), shape = 21, color = "black") + 
  scale_x_continuous( breaks = seq(2002.5, 2016.5, by = 2), 
                      labels = seq(2002, 2016, by = 2), 
                      limits = c(2001, 2017))  + 
  theme( axis.ticks.x  = element_blank() ) 


rectdf <- data.frame( xstart = seq( 2002, 2017, 2), 
                      xend   = seq( 2003, 2017, 2))

g22 <- g2 + geom_rect( data = rectdf, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), 
                       fill = "gray", alpha = 0.2, inherit.aes=FALSE)



## figures TRUNK ----------------


trunk.232 = read.table( "./eco/PACT/out.eco_all_232_h5pika-S1000.skylines", header = T, stringsAsFactors = F )
trunk.234 = read.table( "./eco/PACT/out.eco_all_234_h5-S1000.skylines", header = T, stringsAsFactors = F )


t1 <- 
  ggplot( data = trunk.232, aes( x = time, y = mean )) + 
  geom_area( aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2003,2015, by = 2), limits = c(2004.1, 2016), 
                      expand = c(0, 0) ) +
  scale_y_continuous( expand = c(0,0) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks       = element_blank(), 
         panel.border     = element_blank(), 
         axis.title       = element_text(size = 12),
         axis.text.y      = element_text(size = 8),
         axis.text.x      = element_text(size = 10),
         #axis.text.x = element_blank(),
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )



t2 <- 
  ggplot( data = trunk.234, aes( x = time, y = mean )) + 
  geom_area( size = 1, aes( fill = statistic) ) + 
  xlab("") + ylab("Trunk %") +
  scale_x_continuous( breaks = seq(2003,2015, by = 2), limits = c(2004.1, 2016), 
                      expand = c(0, 0) ) +
  scale_y_continuous( expand = c(0,0) ) +
  theme_bw() + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks       = element_blank(), 
         panel.border     = element_blank(), 
         axis.title       = element_text(size = 12),
         axis.text.y      = element_text(size = 8),
         axis.text.x      = element_text(size = 10),
         legend.position  = "none") + 
  scale_fill_manual( values = c( "#8c564b", "#17becf" ) )

ggarrange( g11, g22, t1, t2, ncol = 2, nrow = 2, heights = c(3.5,1) ) #8*8



## figures Persistence ----------------

per_232 <- read.table("./eco/PACT/out.eco_2012_232_h5pika-S1000.persist", header = TRUE)
per_234 <- read.table("./eco/PACT/out.eco_all_234_h5-S1000.persist", header = TRUE)

# 232
p1 <- 
  ggplot(data = per_232[-1,], aes( x= statistic, y = mean, color = statistic)) +
  geom_point(size = 5) +
  geom_errorbar( aes( ymin = lower, ymax = upper ), width = 0, size = 1) +
  coord_flip() + 
  xlab("") + ylab("Persistence (year)") +
  scale_color_manual( values = c( "#8c564b", "#17becf" ) ) + 
  scale_x_discrete( labels = c("Domestic", "Wild") ) +
  scale_y_continuous( limits = c(0,3.5) ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y     = element_blank(),
         panel.border     = element_blank(),
         axis.text.y      = element_text(size = 10),
         axis.text.x      = element_text(size = 8),
         axis.title.x     = element_text(size = 8),
         legend.position  = "none") 

# 234
p2 <- 
  ggplot(data = per_234[-1,], aes( x= statistic, y = mean, color = statistic)) +
  geom_point(size = 5) +
  geom_errorbar( aes( ymin = lower, ymax = upper ), width = 0, size = 1) +
  coord_flip() + 
  xlab("") + ylab("Persistence (year)") +
  scale_color_manual( values = c( "#8c564b", "#17becf" ) ) + 
  scale_x_discrete( labels = c("Domestic", "Wild") ) +
  scale_y_continuous( limits = c(0,3.5) ) +
  theme_bw() + #ggtitle("Inferred MRCA of Nx viruses") + 
  theme( panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(), 
         axis.ticks.y     = element_blank(),
         panel.border     = element_blank(),
         axis.text.y      = element_text(size = 10),
         axis.text.x      = element_text(size = 8),
         axis.title.x     = element_text(size = 8),
         legend.position  = "none") 

multiplot( p1, p2, ncol = 1) #4*4



















