source("functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(ggpubr)


rawtre.c234 <- "./bigtree/raw/raxml_pH5_c234_1699_e0202.tre"
rawtre.c232 <- "./bigtree/raw/raxml_pH5_c232_1007_e0202.tre"
rawtre.c7   <- "./clade7/pH5_c7_80_rmdP_e0802.tre"

refcsv.c234 <- read.csv( "./state/eco_234.csv", header = TRUE, stringsAsFactors = FALSE )
refcsv.c232 <- read.csv( "./state/eco_232.csv", header = TRUE, stringsAsFactors = FALSE )
refcsv.c7   <- read.csv( "./state/eco_7.csv", header = TRUE, stringsAsFactors = FALSE )


# c234 ml ------

tre.c234   <- read.nexus( rawtre.c234 )
refdf.c234 <- fortify( tre.c234 )

refdf.c234[ ncol(refdf.c234) + 1 ]         <- refcsv.c234$states[ match( gsub( "'", "", refdf.c234$label ), refcsv.c234$name ) ]
colnames( refdf.c234 )[ ncol(refdf.c234) ] <- "eco"
refdf.c234$eco[ is.na(refdf.c234$eco) & refdf.c234$isTip ] <- "nonCN"
refdf.c234$eco[ is.na(refdf.c234$eco) ] <- "internal"
refdf.c234$eco                             <- factor( x      = refdf.c234$eco,
                                                      levels = c("D_a", "D_m", "D_e", "W_e", "W_a", "U_e", "nonCN", "internal") )

#n = refdf.c234$label[ refdf.c234$isTip ]
#subtype.c234 <- data.frame( subtype = str_match( n, "\\|_(H5N[1-9])_"  )[,2] )
#rownames( subtype.c234 ) = n

m1 = 
ggtree( tre.c234, ladderize = FALSE, size = 0.4 ) %<+% refdf.c234 + aes( color = eco ) +
  geom_tippoint( aes( fill = eco, shape = eco, color = eco), size = 0.7) +
  scale_shape_manual( values = c( rep( 21, 6), NA) ) +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#01665e", "black", "gray" ) ) +
  scale_color_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#01665e", "black", "gray", "black") ) +
  geom_treescale( width = 0.005, offset = -25, x = 0.045, y = 455) + 
  scale_y_continuous( expand = c(0,5) ) +
  scale_x_continuous( expand = c(0,0.0005) ) + 
  theme_tree( legend.position = c(0.15,0.35), 
              legend.title = element_blank() ) +
  ggtitle("clade 2.3.4")



# c232 ml ------

tre.c232   <- read.nexus( rawtre.c232 )
refdf.c232 <- fortify( tre.c232 )

refdf.c232[ ncol(refdf.c232) + 1 ]         <- refcsv.c232$states[ match( gsub( "'", "", refdf.c232$label ), refcsv.c232$name ) ]
colnames( refdf.c232 )[ ncol(refdf.c232) ] <- "eco"
refdf.c232$eco[ is.na(refdf.c232$eco) & refdf.c232$isTip ] <- "nonCN"
refdf.c232$eco[ is.na(refdf.c232$eco) ] <- "internal"
refdf.c232$eco                             <- factor( x      = refdf.c232$eco,
                                                      levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "nonCN", "internal") )


m2 = 
  ggtree( tre.c232, ladderize = FALSE, size = 0.4 ) %<+% refdf.c232 + aes( color = eco ) +
  geom_tippoint( aes( fill = eco, shape = eco, color = eco), size = 0.7) +
  scale_shape_manual( values = c( rep( 21, 6), NA) ) +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "gray" ) ) +
  scale_color_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "gray", "black") ) +
  geom_treescale( width = 0.005, offset = -15, x = 0.045, y = 260) + 
  scale_y_continuous( expand = c(0,5) ) + 
  scale_x_continuous( expand = c(0,0.0005) ) + 
  theme_tree( legend.position = c(0.15,0.35), 
              legend.title = element_blank() ) +
  ggtitle("clade 2.3.2")



# c7 ml ------

tre.c7   <- read.nexus( rawtre.c7 )
refdf.c7 <- fortify( tre.c7 )

refdf.c7[ ncol(refdf.c7) + 1 ]         <- refcsv.c7$states[ match( gsub( "'", "", refdf.c7$label ), refcsv.c7$name ) ]
colnames( refdf.c7 )[ ncol(refdf.c7) ] <- "eco"
refdf.c7$eco[ is.na(refdf.c7$eco) & refdf.c7$isTip ] <- "nonCN"
refdf.c7$eco[ is.na(refdf.c7$eco) ] <- "internal"
refdf.c7$eco                        <- factor( x     = refdf.c7$eco,
                                              levels = c("D_a", "D_m", "D_e", "W_e", "W_a", "nonCN", "internal") )


m3 = 
  ggtree( tre.c7, ladderize = FALSE, size = 0.4 ) %<+% refdf.c7 + aes( color = eco ) +
  geom_tippoint( aes( fill = eco, shape = eco, color = eco), size = 0.7) +
  scale_shape_manual( values = c( rep( 21, 5), NA) ) +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#01665e", "gray" ) ) +
  scale_color_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#01665e", "gray", "black") ) +
  geom_treescale( width = 0.005, offset = -3) + 
  scale_y_continuous( expand = c(0,5) ) + 
  scale_x_continuous( expand = c(0,0.0005) ) + 
  theme_tree( legend.position = "NULL" ) +
  ggtitle("clade 7")



# c234 dist. ------
table_234_species        <- as.data.frame( prop.table( table( refcsv.c234[,3], refcsv.c234[,5]), margin = 2) )
table_234_species$Var1   <- factor( table_234_species$Var1, levels = c("D_a", "D_m", "D_e", "U_e","W_e", "W_a") )
table_234_species$Var2   <- as.numeric( as.character(table_234_species$Var2) )

d.234 <- 
  ggplot( data =  table_234_species, aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2003.5, 2016.5), 
                      breaks = seq(2004, 2016, by = 2), 
                      labels = seq(2004, 2016, by = 2),
                      expand = c(0,0.1)) +
  scale_y_continuous( breaks = c(0, 0.5, 1), 
                      labels = c(0, 0.5, 1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Seq.") + #ggtitle("clade 2.3.4") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "white", "#c7eae5", "#01665e") ) +
  theme( legend.title    = element_blank(), 
         legend.position = "none")


# c232 dist. ------
table_232_species        <- as.data.frame( prop.table( table( refcsv.c232[,3], refcsv.c232[,5] ), margin = 2) )
table_232_species$Var1   <- factor( table_232_species$Var1, levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a") )
table_232_species$Var2   <- as.numeric( as.character(table_232_species$Var2) )

d.232 <- 
  ggplot( data =  table_232_species, aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2003.5, 2016.5), 
                      breaks = seq(2004, 2016, by = 2), 
                      labels = seq(2004, 2016, by = 2),
                      expand = c(0,0.1)) +
  scale_y_continuous( breaks = c(0, 0.5, 1), 
                      labels = c(0, 0.5, 1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Seq.") + #ggtitle("clade 2.3.2") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e") ) + 
  theme( legend.title    = element_blank(), 
         legend.position = "none")



# c232 dist. ------
table_7_species        <- as.data.frame( prop.table( table( refcsv.c7[,3], refcsv.c7[,5] ), margin = 2) )
table_7_species$Var1   <- factor( table_7_species$Var1, levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a") )
table_7_species$Var2   <- as.numeric( as.character(table_7_species$Var2) )

d.7 <- 
  ggplot( data =  table_7_species, aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2005, 2016), 
                      breaks = seq(2004, 2015, by = 2), 
                      labels = seq(2004, 2015, by = 2),
                      expand = c(0,0.1)) +
  scale_y_continuous( breaks = c(0, 0.5, 1), 
                      labels = c(0, 0.5, 1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Seq.")  + #ggtitle("clade 7") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#01665e") ) + 
  theme( legend.title    = element_blank(), 
         legend.position = "none" )



# combine 
ggarrange( m2, m1, ggarrange( m3,nrow = 2, heights = c(0.6, 1)), widths = c(1,1,0.8), 
           d.232, d.234, d.7, ncol = 3, 
           nrow = 2,  heights = c(6,1) ) #a5



ggarrange( t.232, t.234, d.232, d.234, ncol = 2, nrow = 2, heights = c(6,1) ) #6*8 inch

