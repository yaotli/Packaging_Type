source("functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(ggpubr)


rawtre.c234 <- "./bigtree/raw/raxml_pH5_c234_1699_e0202.tre"
rawtre.c232 <- "./bigtree/raw/raxml_pH5_c232_1007_e0202.tre"
refcsv.c234 <- read.csv( "./state/eco_234.csv", header = TRUE, stringsAsFactors = FALSE )
refcsv.c232 <- read.csv( "./state/eco_232.csv", header = TRUE, stringsAsFactors = FALSE )

# c234 ml 
tre.c234   <- read.nexus( rawtre.c234 )
refdf.c234 <- data.frame( name  = gsub("'", "", fortify(tre.c234)$label[ 1: length(tre.c234$tip.label) ] ),
                          state = NA, stringsAsFactors = FALSE)

refdf.c234$state                            <- refcsv.c234$states[ match( refdf.c234$name, refcsv.c234$name ) ]
refdf.c234$state[ is.na(refdf.c234$state) ] <- "nonCHN"

df.c234 <- findtaxa2( tree     = tre.c234, 
                      targetid = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "U_e", "nonCHN"), 
                      target   = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "U_e", "nonCHN"),
                      Refdf    = refdf.c234 )

df.c234$colorr <- gsub( "black", "U_e", df.c234$colorr )
df.c234$colorr <- factor( x      =  df.c234$colorr, 
                          levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "U_e", "nonCHN")  )

t.234 <- 
ggtree( tre.c234, ladderize = FALSE, size = 0.4 ) %<+% df.c234 + aes( color = colorr ) + 
  geom_treescale( width = 0.005, offset = -25, x = 0.045, y = 455) + 
  scale_y_continuous( expand = c(0,5) ) + 
  scale_x_continuous( expand = c(0,0.0005) ) + 
  scale_color_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#01665e", "black", "gray"), 
                      labels = c("Domestic avian", "Domestic mammal", "Domestic environ.", 
                                 "Wild environ.", "Wild mammal", "Wild avian", "non-China")) + 
  theme( legend.position = "none", legend.title = element_blank() )
  
# c234 dist.
table_234_species        <- as.data.frame( prop.table( table( refcsv.c234[,3], refcsv.c234[,5]), margin = 2) )
table_234_species$Var1   <- factor( table_234_species$Var1, levels = c("D_a", "D_m", "D_e", "U_e","W_e", "W_a") )
table_234_species$Var2   <- as.numeric( as.character(table_234_species$Var2) )

d.234 <- 
  ggplot( data =  table_234_species, aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2003.5, 2016.5), breaks = seq(2004, 2016, by = 2), labels = seq(2004, 2016, by = 2),
                      expand = c(0,0.1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% China Seq.") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "white", "#c7eae5", "#01665e") ) +
  theme( legend.title    = element_blank(), 
         legend.position = "none", 
         axis.text.y     = element_text(size = 8), 
         axis.text.x     = element_text(size = 10))




# c232 ml 
tre.c232   <- read.nexus( rawtre.c232 )
refdf.c232 <- data.frame( name  = gsub("'", "", fortify(tre.c232)$label[ 1: length(tre.c232$tip.label) ] ),
                          state = NA, stringsAsFactors = FALSE)

refdf.c232$state                            <- refcsv.c232$states[ match( refdf.c232$name, refcsv.c232$name ) ]
refdf.c232$state[ is.na(refdf.c232$state) ] <- "nonCHN"

df.c232 <- findtaxa2( tree     = tre.c232, 
                      targetid = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "U_e", "nonCHN"), 
                      target   = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "U_e", "nonCHN"),
                      Refdf    = refdf.c232 )
df.c232$colorr <- gsub( "black", "U_e", df.c232$colorr )
df.c232$colorr <- factor( x =  df.c232$colorr, 
                          levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a", "U_e", "nonCHN")  )

t.232 <- 
ggtree( tre.c232, ladderize = FALSE, size = 0.4) %<+% df.c232 + aes( color = colorr ) + 
  geom_treescale( width = 0.005, offset = -15, x = 0.045, y = 260) + 
  scale_y_continuous( expand = c(0,5) ) + 
  scale_x_continuous( expand = c(0,0.0005) ) + 
  scale_color_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "black", "gray") ) + 
  theme( legend.position = "none", legend.title = element_blank() )

# c232 dist.
table_232_species        <- as.data.frame( prop.table( table( refcsv.c232[,3], refcsv.c232[,5] ), margin = 2) )
table_232_species$Var1   <- factor( table_232_species$Var1, levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a") )
table_232_species$Var2   <- as.numeric( as.character(table_232_species$Var2) )

d.232 <- 
  ggplot( data =  table_232_species, aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2003.5, 2016.5), breaks = seq(2004, 2016, by = 2), labels = seq(2004, 2016, by = 2),
                      expand = c(0,0.1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% China Seq.") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e") ) + 
  theme( legend.title    = element_blank(), 
         legend.position = "none", 
         axis.text.y     = element_text(size = 8), 
         axis.text.x     = element_text(size = 10) )



# combine 
ggarrange( t.232, t.234, d.232, d.234, ncol = 2, nrow = 2, heights = c(6,1) ) #6*8 inch

