source("functions.R")

require(ggtree)
require(ggridges)
require(tidyverse)
require(extrafont)
require(ggpubr)

gsgd_rmdptre      <- taxaInfo( "./clade/pH5_7326_gsgd_rmdP_e0719.tre", useTree = TRUE )
gsgd_rmdptre[[7]][ is.na( gsgd_rmdptre[[7]] ) ] <- "others"
gsgd_rmdptre[[7]] <- c( "others", "234", "232", "7")[ match( gsgd_rmdptre[[7]], unique( gsgd_rmdptre[[7]] ) ) ]

table0 <- data.frame( sero = gsgd_rmdptre[[3]], geo  = gsgd_rmdptre[[2]], 
                      yr   = gsgd_rmdptre[[4]], anno = gsgd_rmdptre[[7]], 
                      name = gsgd_rmdptre[[6]], stringsAsFactors = FALSE )

tb = 
  table0                                                                              %>% 
  select( sero, geo, yr, anno )                                                       %>%
  mutate( China = ifelse( geo == "China" | geo == "Hong_Kong", "China", "others" ) )  %>%
  mutate( N1 = ifelse( sero == "H5N1", "N1", "Nx" ) ) %>%
  select( N1, yr, anno, China )                                                       


# ggridge --------

tb$yr    <- floor( tb$yr )
tbb      <- data.frame( table(tb), stringsAsFactors = FALSE)
tbb$yr   <- as.numeric( as.character( tbb$yr ) )
tbb$anno <- factor( tbb$anno, levels = c("others", "7", "232", "234" ) )

s1 <- 
tbb %>% 
  filter( China == "China" ) %>%
  ggplot( aes( x = yr, y = anno, height = Freq, fill = anno, linetype =  N1, color = anno) ) +
  geom_density_ridges( stat = "identity", scale = 5.2, size = 1.2, alpha = .4) + 
  
  scale_x_continuous( limit  = c(1996, 2017), 
                      breaks = seq(1996, 2016, by = 2), labels = seq(1996, 2016, by = 2),
                      expand = c(0.001, 0)) +
  
  scale_y_discrete( expand = c(0.1, 0) ) + 
  scale_fill_manual( values = c( pyCol( c( "gray", "blue", "green", "red") ) ) ) + 
  scale_color_manual( values = c( pyCol( c( "gray", "blue", "green", "red") ) ) ) + 
  theme_minimal() +
  theme( text = element_text(family = "Helvetica"), 
         panel.grid.minor = element_blank(),
         strip.background = element_rect( fill = "white", color = "white "), 
         legend.title = element_blank(),
         legend.position = c(0.1,0.6) ) +
  labs(x = "", y = "")

# line plot --------  
  
tbp      <- tb
tbp$anno <- ifelse( tbp$anno == "232" | tbp$anno == "234" | tbp$anno == "7", "23", "others")
tbp      <- tbp[ which( tbp$China == "China" ), ]
tbp <- tbp[ ,c(2,3)]

tbp    <- as.data.frame( prop.table( table( tbp ), margin = 1)  )
tbp$yr <- as.numeric( as.character(tbp$yr) )

s2 <- 
ggplot( data = tbp ) + 
  geom_line(aes( x = yr, y = Freq, linetype = anno, group = anno ), size = 1.5, color = "darkblue") + 
  labs(x = "", y = "") +   theme_minimal() +
  scale_x_continuous( limit = c(1996, 2017),
                      breaks = seq(1996, 2016, by = 2), labels = seq(1996, 2016, by = 2), 
                      expand = c(0, 0.001) ) +
  
  scale_y_continuous( expand = c(0.1, 0) ) + 
  theme( panel.grid.minor   = element_blank(), 
         panel.grid.major.x = element_blank(), 
         legend.title = element_blank(),
         legend.position = c(0.2,0.5) )


# gsgd tree --------  

gsgd_rawtr <- read.tree("clade/pH5_7326_gsgd_rmdP_e0627.nwk" )
trefile    <- fortify( gsgd_rawtr )

# geom_text(aes(label = node), size = 0.5)
c234   <- c( 4777, getDes( tre.d = trefile, node = 4777) )
c232   <- c( 6396, getDes( tre.d = trefile, node = 6396) )
c7     <- c( 9284, getDes( tre.d = trefile, node = 9284) )
geo_CN <- grep( "China|Hong_Kong", trefile$label )

c234.CN <- intersect( c234, geo_CN  )
c232.CN <- intersect( c232, geo_CN  )
c7.CN   <- intersect( c7, geo_CN  )
o.CN    <- setdiff( geo_CN, c( c234, c232, c7 ) )

trefile[, 10 ] <- "others"
trefile[, 11 ] <- 0.3

colnames( trefile )[ 10 ] <- "branch.col"
colnames( trefile )[ 11 ] <- "branch.size"

trefile$branch.col[ c234 ]    <- "c234"
trefile$branch.col[ c234.CN ] <- "CN_c234"
trefile$branch.col[ c232 ]    <- "c232"
trefile$branch.col[ c232.CN ] <- "CN_232"
trefile$branch.col[ c7 ]      <- "c7"
trefile$branch.col[ c7.CN ]   <- "CN_7"
trefile$branch.col[ o.CN ]    <- "CN_others"

trefile$branch.size[ c(c234, c232, c7) ] <- 0.6

s3 <- 
ggtree( gsgd_rawtr ) %<+% trefile + aes( color = branch.col, size = I(branch.size) ) + 
  
  theme_tree( legend.position = c(0.15,0.85), 
              legend.title = element_blank(), 
              text = element_text(family = "Helvetica")) +
  
  geom_treescale( wid = 0.01, x = 0, y = 500, offset = 100) +
  scale_y_continuous( expand = c(0.01, 0) ) +
  scale_x_continuous( expand = c(0.01, 0) ) +
  geom_tippoint( aes( fill = branch.col, shape = branch.col, color = branch.col) ) +
  
  scale_fill_manual( values = c("#9be39b", "#eb9394", "#92c7ed",
                                 "#1c641c","#12476d", "#6c1414", 
                                 "#404040", "#bfbfbf" ) ) +
  scale_color_manual( values = c("#9be39b", "#eb9394", "#92c7ed",
                                "#1c641c","#12476d", "#6c1414", 
                                "#404040", "#bfbfbf" ) ) +
  scale_shape_manual( values = rep( 21, 8) )


# combine
ggarrange( s3,  
           ggarrange( s1, s2, ncol = 1, nrow = 2, heights = c(2.5,1), align = "v", labels = c("B", "C")),
           ncol = 2, widths = c(3,5), labels = "A" ) #a5 landscape



# subtype tree --------  

gsgd_rawtr  <- read.tree("clade/pH5_7326_gsgd_rmdP_e0627.nwk" )
trefile.sub <- fortify( gsgd_rawtr )

trefile.sub[ ncol( trefile.sub ) + 1 ] <- str_match( trefile.sub$label, pattern = "_(H5N[0-9])_[12]" )[,2]
colnames( trefile.sub )[ ncol( trefile.sub ) ] = "subtype"

s <- 
ggtree( gsgd_rawtr ) %<+% trefile.sub + 
  
  theme_tree( legend.position = c(0.15,0.78), 
              legend.title = element_blank(), 
              text = element_text(family = "Helvetica")) +
  
  geom_tippoint( aes( fill = subtype, shape = subtype, color = subtype) )+
  geom_treescale( wid = 0.01, x = 0, y = 500, offset = 100) +
  scale_y_continuous( expand = c(0.01, 0) ) +
  scale_x_continuous( expand = c(0.01, 0) ) +
  scale_fill_manual( values = c( "#7f6e85", "#e1c72f", "#48a365", "#ccc197", "#77bedb", 
                                  "#cc79a7", "#d0694a" ) ) + 
  scale_colour_manual( values = c( "#7f6e85", rep( "black", 6)  ) ) +
  scale_shape_manual( values = c( 21, 21, 21, 21, 21, 21, 21 ) ) + 
  geom_cladelabel( node = 9284, label = "7",  barsize = 1, angle = 270, align = T, color = pyCol("blue")) +
  geom_cladelabel( node = 4777, label = "2.3.4",  barsize = 0.6, angle = 270, align = T, color = pyCol("red")) +
  geom_cladelabel( node = 6396, label = "2.3.2",  barsize = 0.6, angle = 270, align = T, color = pyCol("green")) 
#a6por


sub7_rawtr <- read.tree("clade/sub7_pH5_7326_gsgd_e0702.nwk" )
trefile.c7 <- fortify( sub7_rawtr )

trefile.c7[ ncol( trefile.c7 ) + 1 ] <- str_match( trefile.c7$label, pattern = "_(H5N[0-9])_[12]" )[,2]
colnames( trefile.c7 )[ ncol( trefile.c7 ) ] = "subtype"

s.7 <- 
  ggtree( sub7_rawtr ) %<+% trefile.c7 + 
  
  theme_tree( legend.title = element_blank(), 
              text = element_text(family = "Helvetica")) +
  
  geom_tippoint( aes( fill = subtype, shape = subtype, color = subtype) )+
  scale_x_continuous( expand = c(0.2, 0) ) +
  scale_fill_manual( values = c( "#7f6e85", "#e1c72f" ) ) + 
  scale_colour_manual( values = c( "#7f6e85", "black" ) ) +
  scale_shape_manual( values = c( 21, 21 ) )  +
  geom_cladelabel( node = 81, label = "7",  barsize = 1, angle = 270, align = T, color = pyCol("blue")) 
  
  
#a6por 


ggarrange( s, s.7, nrow = 2, heights = c(1,0.15), labels = c( "A", "B" )) 
  