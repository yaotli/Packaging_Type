source("functions.R")

require(ggplot2)
require(ggpubr)
require(tidyverse)
require(HDInterval)


# readin 

mrca_c234 = read.table( "./Nx/mrca/20180711_clade234_mrca", header = TRUE)
mrca_c234 = mrca_c234[ ,c(1,2,4,6,8)]

mrca_c7   = read.table( "./Nx/mrca/20180715_clade7_mrca", header = TRUE)
mrca_c7   = mrca_c7[ ,c(1,2,4,6,8)]


rawtb <- read.csv( "./dynamics/dynamics.csv", header = TRUE, stringsAsFactors = FALSE )
rawtb$note = as.character( rawtb$note )

# violin plot 

df.234hdi = data.frame( m = rep( c( "se", "sr", "ue", "ur" ) ), 
                        u = c( hdi( mrca_c234$se, 0.95 )[[2]], hdi( mrca_c234$sr, 0.95 )[[2]], hdi( mrca_c234$ue, 0.95 )[[2]], hdi( mrca_c234$ur, 0.95 )[[2]] ),
                        l = c( hdi( mrca_c234$se, 0.95 )[[1]], hdi( mrca_c234$sr, 0.95 )[[1]], hdi( mrca_c234$ue, 0.95 )[[1]], hdi( mrca_c234$ur, 0.95 )[[1]] ),
                        stringsAsFactors = FALSE)
v1 <- 
mrca_c234 %>% 
  gather( method, value, se:ur ) %>% 
  ggplot( aes( x = method, y = 2013.989 - value  ) ) + 
  geom_violin( ) + 
  theme_bw() +
  stat_summary( fun.y = median, geom = "point", size = 2) + 
  geom_errorbar( data = df.234hdi, aes( x = m, ymin = 2013.989-l, ymax = 2013.989-u ), inherit.aes = FALSE, width = 0) +
  scale_y_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) + 
  scale_x_discrete( labels = c( "strict+exp", "strict+ride", "relaxed+exp", "relaxed+ride"), breaks = c( "se", "sr", "ue", "ur") ) +
  theme( axis.title.x = element_blank(), 
         axis.title.y = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y  = element_blank() ) +
  coord_flip( ylim = c(2003, 2012) ) + 
  ggtitle("clade 2.3.4")
  

 
df.7hdi = data.frame( m = rep( c( "se", "sr", "ue", "ur" ) ), 
                        u = c( hdi( mrca_c7$se, 0.95 )[[2]], hdi( mrca_c7$sr, 0.95 )[[2]], hdi( mrca_c7$ue, 0.95 )[[2]], hdi( mrca_c7$ur, 0.95 )[[2]] ),
                        l = c( hdi( mrca_c7$se, 0.95 )[[1]], hdi( mrca_c7$sr, 0.95 )[[1]], hdi( mrca_c7$ue, 0.95 )[[1]], hdi( mrca_c7$ur, 0.95 )[[1]] ),
                        stringsAsFactors = FALSE)
v2 <- 
  mrca_c7 %>% 
  gather( method, value, se:ur ) %>% 
  ggplot( aes( x = method, y = 2015.038 - value  ) ) + 
  geom_violin( ) + 
  theme_bw() +
  stat_summary( fun.y = median, geom = "point", size = 2) + 
  geom_errorbar( data = df.7hdi, aes( x = m, ymin = 2015.038-l, ymax = 2015.038-u ), inherit.aes = FALSE, width = 0) +  
  scale_y_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) + 
  scale_x_discrete( labels = c( "strict+exp", "strict+ride", "relaxed+exp", "relaxed+ride"), breaks = c( "se", "sr", "ue", "ur") ) +
  theme( axis.title.x = element_blank(), 
         axis.title.y = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.y  = element_blank() ) +
  coord_flip( ylim = c(2003, 2012) )+
  ggtitle("clade 7")
                


d0 <- 
  rbind( cbind( 2013.989 - mrca_c234, clade = "234" ), cbind( 2015.038 - mrca_c7, clade = "7" ) ) %>% 
  gather( method, value, se:ur ) %>% 
  filter( method == "ur") %>% 
  select( state, value, method, clade) %>% 
  ggplot( aes(  x = value, fill = clade, color = clade) ) + 
  geom_density( alpha = 0.5 ) +
  theme_bw() +
  ylab("Estimated TMRCA") +
  coord_cartesian( xlim = c(2003, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) + 
  scale_color_manual( values = pyCol( c("red", "blue") ) ) +
  theme( axis.title.x = element_blank(), 
         axis.text.y  = element_blank(),
         axis.ticks.y = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(), 
         panel.grid.major.y = element_blank(), 
         legend.position = c(0.8,0.6)) 


# skygrowth 

l1 <- 
rawtb %>%
  filter( 2003 <= time ) %>% 
  filter( time <= 2012 ) %>% 
  filter( type == "rate" ) %>%
  filter( note == 7 | note == 234 )    %>%
  select( e, time, note )   %>%
  ggplot(  ) + theme_bw() + 
  geom_hline( yintercept = 0, linetype = "dashed") +
  geom_line( aes( x = time, y = e, color = note, group = note ), size = 2 ) + 
  xlab("") +   ylab( "Growth rate" ) +
  ggtitle(" ") +
  scale_y_continuous( breaks = 0, labels = 0, 
                      sec.axis = sec_axis( ~.*2, name = "Reported human cases") ) +
  coord_cartesian( ylim = c(-6, 9), xlim = c(2003, 2012) ) +
  scale_color_manual( values = pyCol( c("red", "blue")) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  theme( axis.title.x = element_blank(), 
         axis.text.x  = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         legend.title = element_blank(),
         legend.position = c(0.5, 0.8) )

who_cn <- data.frame( year = seq( 2003, 2012, 1), 
                      case = c( 1, 0, 8, 13, 5,
                                4, 7, 2, 1, 2 ))


l2 = l1 + geom_line( data = who_cn, aes( x = year, y = case/2) ) +
          geom_point( data = who_cn, aes( x = year, y = case/2) ) 
  



# combine 

ggarrange( l2, d0, nrow = 2, heights = c(1,0.8), align = "v" ) #4.5*4
ggarrange( v1, v2, ncol = 2, labels = c("A","B") ) #a6




