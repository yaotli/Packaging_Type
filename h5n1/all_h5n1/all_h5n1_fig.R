source( "functions.R" )
require( tidyverse )
require( ggtree )
require( ggpubr )
require(ape)


# read in 
h5n1.tb  <- read.table( "./all_h5n1/result/20180829_h5n1_big_h5-com", sep = "\t", header = T ) 
h5n1.nwk <- "./all_h5n1/result/20180829_h5n1_big_h5-anno.nwk"
h5n1.tre <- read.beast( "./all_h5n1/result/20180829_h5n1_big_h5-anno.tre" )
nx.tre   <-  "./Nx/result/20180913_mrca_h5nx_234_h5-anno.tre"
nx.nwk   <- "./Nx/result/20180913_mrca_h5nx_234_h5-anno.nwk" 
nx.tb    <- read.table( "./Nx/result/20180913_mrca_h5nx_234_h5-com", sep = "\t", header = T)

# skygrowth 
h5n1.g <- skygrowth_df( ls.nwk = h5n1.nwk, ls.t = 2016.139, name = "allh5n1", n.res = 200)
nx.g   <- skygrowth_df( ls.nwk = nx.nwk, ls.t = 2013.989, name = "allh5n1", n.res = 25)

# plot 
# tree
h0 <- 
  ggtree( h5n1.tre, right = TRUE, mrsd = "2016-02-21", size = 0.3 ) + 
  geom_tippoint( size = 0.3 ) +
  coord_cartesian( xlim = c(1995, 2016) ) +
  scale_x_continuous( breaks = seq( 1995, 2017, by = 2 ) ) + 
  scale_y_continuous( expand = c(0,5) ) +
  theme_tree2( axis.text.x = element_blank() ) 

recdf <- data.frame( xstart = 2006.6, 
                     xend   = 2009 )
  
h0 <- h0 + geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
                      fill = "gray", alpha = 0.2, inherit.aes = FALSE)   


ht <- treelayer( origin_layer = h0, treefile = nx.tre, mostRecTime = 2013.989,
                 branchSize = 0.3, branchColor = "#ff4000", 
                 yPosition = 205, yFactor = 2, tipSize = 0.1)


# population size 
h2 <- 
  h5n1.g %>%
  filter( type == "mcmc" ) %>% 
  select( lb, e, ub, time ) %>% 
  ggplot() + theme_bw() + 
  # h5n1
  geom_line( aes( x= time, y = e ), size = 1.5, alpha = 0.6) +
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  geom_line( data = h5n1.tb, aes( x = Time, y = Median ), size = 1.5 ) + 
  geom_ribbon( data = h5n1.tb, aes( x = Time, ymin = Upper, ymax = Lower ), alpha = 0.1 ) + 
  # nx
  geom_line( data = nx.tb, aes( x = Time, y = Median ), size = 1.5, color = "darkred") + 
  geom_ribbon( data = nx.tb, aes( x = Time, ymin = Upper, ymax = Lower ), alpha = 0.1 ) + 
  geom_line( data = nx.g[ which(nx.g$type == "mcmc"), ], aes( x= time, y = e ), size = 1.5, color = "darkred", alpha = 0.6)  +
  geom_ribbon( data = nx.g[ which(nx.g$type == "mcmc"), ], aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  
  #geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = 0.0001, ymax = 10000),
  #           fill = "gray", alpha = 0.2, inherit.aes = FALSE) +  

  coord_cartesian( ylim = c(0.1, 100), xlim = c(1995, 2016) ) +
  scale_y_log10( breaks = scales::trans_breaks( "log10", function(x) 10^x),
                 labels = scales::trans_format( "log10", scales::math_format(10^.x)) ) +
  annotation_logticks() + 
  ylab( "Population size" ) + 
  scale_x_continuous( breaks = seq( 1995, 2017, by = 2)) + 
  theme( panel.grid.major.y = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.title.x = element_blank(),
         axis.text.x = element_blank() ) 

# growth rate
h3 <- 
  h5n1.g %>%
  filter( type == "rate" ) %>% 
  select( lb, e, ub, time ) %>% 
  ggplot() + theme_bw() + 
  geom_line( aes( x= time, y = e ), size = 1.5) +
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  geom_line( data = nx.g[ which(nx.g$type == "rate"), ], aes( x= time, y = e ), size = 1.5, color = "darkred") +
  geom_ribbon( data = nx.g[ which(nx.g$type == "rate"), ], aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  
  coord_cartesian( xlim = c(1995, 2016) ) +
  scale_x_continuous( breaks = seq( 1995, 2016, by = 2), 
                      labels = seq( 1995, 2016, by = 2) ) + 
  #geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
  #           fill = "gray", alpha = 0.2, inherit.aes = FALSE) + 
  ylab( "Growth rate (skygroth)" ) + 
  geom_hline( yintercept = 0, linetype = "dashed" ) +
  theme( panel.grid.major.y = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.title.x = element_blank())


ggarrange( ht, h2, h3, ncol = 1, nrow = 3, align = "v", labels = c("A", "B", "C"))
