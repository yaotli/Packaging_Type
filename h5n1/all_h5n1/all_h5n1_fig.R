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

write.csv( h5n1.g, file = "./all_h5n1/h5n1.g.csv", row.names = FALSE)
write.csv( nx.g, file = "./all_h5n1/nx.g.csv", row.names = FALSE)

h5n1.g <- read.csv( "./all_h5n1/h5n1.g.csv", header = TRUE, stringsAsFactors = FALSE )
nx.g   <- read.csv( "./all_h5n1/nx.g.csv", header = TRUE, stringsAsFactors = FALSE )



# plot 
# tree
h0 <- 
  ggtree( h5n1.tre, right = TRUE, mrsd = "2016-02-21", size = 0.3 ) + 
  geom_tippoint( size = 1 ) +
  coord_cartesian( xlim = c(1995, 2016) ) +
  scale_x_continuous( breaks = seq( 1995, 2016, by = 5), 
                      labels = seq( 1995, 2016, by = 5) ) + 
  scale_y_continuous( expand = c(0,5) ) +
  theme_tree2() 

#recdf <- data.frame( xstart = 2006.6, 
#                     xend   = 2009 )
# h0 <- h0 + geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
#                       fill = "gray", alpha = 0.2, inherit.aes = FALSE)   

h0 <- 
h0 + geom_vline( xintercept = 2006.3, linetype = "dashed", color = "gray") + 
     geom_vline( xintercept = 2009, linetype = "dashed", color = "gray")


ht <- treelayer( origin_layer = h0, treefile = nx.tre, mostRecTime = 2013.989,
                 branchSize = 0.3, branchColor = "red", 
                 yPosition = 200, yFactor = 3, tipSize = 0.5, tipColor = "black", pseudoRoot = 1)


# population size 
h2 <- 
  h5n1.g %>%
  filter( type == "mcmc" ) %>% 
  select( lb, e, ub, time ) %>% 
  ggplot() + theme_bw() + 
  # h5n1
  geom_line( aes( x= time, y = e ), size = 1.5, alpha = 0.5, linetype = "dashed") +
  #geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  geom_line( data = h5n1.tb, aes( x = Time, y = Median ), size = 1.5 ) + 
  #geom_ribbon( data = h5n1.tb, aes( x = Time, ymin = Upper, ymax = Lower ), alpha = 0.1 ) + 
  # nx
  #geom_line( data = nx.tb, aes( x = Time, y = Median ), size = 1.5, color = "darkred") + 
  # geom_ribbon( data = nx.tb, aes( x = Time, ymin = Upper, ymax = Lower ), alpha = 0.1 ) + 
  #geom_line( data = nx.g[ which(nx.g$type == "mcmc"), ], aes( x= time, y = e ), size = 1.5, color = "darkred", alpha = 0.5, linetype = "dashed") +
  # geom_ribbon( data = nx.g[ which(nx.g$type == "mcmc"), ], aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  # 
  #geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = 0.0001, ymax = 10000),
  #           fill = "gray", alpha = 0.2, inherit.aes = FALSE) +  

  coord_cartesian( ylim = c( 10^(-0.5), 100), xlim = c(1995, 2016) ) +
  scale_y_log10( breaks = scales::trans_breaks( "log10", function(x) 10^x),
                 labels = scales::trans_format( "log10", scales::math_format(10^.x)) ) +
  annotation_logticks() + 
  ylab( "Population size" ) + 
  scale_x_continuous( breaks = seq( 1995, 2017, by = 2)) + 
  theme( panel.grid.major  = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.title.x = element_blank(),
         axis.text.x = element_blank() ) +
  geom_vline( xintercept = 2006.3, linetype = "dashed", color = "gray") + 
  geom_vline( xintercept = 2009, linetype = "dashed", color = "gray")

# growth rate
h3 <- 
  h5n1.g %>%
  filter( type == "rate" ) %>% 
  select( lb, e, ub, time ) %>% 
  ggplot() + theme_bw() + 
  geom_line( aes( x= time, y = e ), size = 1.5) +
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  #geom_line( data = nx.g[ which(nx.g$type == "rate"), ], aes( x= time, y = e ), size = 1.5, color = "darkred") +
  # geom_ribbon( data = nx.g[ which(nx.g$type == "rate"), ], aes( x = time, ymin = lb, ymax = ub), alpha = 0.1 ) +
  
  coord_cartesian( xlim = c(1995, 2016), 
                   ylim = c( -2.5, 3 )) +
  scale_x_continuous( breaks = seq( 1995, 2016, by = 5), 
                      labels = seq( 1995, 2016, by = 5) ) + 
  #geom_rect( data = recdf, aes( xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf),
  #           fill = "gray", alpha = 0.2, inherit.aes = FALSE) + 
  ylab( "Growth rate" ) + 
  geom_hline( yintercept = 0, linetype = "dashed" ) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
         axis.title.x = element_blank()) +
  geom_vline( xintercept = 2006.3, linetype = "dashed", color = "gray") + 
  geom_vline( xintercept = 2009, linetype = "dashed", color = "gray")



ggarrange( h0, ggarrange( h2, h3, ncol = 1, nrow = 2,  align = "v"), ncol = 2,  nrow = 1, 
           widths  = c(0.5,1 ))

ggsave( width = 2.5, height = 3, units = "in", filename = "allh5n1.pdf", useDingbats=FALSE)

