source("functions.R")

require(ggplot2)
require(ggpubr)
require(skygrowth)
require(ape)
require(tidyverse)


c232_D_tre <- "./eco/results/201807_grid_eco/20180702_grid_eco_232_h5_D-anno.nwk"
c232_W_tre <- "./eco/results/201807_grid_eco/20180702_grid_eco_232_h5_W-anno.nwk"
c234_D_tre <- "./eco/results/201807_grid_eco/20180706_grid_eco_234_h5_D-anno.nwk"
c234_W_tre <- "./eco/results/201807_grid_eco/20180706_grid_eco_234_h5_W-anno.nwk"

t_c232_D <- 2011.989
t_c232_W <- 2011.948
t_c234_D <- 2011.953
t_c234_W <- 2009.496

# skygrowth 

ls.nwk  <- c( c232_D_tre, c232_W_tre, c234_D_tre, c234_W_tre)
ls.t    <- c( t_c232_D, t_c232_W, t_c234_D, t_c234_W )
ls.name <- c( "232D", "232W", "234D", "234W")

tb = skygrowth_df( ls.nwk, ls.t, name = ls.name, n.res = 50) 


skygrowth_df <- function( ls.nwk, ls.t,  n.res = 40, name = NA )
{
  require( skygrowth )
  require( ape )
  
  if( length(ls.nwk) != length(ls.t) ){ stop() }
  if( length(name) != length(ls.nwk) ){ label = as.character( seq( 1, length(ls.nwk) ) )   }
      
  ls.fit  <- list()
  ls.mcmc <- list()
  ls.rate <- list()
  
  for( i in 1: length( ls.nwk ) )
  {
    nwk <- read.tree( ls.nwk[i] )
    
    fit     = skygrowth.map( nwk, res = n.res, tau0 = 0.1 )
    mcmcfit = skygrowth.mcmc( nwk, res = n.res, tau0 = 0.1 )
    
    ls.fit[[ i ]]  <- data.frame( fit$ne_ci, time = fit$time + ls.t[i], label = name[i], type = "fit", stringsAsFactors = FALSE)
    ls.mcmc[[ i ]] <- data.frame( mcmcfit$ne_ci, time = mcmcfit$time + ls.t[i], label = name[i], type = "mcmc", stringsAsFactors = FALSE)
    ls.rate[[ i ]] <- data.frame( mcmcfit$growthrate_ci, time = mcmcfit$time + ls.t[i], label = name[i], type = "rate", stringsAsFactors = FALSE)
    
  }
  
  df.fit  <- do.call( rbind, ls.fit )
  df.mcmc <- do.call( rbind, ls.mcmc )
  df.rate <- do.call( rbind, ls.rate )
  
  colnames(df.fit)[ c(1,2,3) ]  <- c( "lb", "e", "ub" )
  colnames(df.mcmc)[ c(1,2,3) ] <- c( "lb", "e", "ub" )
  colnames(df.rate)[ c(1,2,3) ] <- c( "lb", "e", "ub" )
  
  return( do.call( rbind, list( df.fit, df.mcmc, df.rate ) )  )
}


tb[, 7] <- ifelse( startsWith( as.character(tb[, 5]) , prefix = "232" ),  "232", "234" )
colnames(tb)[7] = "Clade"

tb %>%
  filter( type == "rate" ) %>% 
  select( lb, e, ub, time, label, Clade ) %>%
  ggplot() + theme_bw() +
  facet_wrap( ~ Clade) +
  geom_line( aes( x = time, y = e, color = label), size = 2) + 
  geom_ribbon( aes( x = time, ymin = lb, ymax = ub, fill = label), alpha = 0.1) + 
  coord_cartesian( ylim = c(-1.5, 20), xlim = c(2002.5, 2012) ) +
  scale_x_continuous( breaks = seq(2003, 2011, by = 2), labels =  seq(2003, 2011, by = 2) ) +
  xlab( "Year" ) + ylab( "log(Net)" ) + 
  scale_color_manual( values = pyCol( c( "green", "red", "blue" ) ) ) +
  theme( panel.grid.minor.y =  element_blank(), 
         strip.background = element_rect( fill = "white", color = "white"),
         panel.border = element_rect( color = "black", fill = NA, size = 1),
         legend.position = "none" ) 










