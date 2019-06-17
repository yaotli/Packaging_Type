source( "./function.coevo.R" )
require( seqinr )
require( stringr )
require( ggtree )
require( treeio )
require( tidyverse )
require( ape )



summ_trees = c()

fd1 <- grep( "h5_n[1-9]", list.files( "." ), value = TRUE )
for( i in 1: length( fd1 ) )
{
  if( i == 5 ){ fd2 <- grep( "h5n[1-9]_g[01234]_[hn1-9]+", list.files( paste0( "./", fd1[i], "/dS_com/") ), value = TRUE ) }else
  {
    fd2 <- grep( "h5n[1-9]_g[1234]_[hn1-9]+", list.files( paste0( "./", fd1[i], "/dS_r/") ), value = TRUE )  
  }  
  
  for( j in 1: length(fd2) )
  {
    
    # parameters
    aa_node  = c() # is a replacement or node
    d_root   = c() # distance to root 
    l_branch = c()
    aa       = c() # aa replacement - 1: aa replacement 2: no replacement 3: no placement on the outgreoup 4: replacement on outgroup
    
    ntax    = c()
    protein = c() # serotype 
    sero    = c()
    sam     = c() # sample
    
    if( i == 5 ){ fd   = paste0( "./", fd1[i], "/dS_com/", fd2[j], "/" ) }else
    {
      fd   = paste0( "./", fd1[i], "/dS_r/", fd2[j], "/" )
    }
    
    f_aa = paste0( fd, grep( "aa\\.tre", list.files( fd ), value = TRUE ) )
    f_r  = paste0( fd, grep( "[0-9]\\.r\\.tre", list.files( fd ), value = TRUE ) )
    
    # combine trees
    t.a = read.beast( f_aa )
    t.d = read.beast( f_r )
    
    t.d@phylo$tip.label = t.a@phylo$tip.label[ match( str_match( t.d@phylo$tip.label, "^[A-Z0-9]+"  )[,1], 
                                                      str_match( t.a@phylo$tip.label, "^[A-Z0-9]+"  )[,1] )  ]
    
    tre.l = data_frame( node = 1:Nnode2( t.d@phylo ), lth = fortify(t.d)$branch.length )
    t.l   = treeio::treedata( phylo = t.d@phylo, data  = tre.l )
    
    cbtree   <- merge_tree( t.a, t.l ) 
    tre_data <- fortify( cbtree )
    
    # walk along the branches
    ntax.i <- length( cbtree@phylo$tip.label )
    root.i <- ntax.i + 2
    out.i  <- min( which( tre_data$parent == root.i-1 ) )
    
    for( k in (root.i+1): dim(tre_data)[1] )
    {
      # node ancestral to 2 identical seq would not be considered as an internal node
      if( ( tre_data$lth[ which( tre_data$parent == k ) ][1] == tre_data$lth[ which( tre_data$parent == k ) ][2] ) |
          ( tre_data$lth[ which( tre_data$parent == k ) ][1] == 0.000001 ) ){ next() }
      
      # no replacment 
      if( tre_data$aa.change[k] == " " )
      {
        aa_node  = c( aa_node, 0 )
        aa       = c( aa, 0 )
        d_root   = c( d_root, ape::dist.nodes( cbtree@phylo )[ root.i, k ] )
        l_branch = c( l_branch, tre_data$branch.length[k] )
        
      # with replacement
      }else
      {
        aa.i = unlist( str_split( tre_data$aa.change[k], " " ) )
        n.aa = length( aa.i )
        
        aa_node  = c( aa_node, 0, aa.i )
        aa       = c( aa, 0, rep( 1, n.aa ) )
        d_root   = c( d_root, rep( ape::dist.nodes( cbtree@phylo )[ root.i, k ], (n.aa+1) ) ) 
        l_branch = c( l_branch, rep( tre_data$branch.length[k], (n.aa+1) ) )
      }
    }
    
    
    # outgroup / root 
    if( tre_data$aa.change[ out.i ] == " " )
    {
      aa_node  = c( aa_node, 0 )
      aa       = c( aa, 3 )
      d_root   = c( d_root, ape::dist.nodes( cbtree@phylo )[ root.i, out.i ] )
      l_branch = c( l_branch, tre_data$branch.length[k] )
      
    }else
    {
      aa.i = unlist( str_split( tre_data$aa.change[ out.i ], " " ) )
      n.aa = length( aa.i )
      
      aa_node  = c( aa_node, 0, aa.i )
      aa       = c( aa, 3, rep( 4, n.aa ) )
      d_root   = c( d_root, rep( ape::dist.nodes( cbtree@phylo )[ root.i, out.i ], (n.aa+1) ) ) 
      l_branch = c( l_branch, rep( tre_data$branch.length[k], (n.aa+1) ) )
    }
  
    n.row     = length( aa_node )
    
    protein.i = ifelse( str_match( fd2[j], "h5(n[0-9])_(g[0-9])_([hn0-9]+)" )[,4] =="h5", "ha", "na" )
    sero.i    = str_match( fd2[j], "h5(n[0-9])_(g[0-9])_([hn0-9]+)" )[,2]
    sam.i     = str_match( fd2[j], "h5(n[0-9])_(g[0-9])_([hn0-9]+)" )[,3]
    
    protein = rep( protein.i, n.row )
    sero    = rep( sero.i, n.row )
    sam     = rep( sam.i, n.row ) 
    ntax    = rep( ntax.i, n.row )
    
    summ_trees[[ length( summ_trees ) +1 ]] = data.frame( aa_node  = aa_node,
                                                          aa       = aa,
                                                          d_root   = d_root,
                                                          l_branch = l_branch, 
                                                          ntax     = ntax, 
                                                          protein  = protein, 
                                                          sero     = sero, 
                                                          sam      = sam, stringsAsFactors = FALSE )
    names( summ_trees )[ length( summ_trees ) ] = fd2[j]
    cat( paste0( fd2[j], "\n" ) )
  }
}


# count -------------------------------------------------------------------

# treedf_count = data.frame( samples = sub( "_[hn1-9]{2}$", "",  names( summ_trees ) ),
#                            count   = as.vector( sapply( summ_trees, function(x) x %>% filter( aa == 1 ) %>% select( aa ) %>% sum()  ) ),
#                            protein = ifelse( startsWith( str_match( names( summ_trees ), "_([hn1-9]{2})$")[,2], "h" ), "ha", "na" ),
#                            root    = as.vector( sapply( summ_trees, function(x) x %>% filter( aa == 4 ) %>% select( aa ) %>% nrow  ) ),
#                            stringsAsFactors = FALSE )
# 
# ggplot( treedf_count ) +
#   geom_bar( aes(x = samples, y = count, fill = protein ), color = "black", stat = "identity", position = position_dodge() ) +
#   scale_x_discrete( label = gsub( "h5|g", "", unique(treedf_count$samples) ) ) + 
#   ylab( "aa replacement" ) + xlab("groups") +
#   theme( legend.position = c( 0.1 , 0.9),
#          axis.title = element_text(size = 20) )
# 
# treedf_count2 = spread( treedf_count[,c(1,2,3)], key = protein, value = count )
# ggplot( treedf_count2 ) + geom_point( aes( x = ha, y = na)  ) + geom_abline( slope = 1, color = "gray", linetype = "dashed") +
#   annotate( geom = "text", x = 30, y = 200, label = "R^2 = 0.98", size = 8 ) + 
#   theme( axis.title = element_text(size = 20) )


#replacement on the root

# ggplot( treedf_count ) +
#   geom_bar( aes(x = samples, y = root, fill = protein ), color = "black", stat = "identity", position = position_dodge() ) +
#   scale_x_discrete( label = gsub( "h5|g", "", unique(treedf_count$samples) ) ) + ylab( "replacement on the RCB" ) + xlab("") +
#   theme( legend.position = c( 0.9, 0.8 ) )




# dS mean -------------------------------------------------------------------

# treedf_dS <- lapply( summ_trees,
#                     function(x)
#                     {
#                       x %>%
#                         filter( aa == 1 ) %>%
#                         select( d_root, sero, sam, protein )
#                     } )
# 
# treedf_dS <- do.call( "rbind", treedf_dS )
# treedf_dS <- unite( treedf_dS, col = group, ... = c(sero, sam))

# ggplot( treedf_dS ) +
#   geom_jitter( aes( x = group, y = d_root, color = protein ), width = 0.1, alpha = 0.5 ) + 
#   theme_bw() +
#   theme( legend.position = "none" ) + xlab( "" ) + ylab( "dS of nodes" ) + coord_flip() 
  


# convergent evo -------------------------------------------------------------









