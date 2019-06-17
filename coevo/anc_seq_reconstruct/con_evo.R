source( "./function.coevo.R" )
require( phangorn )
require( stringr )

# reconstruction of  HA ---------------------------------------------------

# detect convergent evo. on HA 

n615_tre   <- "./gsgd/processed/tree/raxml_c2344_615/RAxML_bipartitions.pH5_615"
n615_tre_r <- "./gsgd/processed/tree/raxml_c2344_615/raxml_pH5_c2344_615.tre"
n615_seq   <- "./gsgd/processed/pH5_c2344_1547ts.fasta"

tre   = treeio::read.tree( n615_tre )
seq   = phyDat( read.FASTA( n615_seq ,type = "DNA" ) )
tre_r = treeio::read.beast( n615_tre_r )

# View( fortify( tre )[, c(1,2,4)] )
# rooting
tre = root( tre, node =  758, resolve.root = TRUE )

tre$tip.label = gsub( "'", "", tre$tip.label )

fit  <- pml( tre, seq )
fit  <- update( fit, k = 4 )
fit  <- optim.pml( fit, model = "GTR", optGamma = TRUE, optEdge = TRUE, control = pml.control(trace = 0), optRooted = TRUE)
anc  <- ancestral.pml( fit, type =  "ml" )

seq_inferred <- 
  lapply( anc, 
          function(x)
          { 
            y   = attr( anc, "levels" )[ apply( x[ attr( anc, "index" ), ], MARGIN = 1, which.max ) ]
            amb = apply( x[ attr( anc, "index" ), ], MARGIN = 1, function(x) TRUE %in% duplicated( x[ x == max(x) ]  ) )
            y[ amb ] = "-"
            return(y)
          })

tre.fit = fit$tree
tre.n   = data_frame( node = 1: treeio::Nnode2( tre.fit ), nt = sapply( seq_inferred, function(x) c2s(x) ) )
tre.att = treeio::treedata( phylo = tre.fit, data  = tre.n )

# merge the two trees
tre.merge = merge_tree( tre_r, tre.att )

tre.merge_d <- fortify( tre.merge )
aa.change <- 
  sapply( as.list( tre.merge_d$node ), 
          function( x )
          {
            mx <- sapply( list( tre.merge_d$nt[[ tre.merge_d$parent[x] ]], tre.merge_d$nt[[x]] ),
                          function(x) translate( s2c(x) ) )
            
            pos <- which( ! mx[,1] == mx[,2] ) 
            pos <- pos[  mx[,1][pos] != "X" &  mx[,2][pos] != "X"  ]
            return( paste( paste0( mx[,1][pos], pos, mx[,2][pos] ), collapse = " " ) )
          })
aa.change[ aa.change == "" ] = " "

tre.aa   = data_frame( node = 1:Nnode2( tre.merge@phylo ), aa.change = aa.change )
tree.fin = treeio::treedata( phylo = tre.merge@phylo, data  = tre.aa )
treeio::write.beast( tree.fin, file = "./anc_seq_reconstruct/pH5_c2344_615_aa.tre" )


# summary n615 tree -------------------------------------------------------

n615_aa_tre <- read.beast( "./anc_seq_reconstruct/pH5_c2344_615_aa.tre" )
n615_aa_d   <- fortify( n615_aa_tre )

n615_aa   <- unlist( str_split( n615_aa_d$aa.change[ which( !n615_aa_d$isTip ) ], " ") )
n615_aa   <- n615_aa[ -which( n615_aa == "" ) ]
n615_aadf <- as.data.frame( sort(table(n615_aa)) )

# posLargerThan2 <- str_extract( n615_aadf$n615_aa, "\\d+" )[ which( n615_aadf$Freq > 2 ) ]
# for( p in 1: length(posLargerThan2) )
# {
#   t = 
#     viewResi( trefile = "./gsgd/processed/tree/raxml_gsgd_5180/raxml_pH5_gsgd_5180_e1026.tre",
#               fasfile = "./gsgd/processed/pH5_8334_trim2.2_7245_rmd_rmd2.fasta",
#               pos     =  as.numeric( posLargerThan2[p] )  )
#   ggsave( plot = t, filename = paste0( "resi_", posLargerThan2[p], ".pdf"), device = "pdf", width = 4, height = 5, units = "in")
# }


# gross detection of NA con. evo ------------------------------------------

source( "./summ_trees_r.R" )

summ_trees <- summ_trees[ seq(2, 24, by = 2) ]

# summarize the aa replacement, change to N2 numbering, and calculate freq
na_rep     <- lapply( summ_trees, 
                      function(x)
                      {
                        y = 
                        x %>% 
                          filter( aa == 1 | aa == 4 ) %>%
                          select( aa_node, sero, sam, aa )
                        
                        y$from    = str_match( y$aa_node, "([A-Z])(\\d+)([A-Z])" )[,2]
                        y$to      = str_match( y$aa_node, "([A-Z])(\\d+)([A-Z])" )[,4]
                        y$aa_node = str_match( y$aa_node, "([A-Z])(\\d+)([A-Z])" )[,3]
                        
                        tem  = y$from[ which( y$aa == 4 ) ]
                        y$from[ which( y$aa == 4 ) ] = y$to[ which( y$aa == 4 ) ]
                        y$to[ which( y$aa == 4 ) ] = tem
                        
                        
                        y$perc = 0 
                        for( r in 1: dim( y )[1] )
                        {
                          file_dir = paste0( "./h5_", y$sero[r], "/NA_h5", y$sero[r], "_", y$sam[r], ".fasta" )
                          fas_seq  = fastaEx( file_dir )$seq
                          fas_mx   = do.call( rbind, lapply( fas_seq, translate ) )
                          
                          po    = as.numeric( y$aa_node[r] )
                          repTo = y$to[r]
                          y$perc[r] = length( which( fas_mx[,po] == repTo ) )/length( fas_mx[,po] )
                        }
                        
                        y$aa_node = 
                        sapply( Nx_num( align_file = "./aa_numbering/nx_paring_data_ref.fasta",
                                        ref_file   = "./aa_numbering/nx_aa_ref.fasta",
                                        data_pos   = y$aa_node, 
                                        input_na   = str_extract( y$sero, "\\d" ),
                                        out_na     = rep( 2, length( y$aa_node ) ) ),
                                function(x) x[2,2] )
                        
                        return(y)
                      } 
                      )


y = 
do.call(rbind, na_rep)   %>% 
  filter( perc > 0.4 )   %>% 
  filter( sero == "n2" ) %>% 
  select( aa_node )      %>% table() %>% sort()

point = 412
sapply( na_rep, 
        function(x)
        {
          if( point %in% x$aa_node ){ return( paste0( x$sero[1], x$sam[1] ) ) }else
          {
              return( "" )
            }
        }
        )


na_rep$h5n2_g4_n2











ls_sum_trees <- lapply( summ_trees,
                        function(x)
                        {
                          y =
                          x %>%
                            filter( protein == "na" ) %>%
                            filter( aa == 1 | aa == 4 )

                          z = unique( as.numeric( str_extract( y$aa_node, "\\d+" ) ) )

                          if( length(z) < 1 ){ r = NA }else
                          {
                            r = Nx_num( align_file = "./aa_numbering/nx_paring_data_ref.fasta",
                                        ref_file   = "./aa_numbering/nx_aa_ref.fasta",
                                        data_pos   = z,
                                        input_na   = rep( as.numeric( str_extract( unique(x$sero), "\\d" ) ), length(z) ),
                                        out_na     = rep( 2, length(z) ) )

                            r <-  sapply( r, function(x){ return(x[2,2]) } )

                          }
                          return( r )
                        })

df_sum_trees      = data.frame( sort( table( unlist( ls_sum_trees ) ) ), stringsAsFactors = FALSE )
df_sum_trees$Var1 = as.numeric( as.character( df_sum_trees$Var1 ) )

df_sum_trees$nx   = sapply( as.list( df_sum_trees$Var1 ),
                            function(x)
                            {
                              nx = c()
                              for( i in 1: length( ls_sum_trees ) )
                              {
                               if( x %in% ls_sum_trees[[i]] )
                               {
                                 nx = c( nx, str_extract( names(ls_sum_trees)[ i ], "g\\d_n\\d" ) )
                               }
                              }
                              return( paste0( nx, collapse = "|" ) )
                            })



Nx_num( align_file = "./aa_numbering/nx_paring_data_ref.fasta",
        ref_file   = "./aa_numbering/nx_aa_ref.fasta",
        data_pos   = 416,
        input_na   = 2, 
        out_na     = 8)


viewResi( trefile = "./h5_n2/dS/h5n2_g4_n2/pNA_h5n2_g4_aa.tre",
          fasfile = "./h5_n2/pNA_h5n2_g4.fasta",
          pos     =  412  )



summ_mx <- do.call( rbind, summ_trees )

nx_mx = 
summ_mx %>% 
  filter( protein == "na" ) %>%
  filter( aa != 0  ) %>%
  filter( aa != 3  ) 
nx_mx$d_root[ which( nx_mx$aa == 4) ] = 0
nx_mx$rep = str_extract( nx_mx$aa_node, "\\d+")
nx_mx$dup = "0" 
nx_mx$dup[ duplicated( nx_mx$rep ) ] = "1"

nx_mx$size = 1

for( k in 1:dim(nx_mx)[1] )
{
  file_dir = paste0( "./h5_", nx_mx$sero[k], "/dS_r/h5", nx_mx$sero[k], "_", nx_mx$sam[k], "_", nx_mx$sero[k], "/", "NA_h5", nx_mx$sero[k], "_", nx_mx$sam[k], ".fasta" )
  fas_seq  = fastaEx( file_dir )$seq
  fas_mx   = do.call( rbind, lapply( fas_seq, translate ) )
  
  po    = as.numeric( nx_mx$rep[k] )
  repTo = str_match( nx_mx$aa_node[k], "([A-Z])[0-9]+([A-Z])" )[ ifelse( nx_mx$aa[k]==1, 3, 2 ) ]
  
  nx_mx$size[k] = length( which( fas_mx[,po] == repTo ) )/length( fas_mx[,po])
  print(k)
}



ggplot( nx_mx ) + 
  geom_point( aes( x = d_root, y = as.numeric(rep), color = sero, size = size, shape = sam, alpha = dup) ) + 
  scale_x_continuous( limit = c(0, 0.12) ) +
  scale_alpha_manual( values = c( 1, 0.1 ) )





ggplot( nx_mx ) + 
  geom_text( aes( x = d_root, y = as.numeric(rep), color = sero, size = size, alpha = dup, label = aa_node) ) + 
  scale_x_continuous( limit = c(0, 0.12) ) +
  scale_alpha_manual( values = c( 1, 0.1 ) ) 





