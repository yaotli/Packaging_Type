# sequence 
# raxml tree
# dNdS output


.aa_recon = function( folderdir )
{
  require( tidyverse )
  require( ggtree )
  require( phangorn )
  require( seqinr )
  require( treeio )

  dNdSfile = paste0( folderdir, grep( "dNdS", list.files( folderdir ), value = TRUE ) )
  seqfile  = paste0( folderdir, grep( "fasta$", list.files( folderdir ), value = TRUE ) )
  trefile  = paste0( folderdir, grep( "^raxml", list.files( folderdir ), value = TRUE ) )
  
  if(  folderdir %in% c( dNdSfile, seqfile, trefile ) ){ stop( 'lack dNdS or sequence or tree file...' ) }
  
  
  # dS ----------------------------------------------------------------------
  cat( 'parse dS tree' )
  
  textin = read_file( dNdSfile )
  ds     = str_match( textin, "dS tree\\: ([\\(\\):;0-9A-Za-z_\\-\\.,]+) dN tree" )[,2]
  
  write( ds,  gsub( "dNdS$", "dS\\.tre", dNdSfile )  )
  system( paste0( "open ", gsub( "dNdS$", "dS\\.tre", dNdSfile ) )  ) 
  
  m0 = readline( prompt = "export as dS_nex.tre (enter to continue)" ) # export nexus file
  if( !TRUE %in% grepl( "dS_nex\\.tre", list.files( folderdir ) ) )
  {
    stop( 'inaccurate tree naming' )
    
  }else{ ds_trefile0 = paste0( folderdir, grep( "dS_nex\\.tre", list.files( folderdir ), value = TRUE ) ) }
  
  ds0.tre = treeio::read.nexus( ds_trefile0 )
  seq     = phyDat( read.FASTA( seqfile ,type = "DNA" ) )
  
  # fix taxa name 
  if( NA %in% match( attributes( seq )$names, ds0.tre$tip.label ) )
  {
    ds0.tre$tip.label <- 
      attributes( seq )$names[ match( str_match( ds0.tre$tip.label, "^[A-Z0-9]+")[,1], str_match( attributes( seq )$names, "^[A-Z0-9]+")[,1] ) ]
    
    if( NA %in% match( attributes( seq )$names, ds0.tre$tip.label ) ){ stop( "sequence and tree do not match" ) }
  }
  treeio::write.tree( ds0.tre, file = ds_trefile0  )
  system( paste0( "open ", gsub( "dNdS$", "dS_nex\\.tre", dNdSfile ) )  ) 
  
  m1 = readline( prompt = "rooted the tre and export as dS_nex.r.tre (enter to continue)" ) # export rooted nexus file
  if( !TRUE %in% grepl( "dS_nex\\.r", list.files( folderdir ) ) )
  {
    stop( 'inaccurate tree naming' )
  }else
  { ds_r_trefile = paste0( folderdir, grep( "dS_nex\\.r", list.files( folderdir ), value = TRUE ) )
    }
  
  ds.tre = treeio::read.beast( ds_r_trefile )
  

  # aa reconstruction -------------------------------------------------------
  cat( 'fitting and reconstruction' )
  tre <- treeio::read.tree( trefile )
  
  View( fortify( tre )[, c(1,2,4)] )
  m2  =  as.numeric( readline( prompt = "choose a outgroup/node as root" ) )
  
  if( m2 > treeio::Ntip( tre ) )
  {
    tre = root( tre, node =  m2, resolve.root = TRUE )
  }else
  {
    tre = root( tre, outgroup = m2, resolve.root = TRUE )
    }
  
  if( !is.rooted( tre ) ){ stop( 'Sth. wrong in rooting' ) }
  
  if( NA %in% match( attributes( seq )$names, tre$tip.label ) )
  {
    tre$tip.label <- 
      attributes( seq )$names[ match( str_match( tre$tip.label, "^[A-Z0-9]+")[,1], str_match( attributes( seq )$names, "^[A-Z0-9]+")[,1] ) ]
    
    if( NA %in% match( attributes( seq )$names, tre$tip.label ) ){ stop( "sequence and tree do not match" ) }
  }
  
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
  tre.merge = merge_tree( ds.tre, tre.att )
  #write.beast( tre.merge, file = gsub( "\\.fasta", "_merged.tre", seqfile ) )
  
  # aa change 
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
  
  cat( 'export aa annotated tree' )
  treeio::write.beast( tree.fin, file = gsub( "\\.fasta", "_aa.tre", seqfile ) )
  
  #v20190114
}


# find a sister strain as root --------------------------------------------

.root_seq <- function( seqfile = file.choose(), treefile, originfas )
{
  
  getDes_ <- function( node, curr = NULL )
  {
    if( is.null(curr) ){ curr <- vector() }
    
    edgemax   <- tree.d[ c(1,2) ]
    daughters <- edgemax[which( edgemax[,1] == node ), 2]
    
    curr <- c(curr, daughters)
    nd   <- which( daughters >= length( which( tree.d$isTip )) )
    
    if( length(nd) > 0)
    {
      for(i in 1:length(nd) ){ curr <- getDes_( daughters[ nd[i] ], curr ) }
    }
    return(curr)
  }
  
  tiplist = str_match( fastaEx( seqfile )$id, "^[A-Z0-9]+" )[,1] 
  
  tre    = treeio::read.nexus( treefile )
  tre$tip.label = gsub( "'", "", tre$tip.label )
  tree.d = as.data.frame( fortify( tre ) )
  
  tiplist_adj = tre$tip.label[ match( tiplist, str_match( tre$tip.label, "^[A-Z0-9]+" )[,1] ) ]
  
  # mrca
  mrca_c = ggtree::MRCA( obj = tre,  
                         tip = tiplist_adj )
  if( mrca_c == ( length(tre$tip.label) + 1) ){ stop() }
  
  # get the most closely related seq
  mrca_s = setdiff( which( tree.d$parent == tree.d$parent[ mrca_c ] ), mrca_c )
  
  if( mrca_s > length(tre$tip.label) )
  {
    taxa_s = getDes_( node = mrca_s )[ getDes_( node = mrca_s ) < length(tre$tip.label) ]
    outgp  = taxa_s[ which.min( tree.d$x[ taxa_s ] ) ]
    
  }else
  {
    outgp = mrca_s  
  }
  
  m0 = readline( prompt = paste0( tree.d$label[ outgp ], " (Y/N)\n") ) # export nexus file
  
  if( grepl( "y", m0, ignore.case = TRUE ) )
  {
    
    s.i   = which( fastaEx( originfas )$id == tree.d$label[ outgp ] )
    id.i  = c( fastaEx( seqfile )$id, tree.d$label[ outgp ] )
    seq.i = c( fastaEx( seqfile )$seq, fastaEx( originfas )$seq[s.i] )
    
    write.fasta( sequences = seq.i, names = id.i, file.out = gsub( "\\/p", "/", seqfile )  )
    
  }else
  {
    stop()
  }
  
}




  