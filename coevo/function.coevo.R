### fastaEx --------------------------------

fastaEx <- function(filedir = file.choose())
{
  require(seqinr)
  
  file     <- seqinr::read.fasta(filedir)
  file_seq <- getSequence(file)
  file_id  <- attributes(file)$names
  
  return( list(seq = file_seq, 
               id  = file_id ) )
  #v201706
}

### keepLongSeq --------------------------------

keepLongSeq <- function(seq_0, 
                        id_0, 
                        showRemain = FALSE)
{
  require(seqinr)
  
  if( length( which( duplicated(id_0) == TRUE) ) > 0  )
  {
    toberemove <- c()
    dup        <- which( duplicated(id_0) == TRUE)
    
    for( k in 1: length(dup) )
    {
      id_dup_k <- which(id_0 %in% id_0[ dup[k] ] == TRUE)
      SeqL     <- which.max(
        
        sapply(seq_0[id_dup_k], function(x)
        {
          
          y = grep( pattern     = "a|t|c|g",
                    x           = x,
                    ignore.case = TRUE, 
                    value       = TRUE)
          
          l = length( y )
          
          return(l)
          
        })
      )
      
      toberemove <- c(toberemove, id_dup_k[-SeqL])  
    }
    
    toberemove <- unique( sort(toberemove) )
    remain     <- seq(1, length(seq_0))[-toberemove]
    
    if (showRemain == TRUE)
    {
      return(remain)
      
    }else
    {
      return( list(seq = seq_0[remain], 
                   id  = id_0[remain]) ) 
    }
    
  }else
  {
    print("No identical ID here")
  }
  
  #v20171124
}


### idInfo --------------------------------

idInfo <- function( rawid, 
                    datasource = "n",
                    g.csv      = "" )
{
  # format: 
  # N:  >{accession}_{strain}_{serotype}_|{country}|_{year}-{month}-{day}
  # G:  Isolate name Type Collection date Isolate ID
  # Both need replace the blank with underline 
  
  require(seqinr)
  require(stringr)
  
  # g
  a.string.g <- "EPI_ISL_([0-9]+)"
  s.string.g <- "_A_/_(H[0-9]{1,2}N[0-9xX]{1,2})_"
  y.string.g <- "_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)" 
  r.string.g <- "^([0-9A-Za-z\\.\\/\\-_\\(\\)\\?]+)_A_/_"
  
  # n 
  a.string.n   <- "^[A-Za-z0-9]+"
  s.g.string.n <- "(H[0-9]{1,2}[N0-9xX]{0,2})_\\|([a-zA-Z_\\']+)\\|"
  y.string.n   <- "_[0-9]{4}-[0-9]{2}-[0-9]{2}|_[0-9]{4}-[0-9]{2}-|_[0-9]{4}--|_--"
  n.Nstring.n  <- "^[A-Za-z0-9]+_|(H[0-9]{1,2}[N0-9xX]{0,2})_\\|([a-zA-Z_\\']+)\\||_([0-9-]+)$"
  
  
  if( datasource == "g")
  {
    id.a <- gsub( "_ISL_", "", str_match( rawid, a.string.g )[, 1] )
    id.s <- str_match( rawid, s.string.g )[,2] 
    id.s[ which( id.s == "H5"| id.s == "H5Nx"| id.s == "H5NX")  ] = "H5N0"
    id.y <- str_match( rawid, y.string.g )
    id.y <- gsub( "^_", "", x = id.y)[,1]
    
    id.n <- str_match( rawid, r.string.g )[,2]
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- gsub( "_A/", "A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ] )
    
    id.n <- gsub( "\\?|\\(|\\)|\\[|\\]|\\.|:|-|/", "_", id.n )
    id.n <- gsub( "\\'|\\?|>", "", id.n )
    id.n <- gsub("A_", "", id.n)
    id.n <- gsub( "[_]+", "_", id.n )
    id.n <- gsub( "_$", "", id.n )
    id.n <- gsub( "^_", "", id.n )
    
    g <- gsub( " ", "_", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Location )
    g <- gsub( "_$", "",  str_match( g, "([A-Za-z_]+)_/_([A-Za-z_]+)" )[,3] ) 
    
    g[ which( is.na(g) == TRUE ) ] = "Unknown"
    g[ which(g == "Russian_Federation") ] = "Russia"
    g[ which(g == "United_States") ] = "USA"
    g[ which(g == "Korea") ] = "South_Korea"
    
    id.g <- g[ match( id.a, gsub("_ISL_", "", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Isolate_Id ) ) ]
    
    
  }else
  {
    id.a <- str_match( rawid, a.string.n)[,1]
    id.s <- str_match( rawid, s.g.string.n)[,2]
    id.s[ which( id.s == "H5"| id.s == "H5Nx"| id.s == "H5NX")  ] = "H5N0"
    id.g <- str_match( rawid, s.g.string.n)[,3]
    id.g[ which( id.g == "Viet_Nam") ] = "Vietnam"
    id.g[ which( id.g == "Cote_d'Ivoire") ] = "Cote_dIvoire"
    
    id.y <- str_match( string = rawid, y.string.n)
    id.y <- gsub( "_--", "1900-01-01", id.y)
    id.y <- gsub( "^_", "", id.y)
    
    id.n <- gsub( n.Nstring.n, "", rawid)
    
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- 
      paste0("A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ])
    
    id.n <- gsub("\\?|\\(|\\)|\\[|\\]|\\.|:|-|/|__", "_", id.n)
    id.n <- gsub("\\'|\\?|>", "", id.n)
    id.n <- gsub("A_", "", id.n)
    id.n <- gsub("[_]+", "_", id.n)
    id.n <- gsub("_$", "", id.n)
    id.n <- gsub("^_", "", id.n)
    
  }
  
  infolist = list(id.a, id.s, id.g, id.y, id.n)
  
  e = 
    which(
      sapply( infolist, 
              function(x)
              {
                TRUE %in% is.na(x)
                
              })  == TRUE )
  
  
  print( paste("ERROR in ", c("ac", "sero", "geo", "year", "name")[e] )  )
  
  return(infolist)
  
  #v20181002e
}



### strainSelect --------------------------------

strainSelect <- function( infolist )
{
  
  infolist.n <- infolist[[ length(infolist) - 1 ]]
  infolist.q <- infolist[[ length(infolist) ]]
  infolist.y <- infolist[[ length(infolist) - 2 ]]
  infolist.a <- infolist[[1]]
  
  toberemove <- c() 
  dup        <- which( duplicated( infolist.n ) )
  
  for(i in 1: length(dup) )
  {
    id_dup_ii <- which( infolist.n %in% infolist.n[ dup[i] ] == TRUE )
    lth_ii    <- sapply( infolist.q[id_dup_ii],
                         
                         function(x)
                         {
                           
                           z = grep( pattern = "a|t|c|g", 
                                     x = x, 
                                     ignore.case = TRUE, value = TRUE )
                           
                           l = length( z )
                           
                           return(l)
                           
                         } )
    
    SeqL      <- which.max( lth_ii ) 
    
    if ( length(  which( lth_ii == max(lth_ii) )  ) > 1 )
    {
      
      id_dup_jj  <- id_dup_ii[ which( lth_ii == max(lth_ii) ) ] 
      
      nchar_jj   <- nchar( gsub( pattern     = "[-\\(\\)A-Za-z]+", 
                                 replacement = "",
                                 x           = infolist.y[id_dup_jj] ) )
      
      id_dup_j   <- which.max( nchar_jj )
      SeqL       <- which( lth_ii == max(lth_ii) )[id_dup_j]
      
      
      if (  length( which( nchar_jj == max(nchar_jj) ) ) > 1  )
      {
        
        id_dup_kk <- id_dup_jj[ which(nchar_jj == max(nchar_jj) ) ]
        
        ac        <- infolist.a[id_dup_kk]
        ac.a      <- nchar( gsub( pattern = "[0-9]+", replacement = "", x = ac) )
        ac.d      <- as.numeric( gsub( pattern = "[a-zA-Z]+", replacement = "", x = ac ) )
        ac.df     <- data.frame( id_dup_kk, ac.a, ac.d )
        
        SeqL      <- which( id_dup_ii == ac.df[order( ac.df[,2], ac.df[,3] ),][1,1] )
        
      }
    }
    
    toberemove = c( toberemove, id_dup_ii[-SeqL] )
    
  }
  
  remain <- seq( 1, length(infolist.q) )[- toberemove]
  
  newlist = list()
  for(l in 1 : length(infolist) )
  {
    newlist[[l]] <- infolist[[l]][remain]
    
  }
  
  newlist[[ length(newlist) + 1 ]] <- ifelse( grepl( pattern = "--$|Month", x = newlist[[4]] ), 1, 0)
  
  if( TRUE %in% is.na( unlist(newlist) ) ){ print("ERROR") }
  
  return( newlist )
  
  #v20171108
}

### seqDate --------------------------------

seqDate <- function( rawdata )
{
  require(stringr)
  require(lubridate)
  
  # gisaid
  
  rawdata.1 <- gsub( "_\\(Day_unknown\\)", "-15", 
                     gsub( "_\\(Month_and_day_unknown\\)", "-07-01", rawdata ) )
  
  # ncbi
  
  rawdata.2 <- gsub( "-$", "-15", 
                     gsub( "--$", "-07-01", rawdata.1) )
  
  # parse into numeric
  
  d = "([0-9]{4})-([0-9]{1,2})-([0-9]{2})"
  
  yr   <- as.numeric( str_match(rawdata.2, d)[,2] )
  yr.0 <- paste0(yr, "-01-01")
  yr.e <- paste0(yr, "-12-31")
  
  daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d", tz = "GMT"),
                                         strptime( yr.0, "%Y-%m-%d", tz = "GMT"), 
                                         units = "days") 
  )/yday(yr.e)
  
  # bug?
  if ( TRUE %in% is.na(daydifference) )
  {
    
    rawdata.2[ which(is.na(daydifference)) ] <- 
      sub("01$", "02", rawdata.2[ which(is.na(daydifference)) ] )
    
    
    daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d", tz = "GMT"),
                                           strptime( yr.0, "%Y-%m-%d", tz = "GMT"), 
                                           units = "days") 
    )/yday(yr.e)
  }
  
  yr.daydifference <- yr + daydifference
  yr.daydifference <- format( round( yr.daydifference, 3 ), nsmall = 3)
  
  return(yr.daydifference)
  
  #v20180409
}


### seqSelect --------------------------------

seqSelect <- function( minlth  = 1000, 
                       maxamb  = 1,
                       rmdup   = TRUE,
                       seqlist )
{
  df   <- data.frame( a = seqlist[[1]], 
                      s = seqlist[[2]],
                      y = seqlist[[4]],
                      i = seqlist[[7]], stringsAsFactors = FALSE) 
  
  df.s <- df[ order(df$i, df$y, df$a, df$s), ]
  
  idx  <- as.numeric( rownames(df.s) )
  q    <- seqlist[[6]][ idx ]
  
  
  # length and ambiguous nucleotide
  lth_amb <- which( sapply( q, 
                            
                            function(x)
                            {
                              ATCG  <-  c("a", "t", "c", "g")
                              
                              x.s   <- gsub( "-|~", "", c2s( x ) )
                              x.l   <- length( s2c(x.s) )
                              
                              x.c.s <- grep( "a|t|c|g", s2c(x.s) )[1]
                              x.c.e <- grep( "a|t|c|g", s2c(x.s) )[ length( grep( "a|t|c|g", s2c(x.s) ) ) ]
                              
                              x.a   <- length( which(! s2c(x.s)[x.c.s: x.c.e] %in% ATCG ) )
                              
                              return( x.l < minlth | x.a > (maxamb/100)*x.l )
                              
                            } ) ) 
  # duplicated sequence
  if (rmdup)
  {
    dup     <- which( duplicated( sapply( q,  
                                          function(x)
                                          {
                                            x.s <- gsub( "~|-", "", c2s(x) )
                                            return(x.s)
                                          }) 
    ) )
    
  }else{
    dup = c()
  }
  
  if ( ( length(dup) + length(lth_amb) ) > 0 )
  {
    remain  <- seq(1, length( seqlist[[6]] ) )[ - unique( sort( c(dup, lth_amb) )) ]
    
  }else
  {
    remain  <- seq(1, length( seqlist[[6]] ) )  
  }
  
  
  newlist = list()
  for(l in 1 : length(seqlist) )
  {
    newlist[[l]] <- seqlist[[l]][idx][remain]
  }
  
  print( paste0("n = ", length( newlist[[1]] )) )
  
  return(newlist)
  
  #v20170921b
}




### gapFill --------------------------------

gapFill <- function( file.dist = file.choose(), 
                     s.start   = 1000, 
                     s.end     = 1110 )
{
  library(seqinr)
  
  file       <- read.fasta( file.dist )
  seq.name0  <- attributes( file )$names
  seq0       <- getSequence( file )
  
  seq <- 
    lapply( seq0, 
            function(x)
            {
              x1 <- grep( "-", x[ s.start: s.end], invert = TRUE, value = TRUE ) 
              x2 <- c( x1[c(1,2,3)], rep( "-", ( s.end - s.start + 1 - length(x1) ) ), x1[4: length(x1) ] )
              x[ s.start: s.end] <- x2
              return(x)
            }
    )
  
  write.fasta( sequences = seq, file.out = gsub( "\\.fasta|\\.fas", "_gapfilled.fasta", file.dist ),
               names     = seq.name0 )
  
  #v20180920
}

### tagExtra --------------------------------

tagExtra <- function( filedir = file.choose() )
{
  require(stringr)
  
  anno.tre <- read.csv( filedir, stringsAsFactors = FALSE)
  taxa.s   <- grep( x = anno.tre[,1], pattern = "taxlabels" ) + 1
  ntax     <- as.numeric( str_match( grep( x       = anno.tre[,1], 
                                           pattern = "ntax", value = TRUE) , "(ntax=)([0-9]+)")[,3] 
  )
  
  taxa.e <- taxa.s + ntax - 1
  
  GsGDlike_name <- str_match( string = anno.tre[, 1][taxa.s: taxa.e], 
                              pattern = "\'([0-9A-Za-z_\\|.]+)\'" )[,2]
  
  GsGDlike_tag  <- str_match( string = anno.tre[, 1][taxa.s: taxa.e], 
                              pattern = "color=#([a-z0-9]{6})")[, 2]
  
  return( df = data.frame(id  = GsGDlike_name, 
                          tag = GsGDlike_tag, stringsAsFactors = FALSE))
  
  #v20180127
}

### leafEx --------------------------------

leafEx <- function( filedir   = file.choose(), 
                    leaflist  = c(), 
                    consenseq = FALSE, 
                    seq.out  = NA )
{
  require( seqinr )
  
  leaflist <- unlist( strsplit( gsub( " ", "", leaflist ), "\n" ) )
  
  fas   <- seqinr::read.fasta( filedir )
  fas.s <- getSequence( fas )
  fas.n <- attributes( fas )$names
  
  m <- match( leaflist, fas.n )
  
  if( TRUE %in% is.na(m) ){ stop() }else
  {
    
    if( !is.na(seq.out) )
    {
      write.fasta( fas.s[ m ], fas.n[ m ], file.out = seq.out )  
      
    }else
    {
      write.fasta( fas.s[ m ], fas.n[ m ], gsub( pattern = ".fasta", 
                                                 replacement = paste0( "_", length(m), ".fasta"), 
                                                 x = filedir ) )    
      
      }
    
    
  }
  
  mx <- do.call( rbind, fas.s[ m ] )
  
  if( consenseq == TRUE )
  {
    Conseq = 
      apply( mx, 2, 
             function(x){
               
               if ( length( grep( "a|t|c|g", x, value = TRUE, ignore.case = TRUE) )  == 0 ){ most = "-" }else
               {
                 most <- as.data.frame( table( grep( "a|t|c|g", x, value = TRUE, ignore.case = TRUE) ), 
                                        stringsAsFactors = FALSE)[1,1]
               }
               
               return( most )
               
             })
    
    write.fasta( list(Conseq), "Conseq", gsub( pattern = ".fasta", 
                                               replacement = paste0( "_con", length(m), ".fasta"), 
                                               x = filedir ) )
  }
  
  #v20190121
}

### rmDup --------------------------------

rmDup <- function( fasfile     = file.choose(), 
                   year        = c(1000,3000),
                   geo         = c(),
                   sero        = "",
                   rmdup       = TRUE, 
                   sero.invert = FALSE )
{
  require(seqinr)
  require(stringr)
  
  readin <- read.fasta( fasfile )
  seq    <- getSequence( readin )
  id     <- attributes( readin )$names
  
  if( rmdup )
  {
    
    # order: time ( data completeness ), accession number ( data source )
    
    id.y  <- as.numeric( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2] )
    id.a  <- str_match( id, "EPI[0-9]+|[A-Z]{1,2}[0-9]{5,6}" )[,1]
    
    id.d  <- 
      ifelse( ( endsWith( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], ".496" )|
                  endsWith( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], ".497" ) )
              , 1, 0)
    
    id.a.c <- nchar( gsub("[0-9]+", "", id.a) )
    id.a.d <- as.numeric( gsub("[A-Za-z]+", "", id.a) )
    
    
    df <- data.frame( id.y, id.d, id.a.c, id.a.d )  
    df <- df[ order( df[,2], df[,1], df[,3], df[,4] ), ]
    
    seq <- seq[ as.numeric( rownames(df) ) ]
    id  <- id[ as.numeric( rownames(df) ) ]
    
    
    dup <- which( duplicated( sapply( seq, 
                                      function(x)
                                      {
                                        x.s <- gsub( "~|-", "", c2s(x) )
                                        return(x.s)
                                      } ) 
    ) )
    
  }else
  {
    dup = NA
  } 
  
  if( length(dup) > 1 )
  {
    remain = seq( 1, length(seq) )[ - sort( unique(dup) ) ]
    
  }else
  {
    remain = seq( 1, length(seq) )
    
  }
  
  # year
  y = which( (str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2] > year[1] & 
                str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2] < year[2]) )
  
  # geo
  geo.p <- paste0( geo, collapse = "|" )
  g     <- grep( geo.p, str_match( id, "\\|([A-Za-z_]+)\\|" )[,2] )
  
  # sero
  
  s     <- grep( sero, str_match( id, "\\|_(H[0-9]{1,2}N[0-9]{1,2})_" )[,2], invert = sero.invert )
  
  if ( TRUE %in% is.na( 
    c(str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], 
      str_match( id, "\\|([A-Za-z_]+)\\|" )[,2], 
      str_match( id, "\\|_(H[0-9]{1,2}N[0-9]{1,2})_" )[,2]) ) )
  {
    stop( 
      c("Year", "Geo", "Serotype")[ ceiling( 
        which( is.na( 
          c(str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], 
            str_match( id, "\\|([A-Za-z_]+)\\|" )[,2], 
            str_match( id, "\\|_(H[0-9]{1,2}N[0-9]{1,2})_" )[,2]) ))/3) ] 
    )
  }  
  
  remain <- sort( Reduce( intersect, list( remain, y, g, s ) ) )
  
  write.fasta( seq[ remain ], 
               id[remain], 
               file.out = sub(".fasta", "_rmd.fasta", fasfile) )
  
  print( length( remain ) )
  
  #v20190106
}

### rmdup_plus --------------------------------

rmdup_plus <- function( fasdir = file.choose() )
{
  require(seqinr)
  require(stringr)
  
  fas = read.fasta( fasdir )
  
  #readin
  fas.s0 <- getSequence(fas)
  fas.i0 <- attributes(fas)$names
  
  # sort1
  o.t <- 
    sort( as.numeric( str_match( fas.i0, "_([0-9]{4}.[0-9]{3})$" )[,2] ), index.return = TRUE )$ix
  
  if( TRUE %in% is.na(as.numeric( str_match( fas.i0, "_([0-9]{4}.[0-9]{1,3})$" )[,2] )) ){ stop() }
  
  fas.s1 <- fas.s0[ o.t ] 
  fas.i1 <- fas.i0[ o.t ]
  
  # sort2
  o.lth <- 
    sort( sapply(fas.s1, 
                 function(x)
                 {
                   length( grep( pattern = "a|t|c|g", x = x, ignore.case = TRUE, value = TRUE ) )
                 }
    ), index.return = TRUE )$ix
  
  fas.s1 <- fas.s1[ o.lth ] 
  fas.i1 <- fas.i1[ o.lth ]
  
  
  
  m   <- matrix( unlist(fas.s1), ncol = length( fas.s1[[1]] ), byrow = TRUE )
  
  todel <- c()
  for(i in 1: ( length(fas.s1) -2 ) )
  {
    m_i <- grep( pattern = "a|t|c|g", x = m[i,], ignore.case = TRUE)
    
    if( TRUE %in% grepl( c2s( m[ i, m_i] ), apply( m[ (i+1) : nrow(m), m_i], 1, c2s) ) ){ todel[ length(todel) + 1 ] <-  i }
    
    print( i )
  }
  
  remain <- seq(1, length(fas.s1) )[- todel ]
  
  print( paste0( "Done:", length(remain) ) )
  
  write.fasta( sequences = fas.s1[remain], names = fas.i1[remain], file.out = gsub( ".fasta", "_rmd2.fasta", fasdir) )
  
  #v20190106
}

### getDescendant --------------------------------

getDes <- function( tre.d = trefile, node, curr = NULL )
{
  if( is.null(curr) ){ curr <- vector() }
  
  edgemax   <- tre.d[ c(2,1) ]
  daughters <- edgemax[which( edgemax[,1] == node ), 2]
  
  curr <- c(curr, daughters)
  nd   <- which( daughters >= length( which( tre.d$isTip )) )
  
  if( length(nd) > 0)
  {
    for(i in 1:length(nd) ){ curr <- getDes( tre.d = tre.d, daughters[ nd[i] ], curr ) }
  }
  return(curr)
  
  #v20180628    -> modified from unknown source  
}


### cladesampling --------------------------------   


cladeSampling <- function( trefile      = file.choose(),
                           fasfile      = file.choose(),
                           listinput    = list(),
                           seed         = 666,
                           grid         = 1,
                           minBranchlth = TRUE, 
                           showTree     = FALSE, 
                           saveFasta    = FALSE,  
                           suppList     = FALSE,
                           list.x       = c("id", "y", "info1", "info2") )
{
  require( ape )
  require( seqinr )
  require( ggtree )
  require( stringr )
  
  # getdescendants
  getDes <- function( node, curr = NULL )
  {
    if( is.null(curr) ){ curr <- vector() }
    
    edgemax   <- tre.d[ c(1,2) ]
    daughters <- edgemax[which( edgemax[,1] == node ), 2]
    
    curr <- c(curr, daughters)
    nd   <- which( daughters >= length( which( tre.d$isTip )) )
    
    if( length(nd) > 0)
    {
      for(i in 1:length(nd) ){ curr <- getDes( daughters[ nd[i] ], curr ) }
    }
    return(curr)
  }
  
  
  nex <- read.nexus( trefile )
  fas <- read.fasta( fasfile )
  seq <- getSequence( fas )
  id  <- attributes( fas )$names
  
  tre.d   <- as.data.frame( fortify( nex ) )
  N.tip   <- length( which( tre.d$isTip ) )
  N.node  <- nex$Nnode
  edgemax <- tre.d[ c(1,2) ]
  
  
  t.id <- gsub("'", "", tre.d$label)
  
  if( suppList )
  {
    
    tem.m  <- match(   str_match(t.id[ 1:  N.tip], "^[A-Za-z0-9]+")[,1]  , 
                       str_match(listinput[[ list.x[1] ]], "^[A-Za-z0-9]+")[,1]   )
    
    if( TRUE %in% is.na(tem.m) ){stop()}
    
    i.id.y <- listinput[[ list.x[2] ]][ tem.m ]
    i.id.g <- listinput[[ list.x[3] ]][ tem.m ]
    
    if( length( list.x ) >3 )
    {
      for( k in 1: ( length(list.x) -3 ) )
      {
        i.id.s <- listinput[[ list.x[ 3 + k ] ]][ tem.m ]
        i.id.g <- paste0( i.id.g, "-", i.id.s )
      }
    }
    
  }else
  {
    i.id.y <- as.numeric( str_match( t.id, "_([0-9]{4}\\.[0-9]+)$")[,2] )
    i.id.g <- str_match( t.id, "\\|([A-Za-z_]+)\\|")[,2]  
  }
  
  # 1st search for node with homogeneous descendants 
  inner.node <- seq( 1, dim( tre.d )[1])[ - seq(1, N.tip+1) ]
  c1.node    <- inner.node[ which( sapply( as.list(inner.node), 
                                           function(x)
                                           {
                                             all <- getDes(x)[ getDes(x) <= N.tip ]
                                             
                                             if( length( all[ !is.na( i.id.g[all] ) ] ) < 1 ){ return( TRUE ) }else
                                             {
                                               r   <- range( i.id.y[ all[ !is.na( i.id.g[all] ) ] ] )[2] - range( i.id.y[ all[ !is.na( i.id.g[all] ) ] ] )[1]
                                               g   <- unique( i.id.g[ all[ !is.na( i.id.g[all] ) ] ] )
                                               
                                               return( ( r <= grid ) & ( length( g ) == 1 ) )    
                                               
                                             } 
                                           } )) 
                            ]
  # reduce redndant nodes
  c2.node    <- c1.node[ which( sapply( as.list(c1.node), 
                                        function(x)
                                        {
                                          if( edgemax[,1][ which( edgemax[,2] == x) ] %in% c1.node )
                                          {
                                            return( FALSE )
                                            
                                          }else
                                          {
                                            return( TRUE )
                                          }
                                          
                                        } )
  )  
  ] 
  
  c2.node_des <- c( c2.node, unlist( sapply( as.list(c2.node), getDes) ) )
  c2.node_tip <- c2.node_des[ c2.node_des <= N.tip ]
  nogroup_tip <- seq(1, N.tip)[ - c2.node_tip ]
  nogroup_tip <- nogroup_tip[ !is.na( i.id.g[ nogroup_tip ] ) ]
  
  # sample within a group
  if( minBranchlth )
  {
    selected_tip <- sapply( as.list(c2.node), 
                            function(x)
                            {
                              alltip <- getDes(x)[ getDes(x) <= N.tip ]
                              
                              if( TRUE %in% is.na( i.id.g[ alltip ] )  )
                              {
                                alltip <- alltip[ - which( is.na( i.id.g[ alltip ] ) ) ]
                              }
                              
                              m      <- which.min( tre.d$x[ alltip ] )
                              
                              if( length( which( tre.d$x[ alltip ] == tre.d$x[ alltip ][m] ) )  > 1 )
                              {
                                set.seed( seed ) 
                                s = sample( which( tre.d$x[ alltip ] == tre.d$x[ alltip ][m] ), 1 )
                                
                              }else
                              {
                                s = m
                              }
                              return( alltip[s] )
                            })
    
  }else
  {
    selected_tip <- sapply( as.list(c2.node), 
                            function(x)
                            {
                              alltip <- getDes(x)[ getDes(x) <= N.tip ]
                              
                              if( TRUE %in% is.na( i.id.g[ alltip ] )  )
                              {
                                alltip <- alltip[ - which( is.na( i.id.g[ alltip ] ) ) ]
                              }
                              
                              set.seed( seed ) 
                              s = sample( alltip, 1)
                              return(s)
                              
                            } 
    )
  }  
  
  if( showTree )
  {
    # view the result
    tre.d[, ncol(tre.d) + 1 ] = "gray"
    colnames(tre.d)[ ncol(tre.d) ] = "colorr"
    tre.d$colorr[c2.node_des] = "red"
    
    tre.d[, ncol(tre.d) + 1 ] = NA
    colnames(tre.d)[ ncol(tre.d) ] = "shapee"
    tre.d$shapee[selected_tip] = 16
    
    g1 <- ggtree( nex ) %<+% tre.d + aes(color = I(colorr)) + geom_tiplab(size = 1)
    print( g1 + geom_tippoint(aes( shape = factor(shapee) ), size = 2) )
  }
  
  remain  <- c( nogroup_tip, selected_tip )
  
  f.match <- match( str_match( t.id[ sort(remain) ], "^[A-Za-z0-9]+")[,1], 
                    str_match( id, "^[A-Za-z0-9]+")[,1])
  
  seq.o   <- seq[ f.match ]
  id.o    <- id[ f.match ]
  
  if( saveFasta )
  {
    write.fasta( seq.o, id.o, 
                 file.out = gsub( ".fasta", "_s.fasta", fasfile) )
  }
  
  print( paste0("sampled n = ", length(remain), " from ", length(id) ) )
  print( table( floor( i.id.y[remain] ), i.id.g[remain] ) )
  
  
  #v20181116
}



### taxaInfo --------------------------------

taxaInfo <- function( file     = file.choose(), 
                      useTree  = FALSE, 
                      makecsv  = FALSE, 
                      root2tip = FALSE)
{
  # input: 
  # 1 colored .tre file
  # 2 .fas file with clean id 
  
  require(seqinr)
  require(stringr)
  require(ape)
  require(ggtree)
  
  if ( useTree )
  {
    anno.tre <- read.csv( file, stringsAsFactors = FALSE)
    nex      <- read.nexus( file )
    
    taxa.s   <- grep( "taxlabels", anno.tre[,1] ) + 1
    
    ntax     <- as.numeric( str_match( grep( "ntax", anno.tre[,1],  value = TRUE ), 
                                       "(ntax=)([0-9]+)" )[,3] )
    taxa.e   <- taxa.s + ntax - 1
    
    id  <- str_match( anno.tre[, 1][taxa.s: taxa.e], "\'([0-9A-Za-z_\\|.]+)\'" )[,2]
    tag <- str_match( string = anno.tre[, 1][taxa.s: taxa.e], 
                      pattern = "color=#([a-z0-9]{6})")[, 2]
    
    id.a <- str_match( id, "[A-Z]{1,2}[0-9]{5,6}|EPI[0-9]+" )[,1]
    id.g <- str_match( id, "\\|([A-Za-z_]+)\\|")[,2]
    id.s <- str_match( id, "\\|_(H[0-9]{1,2}N[0-9]{1,2})_")[,2]
    id.y <- as.numeric( str_match( id, "_([0-9]{4}.[0-9]{3})$")[,2] )  
    id.n <- gsub("^[A-Z]{1,2}[0-9]{5,6}_|^EPI[0-9]+_|_\\|[A-Za-z_]+\\|_|H[0-9]{1,2}N[0-9]{1,2}_[0-9]{4}.[0-9]{3}$", "", id)
    
    ls <- list( id.a, id.g, id.s, id.y, id.n, id, tag)
    df <- data.frame( id.a, id.g, id.s, id.y, id.n, id, tag )
    
    if( TRUE %in% is.na(unlist( ls[c(1:6)] ) ) ){ stop() }
    
    if( root2tip )
    {
      # distance max
      
      nexdata   <- fortify( nex )
      root.node <- length( nex$tip.label ) + 1
      root.dist <- dist.nodes( nex )[ root.node, 1: (root.node - 1) ]
      tre.id    <- gsub("'", "", nex$tip.label[ 1: root.node - 1])
      
      m <- match( tre.id, id )
      
      dist.df <- data.frame( name = id[m], geo = id.g[m], sero = id.s[m], year = id.y[m], 
                             root.dist, stringsAsFactors = FALSE)
      
      lm.tre <- lm( dist.df$root.dist ~ dist.df$year) 
      
      plot( dist.df$year, dist.df$root.dist ) 
      abline(lm.tre) 
      text( x = dist.df$year, y = dist.df$root.dist, labels = rownames(dist.df), cex = 0.5, pos = 3)
      
      out <- sort( lm.tre$residuals, index.return = TRUE, decreasing = TRUE )$ix[1: floor(1/10*( ntax) ) ]
      
      nexdata[, ncol(nexdata) + 1]             = NA
      colnames( nexdata )[ length( nexdata ) ] = "shapee"
      nexdata$shapee[ out ] = colorRampPalette( c("red", "white") )( length(out) )
      
      p = ggtree(nex, right = TRUE) %<+% nexdata + geom_tippoint( aes(color = I(shapee)), size = 3, alpha = 0.9) 
      print(p)
      lm.max <- list( R.sq    = summary( lm.tre )$r.squared, 
                      slope   = summary( lm.tre )$coefficients[2,1],
                      outlier = data.frame(o1 = dist.df$name[out], o2 = out),
                      df      = dist.df
      )
      
      ls[[ length(ls) + 1 ]] = lm.max
      
    }
    
    
  }else
  {
    fas <- read.fasta( file )
    id  <- attributes( fas )$names
    seq <- getSequence( fas )
    
    id.a <- str_match( id, "[A-Z]{1,2}[0-9]{5,6}|EPI[0-9]+" )[,1]
    id.g <- str_match( id, "\\|([A-Za-z_]+)\\|")[,2]
    id.s <- str_match( id, "\\|_(H[0-9]{1,2}N[0-9]{1,2})_")[,2]
    id.y <- as.numeric( str_match( id, "_([0-9]{4}.[0-9]{3})$")[,2] )  
    id.n <- gsub("^[A-Z]{1,2}[0-9]{5,6}_|^EPI[0-9]+_|_\\|[A-Za-z_]+\\|_|H[0-9]{1,2}N[0-9]{1,2}_[0-9]{4}.[0-9]{3}$", "", id)
    
    seq.l <- sapply( seq,
                     function(x)
                     {
                       z = grep( pattern = "a|t|c|g", 
                                 x = x, 
                                 ignore.case = TRUE, value = TRUE )
                       
                       l = length( z )
                       
                       return(l)
                     } )
    
    ls <- list( id.a, id.g, id.s, id.y, id.n, id, seq.l)
    df <- data.frame( id.a, id.g, id.s, id.y, id.n, id, seq.l )
    
    if( TRUE %in% is.na(unlist( ls[c(1:6)] ) ) ){ stop() }
    
  }
  
  if ( makecsv ){ write.csv(df, file = sub( ".fasta", "_info.csv", file) , row.names = FALSE) }
  
  return( ls )
  
  #v20181016
}


### timeDice --------------------------------

timeDice <- function( fas.dir, ecotab.dir, oldfas.dir, ecotable = TRUE )
{
  require(seqinr)
  require(tidyverse )
  
  id    <- attributes( read.fasta( fas.dir ) )$names
  id.ac <- str_match( id, "^[A-Z0-9]+")[,1]
  id.yr <- as.numeric( str_match( id, "_([0-9.]+)$")[,2] )
  
  id.0    <- unlist( sapply( as.list( oldfas.dir ), function(x) attributes( read.fasta(x) )$names ) )
  id.0.ac <- sapply( as.list( id.0), 
                     function(x)
                     {
                       if( grepl( "_EPI_ISL_", x ) )
                       { y = gsub( "_ISL_", "", str_match( x, "EPI_ISL_[0-9]+")[,1] )
                       }else
                       { y = str_match( x, "^[A-Z0-9]+")[,1] }
                       
                       return(y)
                     })
  
  id.0.yr <- ifelse( grepl( "_\\(Month_and_day_unknown\\)_|[0-9]--$", id.0 ), 1, 0 )
  
  if( ecotable )
  {
    df_ecotab <- read.table( ecotab.dir, header = TRUE, stringsAsFactors = FALSE )
    
  }else
  {
    df_ecotab <- data.frame( id = id, states = ".", stringsAsFactors = FALSE )
  }
  
  if( TRUE %in% is.na( match( id.ac, id.0.ac ) ) ){ stop( "mismatch" ) }
  if( TRUE %in% is.na( match( id, df_ecotab$id ) ) ){ stop( "mismatch" ) }
  
  
  yr   <- id.yr
  yr_0 <- paste0( floor(yr), ".000" )
  delT <- which( id.0.yr[ match( id.ac, id.0.ac ) ] == 1 ) 
  tag  <- df_ecotab$states[ match( id, df_ecotab$id ) ] 
  id_t <- id[delT]
  
  p1       <- paste0( "<taxon id=\"", id, "\">", "\n<date value=\"", yr , "\" direction=\"forwards\" units=\"years\"/>", "\n<attr name=\"states\">", tag ,"</attr>\n</taxon>")
  p1[delT] <- paste0( "<taxon id=\"", id[delT], "\">", "\n<date value=\"", yr_0[delT] , "\" direction=\"forwards\" units=\"years\" precision=\"1.0\"/>", "\n<attr name=\"states\">", tag[delT] ,"</attr>\n</taxon>") 
  p1       <- c( "# following <taxa id=\"taxa\">", p1)
  
  # after 	<taxa id="taxa">
  write.table( p1, file = sub( ".fasta", ".taxa", fas.dir ), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # before the end of treeModel
  q1 <- paste0( "<leafHeight taxon=\"", id_t, "\">\n<parameter id=\"age(", id_t, ")\"/>\n</leafHeight>" )
  q1 <- c( "# before </treeModel>",
           "<!-- START Tip date sampling                                                 -->", 
           q1, 
           "<!-- END Tip date sampling                                                   -->")
  
  write.table( q1, file = sub( ".fasta", ".treeModel", fas.dir ), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # before the end of operators 
  q2 <- paste0( "<uniformOperator weight=\"1\">\n<parameter idref=\"age(", id_t, ")\"/>\n</uniformOperator>" )
  q2 <- c( "# before </operators>", q2 )
  
  write.table( q2, file = sub( ".fasta", ".operator", fas.dir ), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # in write log to file (after ratestatistic)
  q3 <- paste0( "<parameter idref=\"age(", id_t, ")\"/>" )
  q3 <- c( "# in write log to file after ratestatistic",
           "<!-- START Tip date sampling                                                 -->", 
           q3,
           "<!-- END Tip date sampling                                                   -->" )
  
  write.table( q3, file = sub( ".fasta", ".log", fas.dir ), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  #v20181022
}



### branchAA --------------------------------

branchAA <- 
  function( trefile = file.choose(), saveFasta = FALSE , writeTre = TRUE )
  {
    require( ape )
    require( ggtree )
    require( seqinr )
    require( treeio )
    
    
    tre <- read.beast( trefile )
    tab <- fortify( tre )
    tab <- as.data.frame( tab )
    
    print( colnames(tab) )
    cat( "choose the sequence partition" )
    p = as.numeric( readline() )
    
    if( saveFasta )
    {
      
      write.fasta( names = ifelse( is.na(tab$label), paste0( "node_", tab$node ), tab$label  ),
                   sequences = as.list( tab[,p] ), 
                   file.out = gsub( ".tre", "_anc.fasta", trefile ) )
    }
    
    AA.change <- 
      sapply( as.list( tab$node ), 
              function( x )
              {
                
                if( grepl( "\\+", tab[,p][x] ) | grepl( "\\+", tab[,p][ tab$parent[ x ] ] ) )
                {
                  tab[,p][x]                 = strsplit( tab[,p][x], "+", fixed = T )[[1]][1]
                  tab[,p][ tab$parent[ x ] ] = strsplit( tab[,p][  tab$parent[ x ] ], "+", fixed = T )[[1]][1]
                }
                
                mx <- sapply( as.list( tab[,p][ c( tab$parent[ x ], x ) ] ), 
                              function(x) translate( s2c(x) ) ) 
                
                pos <- which( ! mx[,1] == mx[,2] ) 
                return( paste( paste0( mx[,1][pos], pos, mx[,2][pos] ), collapse = " " ) )
                
              })
    
    pool <- 
      unlist( sapply( as.list( tab$node ), 
                      function( x )
                      {
                        
                        if( grepl( "\\+", tab[,p][x] ) | grepl( "\\+", tab[,p][ tab$parent[ x ] ] ) )
                        {
                          tab[,p][x]                 = strsplit( tab[,p][x], "+", fixed = T )[[1]][1]
                          tab[,p][ tab$parent[ x ] ] = strsplit( tab[,p][  tab$parent[ x ] ], "+", fixed = T )[[1]][1]
                        }
                        
                        mx <- sapply( as.list( tab[,p][ c( tab$parent[ x ], x ) ] ), 
                                      function(x) translate( s2c(x) ) ) 
                        
                        pos <- which( ! mx[,1] == mx[,2] ) 
                        return( pos )
                        
                      }))
    
    AA.change <- ifelse( AA.change == "", " ", AA.change )
    
    tre@data[,p-4] = AA.change[ as.numeric( tre@data$node ) ]
    tre@data[,p-2] = tre@data$node
    tre@data[,p-1] = " "
    
    if( writeTre )
    {
      write.beast( tre, 
                   gsub( ".tre", "_aa.tre", trefile ) ) 
    }
    
    print( table( pool ) )
    #v20181123
  } 


### ha_num --------------------------------

ha_num <- function( ref_fas = file.choose(),
                    ref_csv = file.choose(),
                    data_pos = c( 2, 354, 500 ) )
{
  require( seqinr )
  require( stringr )
  
  csv  <- read.csv( ref_csv, stringsAsFactors = FALSE, header = TRUE )
  fas  <- lapply( fastaEx( ref_fas )$seq, toupper )
  
  pos = lapply( as.list(data_pos),
                function(x)
                {
                  if( x < length( which( fas[[1]] == "-" ) ) ){ return( "N-signal_peptide" ) }
                  
                  # y - > read position in the ref_fas
                  if( x > 339 ){ y = x + 7  }else{ y = x }
                  
                  # p1 - > number of h5 (vn1203)
                  p1    <- y - length( which( fas[[1]] == "-" ) )
                  
                  n.row <- which( str_match( csv$A.Vietnam.1203.2004.H5N1, "([0-9]{1,3}) ([A-Z]{1})" )[,2] == as.character(p1) )
                  p1.n  <- str_match( csv$A.Vietnam.1203.2004.H5N1, "([0-9]{1,3}) ([A-Z]{1})" )[,3][ n.row ]
                  
                  # p2 - > number of h3 (aichi)
                  p2   <- as.numeric( str_match( csv$A.AICHI.2.68.H3N2, "([0-9]{1,3}) ([A-Z]{1})" )[,2][ n.row ] )
                  p2.n <- str_match( csv$A.AICHI.2.68.H3N2, "([0-9]{1,3}) ([A-Z]{1})" )[,3][ n.row ]
                  if( is.na(p2.n) ){ p2.n = "del" }
                  
                  if( x <=  339 )
                  {
                    p3i = y - 15
                    p4i = y - 15
                    p3.n <- fas[[2]][ y ]
                    p4.n <- fas[[4]][ y ] 
                    
                    p3 = csv$p2fk0_1[ p3i ]
                    p4 = csv$p4jul_1[ p3i ]
                    
                  }else
                  {
                    p3i = y - 346
                    p4i = y - 346
                    p3.n <- fas[[3]][ y ]
                    p4.n <- fas[[5]][ y ] 
                    
                    p3 = csv$p2fk0_2[ p3i ]
                    p4 = csv$p4jul_2[ p3i ]
                    
                  }
                  
                  out = data.frame( type = c( "data", "H5", "H3", paste0( "2fk0-", ifelse( x <= 339, "A", "B") ),
                                              paste0( "4jul-", ifelse( x <= 339, "A", "B") ) ),
                                    pos  = c(  x, p1, p2, p3, p4),
                                    res  = c( fas[[6]][y], p1.n, p2.n, p3.n, p4.n ) )
                  return( out )
                  
                })
  return(pos)
  
  #v20181109
}



### viewResi --------------------------------


viewResi <- function( trefile  = file.choose(), 
                      fasfile  = file.choose(),
                      pos      = 1 )
{
  require( ggtree )
  require( seqinr )
  require( ggplot2 )
  
  tre       = read.nexus( trefile )
  tre.tip.n = gsub( "'", "", tre$tip.label )
  tre.data  = as.data.frame( fortify( tre ) )
  
  fas.n  = attributes( seqinr::read.fasta( fasfile ) )$name
  fas.s  = lapply( getSequence( seqinr::read.fasta( fasfile ) ), translate )
  fas.mx = do.call( rbind, fas.s )
  
  
  if( NA %in% match( tre.tip.n, fas.n ) ){ stop("seqeucne and tree do not match") }else
  {
    
    tre.data[ ncol( tre.data ) + 1 ]          = fas.mx[ ,pos][ match( gsub( "'", "", tre.data$label ), fas.n ) ]
    colnames( tre.data )[ ncol( tre.data )  ] = "Resi"
    
    resi.df <- as.data.frame( table( tre.data$Resi ), stringsAsFactors = FALSE )
    resi.df <- resi.df[ order(-resi.df$Freq), ]
    
    if( dim(resi.df)[1] > 8 )
    {
      
      rr = paste0( resi.df$Var[ 8: dim(resi.df[1]) ], collapse = "|")
      
      tre.data$Resi[ grep( rr, tre.data$Resi ) ] = "Others"
      tre.data$Resi <- factor( tre.data$Resi, levels = c( resi.df$Var1[1:7], "Others" ) )
      
    }else
    {
      tre.data$Resi <- factor( tre.data$Resi, levels = resi.df$Var1 )  
      
    }
  }
  
  g = 
    ggtree( tre, right = TRUE, size = 0.25, color = "gray") %<+% tre.data + 
    geom_tippoint( aes( color = Resi), alpha = 0.8 )  +
    theme_tree( legend.position = "top",
                legend.title = element_blank()) + 
    scale_color_brewer( palette =  "Dark2") +
    ggtitle( label = paste0("Resi.", pos))
  
  return( g )
  
}



### parse_hyphy --------------------------------

parse_hyphy <- function( n_sample = 1,
                         omega    = TRUE,
                         SLAC     = TRUE, 
                         MEME     = TRUE, 
                         folerdif = getwd(),
                         p        = 0.01 )
{
  require( stringr )
  require( readr )
  
  filelist = list.files( folerdif )
  
  out <-  lapply( seq(1, n_sample), list)
  
  if( omega )
  {
    if( length( which( grepl( ".wlog", filelist ) ) ) != n_sample ){ stop( "mismatch sample number" ) }
    
    file_w <- grep( "wlog", filelist, value = TRUE )
    for( i in 1: n_sample )
    {
      s1.i <- read_file( paste0( folerdif, file_w[i] ) )
      
      out[[i]][[1]] = c( as.numeric( str_match( s1.i, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3] ),
                         as.numeric( str_match( s1.i, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2] ),
                         as.numeric( str_match( s1.i, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4] ) )
    }
  }  
  
  
  if( SLAC )
  {
    if( length( which( grepl( ".slaclog", filelist ) ) ) != n_sample ){ stop( "mismatch sample number" ) }
    
    file_s <- grep( "slaclog", filelist, value = TRUE )
    for( i in 1: n_sample )
    {
      s2.i <- read_file( paste0( folerdif, file_s[i] ) )
      
      slac.str = "For partition 1 these sites are significant at p \\<\\=0.1 ([ 0-9A-Za-z-\\|.\\=:<>?\\<\\>]+)"
      
      slac.row.i <- strsplit( str_match( s2.i, slac.str )[,2], split = "\\| \\|" )[[1]][-2]
      
      out[[i]][[2]] = 
        unlist( sapply( as.list(slac.row.i[-1]), 
                        function(x)
                        {
                          codon  = as.numeric( str_match( x, "^ ([0-9]+) [\\|0-9\\. ]+ ([A-Za-z]+). p \\= ([0-9\\.]+) " )[,2] )
                          po_ne  = str_match( x, "^ ([0-9]+) [\\|0-9\\. ]+ ([A-Za-z]+). p \\= ([0-9\\.]+) " )[,3]
                          pvalue = as.numeric( str_match( x, "^ ([0-9]+) [\\|0-9\\. ]+ ([A-Za-z]+). p \\= ([0-9\\.]+) " )[,4] )
                          
                          if( (po_ne == "Pos") & (pvalue < p) ){ return(codon) }
                        } 
        ) )
      out[[i]][[3]] = slac.row.i
    }
  } 
  
  
  if( MEME )
  {
    if( length( which( grepl( ".memelog", filelist ) ) ) != n_sample ){ stop( "mismatch sample number" ) }
    
    file_m <- grep( "memelog", filelist, value = TRUE )
    for( i in 1: n_sample )
    {
      s3.i <- read_file( paste0( folerdif, file_m[i] ) )
      
      meme.str = "For partition 1 these sites are significant at p \\<\\=0.1 ([ 0-9A-Za-z-\\|.\\=:<>?\\<\\>\\+#,]+) ###"
      
      meme.row.i <- strsplit( str_match( s3.i, meme.str )[,2], split = "\\| \\|" )[[1]][-2]
      
      out[[i]][[4]] = 
        unlist( sapply( as.list(meme.row.i[-1]), 
                        function(x)
                        {
                          codon  = as.numeric( str_match( x, "^ ([0-9]+) [\\|0-9\\., A-Za-z]+ p \\= ([0-9\\.]+) \\|" )[,2] )
                          pvalue = as.numeric( str_match( x, "^ ([0-9]+) [\\|0-9\\., A-Za-z]+ p \\= ([0-9\\.]+) \\|" )[,3] )
                          
                          if(  pvalue < p  ){ return(codon) }
                        } 
        ) )
      out[[i]][[5]] = meme.row.i
    }  
  }  
  
  attr(out, "item") = c( "omega with CI", 
                         "sign. Pos. selected site (SLAC)", "full SLAC result", 
                         "sign. Pos. selected site (MEME)", "full MEME result" )
  
  names( out ) = sub( pattern = paste0( ".", c("w", "s", "m")[ which( c( omega, SLAC, MEME ) )[1] ], "log"), replacement = "", 
                      get( paste0( "file_", c("w", "s", "m")[ which( c( omega, SLAC, MEME ) )[1] ] ) ) )
  
  return( out )
  
  #v20181120
}


### idInfo.na --------------------------------

idInfo.na <- function( rawid, 
                       datasource = "n",
                       g.csv      = "",
                       na_subtype = 2 )
{
  # format: 
  # N:  >{accession}_{strain}_{serotype}_|{country}|_{year}-{month}-{day}
  # G:  Isolate name Type Collection date Isolate ID
  # Both need replace the blank with underline 
  
  require(seqinr)
  require(stringr)
  
  # g
  a.string.g <- "EPI_ISL_([0-9]+)"
  s.string.g <- "_A_/_(H[0-9]{1,2}N[0-9xX]{1,2})_"
  y.string.g <- "_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)" 
  r.string.g <- "^([0-9A-Za-z\\.\\/\\-_\\(\\)\\?]+)_A_/_"
  
  # n 
  a.string.n   <- "^[A-Za-z0-9]+"
  s.g.string.n <- "(H[0-9]{1,2}[N0-9xX]{0,2})_\\|([a-zA-Z_\\']+)\\|"
  y.string.n   <- "_[0-9]{4}-[0-9]{1,2}-[0-9]{2}$|_[0-9]{4}-[0-9]{1,2}-$|_[0-9]{4}--$|_--$"
  n.Nstring.n  <- "^[A-Za-z0-9]+_|(H[0-9]{1,2}[N0-9xX]{0,2})_\\|([a-zA-Z_\\']+)\\||_([0-9-]+)$"
  
  
  if( datasource == "g")
  {
    id.a <- gsub( "_ISL_", "", str_match( rawid, a.string.g )[, 1] )
    id.s <- str_match( rawid, s.string.g )[,2] 
    id.s[ which( id.s == paste0("N",na_subtype) | id.s == paste0("HxN", na_subtype) | id.s == paste0("HXN", na_subtype) )  ] = paste0("H0N", na_subtype)
    id.y <- str_match( rawid, y.string.g )
    id.y <- gsub( "^_", "", x = id.y)[,1]
    
    id.n <- str_match( rawid, r.string.g )[,2]
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- gsub( "_A/", "A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ] )
    
    id.n <- gsub( "\\+|\\?|\\(|\\)|\\[|\\]|\\.|:|-|/", "_", id.n )
    id.n <- gsub( "\\'|\\?|>", "", id.n )
    id.n <- gsub("A_", "", id.n)
    id.n <- gsub( "[_]+", "_", id.n )
    id.n <- gsub( "_$", "", id.n )
    id.n <- gsub( "^_", "", id.n )
    
    g <- gsub( " ", "_", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Location )
    g <- gsub( "_$", "",  str_match( g, "([A-Za-z_]+)_/_([A-Za-z_]+)" )[,3] ) 
    
    g[ which( is.na(g) == TRUE ) ] = "Unknown"
    g[ which(g == "Russian_Federation") ] = "Russia"
    g[ which(g == "United_States") ] = "USA"
    g[ which(g == "Korea") ] = "South_Korea"
    
    id.g <- g[ match( id.a, gsub("_ISL_", "", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Isolate_Id ) ) ]
    
    
  }else
  {
    id.a <- str_match( rawid, a.string.n)[,1]
    id.s <- str_match( rawid, s.g.string.n)[,2]
    id.s[ which( id.s == paste0("N",na_subtype) | id.s == paste0("HxN", na_subtype) | id.s == paste0("HXN", na_subtype) )  ] = paste0("H0N", na_subtype)
    id.g <- str_match( rawid, s.g.string.n)[,3]
    id.g[ which( id.g == "Viet_Nam") ] = "Vietnam"
    id.g[ which( id.g == "Cote_d'Ivoire") ] = "Cote_dIvoire"
    
    id.y <- str_match( string = rawid, y.string.n)
    id.y <- gsub( "_--", "1900-01-01", id.y)
    id.y <- gsub( "^_", "", id.y)
    
    id.n <- gsub( n.Nstring.n, "", rawid)
    
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- 
      paste0("A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ])
    
    id.n <- gsub("\\+|\\?|\\(|\\)|\\[|\\]|\\.|:|-|/|__", "_", id.n)
    id.n <- gsub("\\'|\\?|>", "", id.n)
    id.n <- gsub("A_", "", id.n)
    id.n <- gsub("[_]+", "_", id.n)
    id.n <- gsub("_$", "", id.n)
    id.n <- gsub("^_", "", id.n)
    
  }
  
  infolist = list(id.a, id.s, id.g, id.y, id.n)
  
  e = 
    which(
      sapply( infolist, 
              function(x)
              {
                TRUE %in% is.na(x)
                
              })  == TRUE )
  
  
  print( paste("ERROR in ", c("ac", "sero", "geo", "year", "name")[e] )  )
  
  return(infolist)
  
  #v20181205e
}


### trimtool --------------------------------

trimtool <- function( propblank = 0.8, 
                      filedir   = file.choose()){
  
  library(stringr)
  library(seqinr)
  
  file       = seqinr::read.fasta(filedir)
  seq_name0  = attributes(file)$names
  seq0       = getSequence(file)
  seq_matrix = do.call(rbind, seq0)
  
  
  coltoberemove = apply(seq_matrix, 2, 
                        function(x)
                        {
                          blank = ( length( which( x == "-") ) + length( which( x == "~") ) ) 
                          fl    = length(x)
                          
                          if ( fl*propblank < blank ){ return(1) }else{ return(0) }    
                        }
  )
  
  cut_matrix = seq_matrix[ ,-which(coltoberemove == 1) ]
  
  seq_cut    = as.list( data.frame(t(cut_matrix), stringsAsFactors = FALSE) )
  
  filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fas)" )[,2]
  
  write.fasta( seq_cut, 
               file.out = paste0("./trim_", filename, ".fasta"),
               names    = seq_name0)
  print("DONE")
  
  #201811207
}


### treeMDR --------------------------------

treeMDR <- function( index, bigmatrix )
{
  lth = length( index )
  mx  = matrix( nrow = lth, ncol = lth )
  
  for( i in 1: lth )
  {
    for( j in 1: lth )
    {
      mx[ i,j ] =  bigmatrix[ index[i], index[j] ]
    }
  }
  
  colnames( mx ) = index
  rownames( mx ) = index
  
  out <- as.data.frame( cmdscale( mx ) )
  colnames( out ) = c( "Dim_1", "Dim_2" )
  
  return( out )
  
  #v20190110
}


