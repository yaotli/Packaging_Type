### fastaEx --------------------------------

fastaEx <- function(filedir = file.choose())
{
  require(seqinr)
  
  file     <- read.fasta(filedir)
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
    id.n <- gsub( "__", "_", id.n )
    id.n <- gsub( "\\'|\\?|>", "", id.n )
    id.n <- gsub("A_", "", id.n)
    id.n <- gsub( "_$", "", id.n )
    
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
    id.n <- gsub("__", "_", id.n)
    id.n <- gsub("\\'|\\?|>", "", id.n)
    id.n <- gsub("A_", "", id.n)
    id.n <- gsub("_$", "", id.n)
    
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
  
  #v20180912e
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
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
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
                    consenseq = FALSE )
{
  require( seqinr )
  
  leaflist <- unlist( strsplit( gsub( " ", "", leaflist ), "\n" ) )
  
  fas   <- read.fasta( filedir )
  fas.s <- getSequence( fas )
  fas.n <- attributes( fas )$names
  
  m <- match( leaflist, fas.n )
  
  if( TRUE %in% is.na(m) ){ stop() }else
  {
    write.fasta( fas.s[ m ], fas.n[ m ], gsub( pattern = ".fasta", 
                                               replacement = paste0( "_", length(m), ".fasta"), 
                                               x = filedir ) )  
    
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
  
  #v20180921
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
  
  s     <- grep( sero, str_match( id, "_(H5N[0-9]{1,2})_" )[,2], invert = sero.invert )
  
  if ( TRUE %in% is.na( 
    c(str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], 
      str_match( id, "\\|([A-Za-z_]+)\\|" )[,2], 
      str_match( id, "_(H5N[0-9]{1,2})_" )[,2]) ) )
  {
    stop( 
      c("Year", "Geo", "Serotype")[ ceiling( 
        which( is.na( 
          c(str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], 
            str_match( id, "\\|([A-Za-z_]+)\\|" )[,2], 
            str_match( id, "_(H5N[0-9]{1,2})_" )[,2]) ))/3) ] 
    )
  }  
  
  remain <- sort( Reduce( intersect, list( remain, y, g, s ) ) )
  
  write.fasta( seq[ remain ], 
               id[remain], 
               file.out = sub(".fasta", "_rmd.fasta", fasfile) )
  
  print( length( remain ) )
  
  #v20181001
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
  
  if( TRUE %in% is.na(as.numeric( str_match( fas.i0, "_([0-9]{4}.[0-9]{3})$" )[,2] )) ){ stop() }
  
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
  
  #v20181001
}
