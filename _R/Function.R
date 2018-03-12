# Functions applied for this project

## phylo_date ----------------

# convert YYYY-MM-DD to YY.date

phylo_date <- function(x)
{
  library(stringr)
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  yr   <- as.numeric(str_match(x, d)[,2])
  mo   <- as.numeric(str_match(x, d)[,3])
  dy   <- as.numeric(str_match(x, d)[,4])
  
  yr.0 <- paste0(yr, "-01-01")
  
  daydifference <- as.numeric( difftime( strptime( x, format = "%Y-%m-%d"),
                                         strptime( yr.0, format = "%Y-%m-%d"), 
                                         units = "days"
                                         ) 
                               )/365
  
  # bug?
  if ( TRUE %in% is.na(daydifference) )
  {
    x   <- sub(pattern = "01$", replacement = "02", x)
    
    daydifference <- as.numeric( difftime( strptime( x, format = "%Y-%m-%d"),
                                           strptime( yr.0, format = "%Y-%m-%d"), 
                                           units = "days"
                                           )
                                 )/365
  }
    
  yr.daydifference <- yr + daydifference
  yr.daydifference <- format( round( yr.daydifference, 3 ), nsmall = 3)
  
  return(yr.daydifference)
  
  #v20170908
}
  

### cleanID --------------------------------

cleanID <- function(filedir = file.choose())
{
  
  # require Function.R / phylo_date   
  
  library(seqinr)
  library(stringr)
  
  # read-in
  
  file      <- read.fasta(filedir)
  seq_name0 <- attributes(file)$names
  seq0      <-  getSequence(file)
  
  # deal with ineligible id
  
  seq_name  <- gsub(" ", "_", seq_name0)
  seq_name  <- gsub("\\(|\\)|\\[|\\]|\\.|:", "_", seq_name)
  seq_name  <- gsub("\\'|\\?|>", "", seq_name)
  
  # for gisaid
  
  seq_name  <-  gsub("_A_/", "", seq_name)
  
  # deal with date
  
  # for gisaid
  seq_name  <-  gsub("__Month_and_day_unknown_", "-07-01", seq_name)
  seq_name  <-  gsub("__Day_unknown_", "-15", seq_name)
  
  # for ncbi
  # w/o all date info
  seq_name[ which( endsWith(seq_name, "_--") == "TRUE") ] <- 
    
    gsub("_--", "_9999-99-99", seq_name[which( endsWith(seq_name, "_--") == "TRUE")])
  
  # w/o month and day
  seq_name[ which( endsWith(seq_name, "--") == "TRUE") ]  <-
    
    gsub("--", "-07-01", seq_name[which( endsWith(seq_name, "--") == "TRUE")])
  
  # w/o day
  seq_name[ which( endsWith(seq_name, "-") == "TRUE" ) ]   <- 
    
    paste0( seq_name[which( endsWith(seq_name, "-") == "TRUE" ) ], "15" )
  
  
  # generate temporal info
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  for (i in 1: ( length(seq0)) )
  {
    seq_name[i] <- gsub( x           = seq_name[i], 
                         pattern     = d, 
                         replacement = phylo_date( 
                           str_match(string = seq_name[i], pattern = d)[,1] ) 
    )
    
  }
  
  # deal with duplicated ID  
  
  duplicated_id = which(duplicated(seq_name) == "TRUE")
  
  # loop for adding pseudocode
  
  if( length(duplicated_id) > 0 )
  {
    print("There're duplicated sequences!")
    # for all the duplicated id
    
    for (i in 1: length(duplicated_id))
    {
      dup0   = which( match(seq_name, seq_name[duplicated_id[i]] ) != "NA" )
      
      # create a null vector to appendex  
      
      app    = c("", paste0(letters, "_"), toupper(paste0(letters, "_") ) )
      ap.id  = seq_name[dup0]
      app.id = c()
      
      
      # loop to deal with multiple replicated
      # find the date info, insert labeling in the middle
      
      d        = "([0-9]{4})\\.([0-9]{3})"
      time_rep <- str_match(ap.id, d)[,1]
      
      for (k in 1: length(ap.id))
      {
        app.id[k] = sub(time_rep[k], paste0(app[k], time_rep[k]), ap.id[k])
      }
      
      # back to seq_name
      
      seq_name[dup0] = app.id  
      
    }
    
  }
  
  duplicated_id_ed = which(duplicated(seq_name) == "TRUE")  
  
  if ( length( duplicated_id_ed ) > 0 )
  {
    print("ERROR!") 
  }else
  {
    filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)" )[,2]
    
    seq_name <- gsub("-|/", "_", seq_name)
    
    # sort by time
    
    year0  <- as.numeric( gsub("_", "", str_match(seq_name, "_([0-9]{4})\\.([0-9]+)")[,1]) )
    
    sortid <- sort( year0, na.last = FALSE, index.return = TRUE )$ix
    
    write.fasta( seq0[sortid], 
                 file.out = paste0("./cleanID_", filename, ".fasta"), 
                 names = seq_name[sortid])
    
    print("DONE")
    
  }
  #20170728    
} 

  


### curateSeq --------------------------------

curateSeq <- function(maxamb  = 5, 
                      minseq  = 1600, 
                      mode    = 1, 
                      vip     = 0, 
                      filedir = file.choose()  )
{
  # should apply after cleanID
  # mode 1 : curation; 2 : duplicated seq; 3 : duplicated id; 4 : similar id
  #      5 : 1x2x3 ;   6 : 1x2x4
  #      7 : 1+2+3 ;   8 : 1+2+4    
  
  library(seqinr)
  library(stringr)
  
  # read-in
  
  file = read.fasta(filedir)
  
  seq_name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  # sequence curation [1]
  
  # for ambiguous nucleotide and seq length
  
  tobedelect1 <- c()
  
  for(i in 1: length(seq0))
  {
    ATCG = c("a", "t", "c", "g")
    
    # convert character to string
    
    seq_0 = c2s( seq0[[i]] )
    seq_i = gsub("-", "", seq_0)
    seq_i = gsub("~", "", seq_i)
    
    # seq length and number of ambiguous nucleotide
    
    seqlth = length( s2c(seq_i) )
    amb    = length( which(! s2c(seq_i) %in% ATCG ) )
    
    if( ( seqlth < minseq ) | ( amb > maxamb ) )
    {
      tobedelect1[length(tobedelect1) + 1 ] = i
    }
    
  }
  print( paste0("curate: ", length(tobedelect1) ) )  
    
  # duplicated sequcnes [2]
  
  dup_seq <-  which(
    
    duplicated( 
      sapply(seq0, 
             function(x)
             {
               y = c2s(x)
               z = gsub("-", "", y)
               z = gsub("~", "", z)
               
               return(z) 
             } 
             
             )) == "TRUE")  
  
  print( paste0("dupSeq: ", length(dup_seq) ) )
  
  # duplicated id (after cleanID) [3]
  
  du     = "_([HN0-9]+)_([a-zA-z])_"
  dup_id <- grep(du, seq_name0)
  
  print( paste0("dupID: ", length(dup_id) ) )
  
  # similar id [4]
  
  m = "a/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)"
  
  sim_id <- which( duplicated(str_match( tolower(seq_name0), m)[,1]) == TRUE)
  
  print( paste0("simID: ", length(sim_id) ) )
  
  # mode 
  
  mode1x2   <- sort( unique( c(tobedelect1, dup_seq) ) )
  mode1x3   <- sort( unique( c(tobedelect1, dup_id ) ) )
  mode1x4   <- sort( unique( c(tobedelect1, sim_id ) ) )
  
  mode1x2x3 <- sort( unique( c(tobedelect1, intersect(dup_seq, dup_id)) ) )
  mode1x2x4 <- sort( unique( c(tobedelect1, intersect(dup_seq, sim_id)) ) )
  
  mode1p2p3 <- sort( unique( c(tobedelect1, dup_seq, dup_id ) ) )
  mode1p2p4 <- sort( unique( c(tobedelect1, dup_seq, sim_id ) ) )
  

  tobedelect_f <- list(tobedelect1, mode1x2, mode1x3, mode1x4, 
                       mode1x2x3, mode1x2x4, 
                       mode1p2p3, mode1p2p4)
  
  remain       <- seq(1: length(seq0))[- tobedelect_f[[mode]] ]
  
  
  # special list
  
  remain       <- sort( unique( c( remain, rev(seq(1: length(seq0) ))[1:vip] ) ) )
  
  seq_name_out <- seq_name0[remain]
  seq_out      <- seq0[remain]
  
  filename     <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)" )[,2]
  
  write.fasta(seq_out,
              file.out = paste0("./curateSeq", "-", mode, "_", filename, ".fasta"), 
              names = seq_name_out)
  
  print( paste0("delete: ", length(tobedelect_f[[mode]]) ) )
  print( paste0("remain: ", length(remain) ) )
  
  #v20170728
}
  
  

### findtaxa --------------------------------

findtaxa <- function(type, 
                     tree, 
                     targetid, 
                     target)
{
  
  # type 1 = branch coloring 
  # type 0 = tip shape
  # default branch color = black
  
  library(ape)
  library(ggtree)
  
  # extract tree data
  
  tree.d                         <- fortify(tree)
  tree.d[, ncol(tree.d) + 1]     <- gsub("'", "", tree.d$label)
  colnames(tree.d)[ncol(tree.d)] <- "taxaid"
  
  # for tip shape
  
  if (type == 0)
  {
    
    shapetaxa <- data.frame( node   = c( 1:length( tree.d$isTip ) ), 
                             shapee = NA)
    
    for (i in 1: length(targetid) )
    {
      shapetaxa$shapee[  grep(  tolower( targetid[i] ), tolower(tree.d$taxaid) ) ] <- target[i]
    }
    
    return(shapetaxa)
    
  }else {
    
    # for branch colorring  
    
    # new column
    tree.d[, ncol(tree.d) + 1]     <- "black"
    colnames(tree.d)[ncol(tree.d)] <- "colorr"
    
    # for branch extension
    edgematrix <- as.matrix( tree.d[,c(2,1)] )
    
    # color grouping 
    group_color <- unique( target )
    
    for ( i in 1: length( group_color ) )
    {
      
      # color as group to combine key word to targetno
      sub_color <- which( target == group_color[i] )
      targetno  <- c()
      
      for ( t in 1: length( sub_color ) )
      {
        targetno <- 
          unique( c( targetno, grep( tolower( targetid[ sub_color[t] ] ), 
                                     tolower( tree.d$taxaid ) )
          ) )
      }
      
      pre_targetno  <- length( targetno )
      post_targetno = 0
      
      # while loop 
      
      while( pre_targetno != post_targetno )
      {
        
        pre_targetno = length( targetno )
        
        for( k in 1: length(targetno) )
        {
          # all sibiling 
          sibs <- as.integer( edgematrix[ 
            which( edgematrix[,1] == 
                     edgematrix[ which( edgematrix[,2] == targetno[k] ), ][1] ),][,2] )
          
          if ( length( sibs ) == 1 )
          {
            targetno = c( targetno, as.integer( edgematrix[which( edgematrix[,2] == targetno[k]), ][1] ) )
            
          }else{
            
            if ( length( which(! sibs %in% targetno ) ) == 0 )
            {
              # add parent no
              targetno = c( targetno, as.integer( edgematrix[which( edgematrix[,2] == targetno[k]),][1] ) )
            }
          }
          
          targetno = unique( targetno )
        }
        
        post_targetno = length(targetno)
        
      }
      
      # coloring
      
      tree.d$colorr[targetno] <- group_color[i]
      
    }
    return(tree.d)    
    
  }
  #v20170929
}


### subtreeseq --------------------------------

subtreseq<-function(findrep      = 0, 
                    outlier      = 0, 
                    originfile   = 0, 
                    seq_filedir  = file.choose(), 
                    list_filedir = file.choose(),
                    ori_filedir  = file.choose(), 
                    no          = "")
{
  library(seqinr) 
  library(stringr)
  
  # readin
  fasta0    = read.fasta( seq_filedir )
  seq.name0 = attributes(fasta0)$names
  seq0      = getSequence(fasta0)
  filename  = str_match( seq_filedir, "([a-zA-Z0-9-_]+)(\\.)(fasta)")[,2] 
  
  sub.tree   = read.table( list_filedir, header = FALSE, stringsAsFactors = FALSE)
  id.subtree = match( sub.tree[,1], seq.name0)
  
  if (originfile == 1)
  {
    fasta_ori     = read.fasta( ori_filedir )
    seq.name0_ori = attributes( fasta_ori )$names
    seq0_ori      = getSequence( fasta_ori )
  }
  
  # error
  if( any(NA %in% id.subtree) == "TRUE")
  {
    print("Mismatch")
    
  }else
  {
    # findrep = 1
    
    if ( findrep == 1 )
    {
      dup_seq <-  which( duplicated( sapply(seq0, 
                                            function(x)
                                            {
                                              y = c2s(x)
                                              z = gsub("-", "", y)
                                              z = gsub("~", "", z)
                                              
                                              return(z) 
                                            } )) == "TRUE")
      
      dup     <- c()
      
      for ( k in 1: length(dup_seq) )
      {
        tomatch <- gsub("~", "", gsub("-", "", c2s( seq0[[ dup_seq[k] ]] ) ))
        
        dup_i   <- which(sapply(seq0, 
                                function(x)
                                {
                                  y = c2s(x)
                                  z = gsub("-", "", y)
                                  z = gsub("~", "", z)
                                  
                                  return(z)
                                } ) %in% tomatch  == TRUE)
        
        
        if ( any( seq.name0[dup_i] %in% sub.tree[,1] == TRUE) )
        {
          dup = c(dup, dup_i)
        }
        
      }
      
      id.subtree        = sort( unique( c( id.subtree, dup )) )
      
      seq.name0_subtree = seq.name0[id.subtree]
      seq0_subtree      = seq0[id.subtree]
      
      write.fasta(seq0_subtree, 
                  file.out = paste0("./", no, "_", filename, "_", "subtree.fasta"), 
                  names    = seq.name0_subtree )  
      
      # originfile = 1
      
      if (originfile == 1)
      {
        seq.name0_subtree = seq.name0_ori[id.subtree]
        seq0_subtree      = seq0_ori[id.subtree]
        
        write.fasta(seq0_subtree, 
                    file.out = paste0("./", no, "_", filename, "_", "Ori_subtree.fasta"), 
                    names = seq.name0_subtree )
      }
      
      if (outlier == 1)
      {
        id.outlier        <- seq(1, length(seq0))[-id.subtree]
        
        seq.name0_subtree = seq.name0[id.outlier]
        seq0_subtree      = seq0[id.outlier]
        
        write.fasta(seq0_subtree, 
                    file.out = paste0("./", no, "_", filename, "_", "subtreeB.fasta"), 
                    names    = seq.name0_subtree )
        
        if (originfile == 1)
        {
          seq.name0_subtree = seq.name0_ori[id.outlier]
          seq0_subtree      = seq0_ori[id.outlier]
          
          write.fasta(seq0_subtree, 
                      file.out = paste0("./", no, "_", filename, "_", "Ori_subtreeB.fasta"), 
                      names    = seq.name0_subtree )
        }
      }
      
      print("Done")    
      
    }else
    {
      seq.name0_subtree = seq.name0[id.subtree]
      seq0_subtree      = seq0[id.subtree]
      
      write.fasta(seq0_subtree, 
                  file.out = paste0("./", no, "_", filename, "_", "subtree.fasta"), 
                  names    = seq.name0_subtree )
      
      if (outlier == 1)
      {
        id.outlier        <- seq(1, length(seq0) )[-id.subtree]
        
        seq.name0_subtree = seq.name0[id.outlier]
        seq0_subtree      = seq0[id.outlier]
        
        write.fasta(seq0_subtree, 
                    file.out = paste0("./", no, "_", filename, "_", "subtreeB.fasta"), 
                    names    = seq.name0_subtree )
      }
      print("Done")
    }
    
  }
  #v20170728
}



## gg_color_hue ----------------

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n] }



### trimtool --------------------------------

trimtool <- function( propblank = 0.8, 
                      filedir   = file.choose()){
  
  library(stringr)
  library(seqinr)
  
  file       = read.fasta(filedir)
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

  #20170924  
}


## keepLongSeq ----------------


keepLongSeq <- function(seq_0, 
                        id_0, 
                        showRemain = FALSE)
{
  library(seqinr)
  
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


### subfastaSeq --------------------------------

subfastaSeq <- function( subtype         = "H5N1",
                         inverse.subtype = FALSE,
                         time_s          = 1000,
                         time_e          = 3000,
                         AC              = FALSE, 
                         AC_list         = NA, 
                         geo             = "all",
                         inverse.geo     = FALSE, 
                         host            = "all",
                         inverse.host    = FALSE, 
                         no              = "",
                         filedir         = file.choose() 
)
{
  require(seqinr)
  require(stringr)
  
  file  = read.fasta(filedir)
  
  seq_name0 <- attributes(file)$names
  seq0      <- getSequence(file)
  
  # terms ----------------
  
  nohuman  <- 
    "avian|bird|wildbird|chicken|duck|dove|pigeon|mallard|goose|environment|water|teal|hawk|magpie|munia|myna|kestrel|peregrine|crow|sparrow|mesia|gull|egret|swan|shrike|buzzard|heron|Ph|quail|pheasant|grebe|starling|swallow|white_eye|swine|tiger"
  
  nomammal <- 
    "avian|bird|wildbird|chicken|duck|dove|pigeon|mallard|goose|environment|water|teal|hawk|magpie|munia|myna|kestrel|peregrine|crow|sparrow|mesia|gull|egret|swan|shrike|buzzard|heron|Ph|quail|pheasant|grebe|starling|swallow|white_eye"
  
  N    <- 
    "Shanxi|Hebei|Beijing|Jilin|Sheny|Liaoning|Heilongjiang|North_China"
  
  E    <-
    "Shandong|Jiangsu|Huadong|Eastern_China|Fujian|Anhui|Shanghai|Zhejiang|Jiangxi|JX|Nanchang|Zaozhuang|Dongtai|Xuzhou|Danyang|Yangzhou"
  
  C    <- 
    "Hunan|Hubei|Henan|Changsha"
  
  S    <- 
    "Hong_Kong|Shantou|Guangdong|GD|ST|Shenzhen|Guangzhou"
  
  SW   <- 
    "Guizhou|Guangxi|Yunnan|Guiyang|Tibet|Sichuan|Chongqing|Anning|Dali|TongHai"
  
  NW   <- 
    "Ningxia|Xinjiang|Gansu|Qinghai|Shaanxi"
  
  
  # main ----------------
  
  host_string <- c( "nohuman", "nomammal", "all")
  host_v      <- c( nohuman, nomammal, "")
  geo_string  <- c( "N", "E", "C", "S", "SW", "NW", "all")
  geo_v       <- c( N, E, C, S, SW, NW, "")
  
  host_key    <- grep( host, host_string )
  geo_key     <- grep( geo, geo_string )
  
  if ( AC )
  {
    ac.file  <- str_match( seq_name0, "^[A-Z0-9]+" )
    ac.i     <- match( AC_list, ac.file)
    
    if ( TRUE %in% is.na(ac.i) )
    {
      print("ERROR")
      
    }else
    {
      filename <- str_match( filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)" )[,2]
      
      write.fasta(sequences = seq0[ac.i], 
                  names     = seq_name0[ac.i],
                  file.out  = paste0("./", "ac", no, "_", filename, ".fasta") )
      
      print( paste0("no: ", length( ac.i ) ) )
      
    }
    
  }else
  {
    
    if( inverse.subtype )
    {
      subtype_i <- grep( pattern = paste0("_", subtype, "_"), 
                         x       = seq_name0, 
                         invert  = TRUE)   
    }else
    {
      subtype_i <- grep( pattern = paste0("_", subtype, "_"), 
                         x       = seq_name0)        
    }
    
    
    # time of isolation
    
    T_iso  <- as.numeric( gsub("_", "", str_match( seq_name0, "_([0-9]{4})\\.([0-9]+)" )[,1] )  )
    
    time_i <- which(T_iso < time_e & T_iso > time_s)
    
    
    # host & geo
    
    if ( inverse.geo )
    {
      g_i    <- grep( pattern = geo_v[ geo_key ], x = seq_name0, invert = TRUE)
      
    }else
    {
      g_i    <- grep( pattern = geo_v[ geo_key ], x = seq_name0)
      
    }
    
    if( inverse.host ){
      
      h_i    <- grep( pattern = host_v[ host_key ], x = seq_name0, invert = TRUE, ignore.case = TRUE)
      
    }else
    {
      h_i    <- grep( pattern = host_v[ host_key ], x = seq_name0, ignore.case = TRUE) 
      
    }
    
    # output 
    
    remain   <- sort( Reduce( intersect, list( subtype_i, time_i, h_i, g_i ) ) )
    
    filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)")[,2]
    
    write.fasta( sequences = seq0[remain], 
                 names     = seq_name0[remain],
                 file.out  = paste0("./", filename, "_", subtype, "-", time_s, "-", time_e,"_", host,"_", geo, "_", no, ".fasta") )
    
    print( paste0("no: ", length( remain ) ) )
    
    
  }
  
  
  #v20171205
}


### nucleotide partition --------------------------------

ntpartition <- function( position   = c(1:100), 
                         filedir    = file.choose(), 
                         no         = NULL)
{
  library(seqinr)
  library(stringr)
  
  file      = read.fasta(filedir)
  seq_name0 = attributes(file)$names
  seq0      = getSequence(file)
  
  seq_matrix <- do.call( rbind, seq0 )
  
  if( nrow(seq_matrix)[1] == 1 )
  {
    seq0[[1]] <- seq0[[1]][position]
    cut_list  <- seq0
    
  }else
  {
    cut_matrix <- seq_matrix[, position]
    cut_list   <- as.list( data.frame( t(cut_matrix), stringsAsFactors = FALSE) )  
  }
  
  
  filename   <- str_match( filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)")[,2]
  
  write.fasta(sequences = cut_list, 
              file.out  = paste0("./p-", filename, no),
              names     = seq_name0)
  
  print( dim(cut_matrix)[2] )
  
  #v20170812
}



### remove ambiguous nucleotide --------------------------------

rmAMBnt <- function( filedir = file.choose() )
{
  require(seqinr)
  
  filein <- read.fasta( filedir )
  name   <- attributes( filein )$names
  seq0   <- getSequence( filein )
  
  seq    <- lapply(seq0, 
                   function(y)
                   {
                     amb    <- grep(pattern = "a|t|c|g|-", 
                                    x = y, invert = TRUE, ignore.case = TRUE)
                     y[amb] <- "-"
                     
                     return(y)
                   })
  
  write.fasta(seq, 
              name, 
              file.out = filedir)
  #v20170812 
}

### random sampling of sequence --------------------------------

rSeq <- function(n       = 10, 
                 seed    = 1, 
                 s_times = 1,
                 filedir = file.choose())
{
  library(seqinr)
  library(stringr)
  
  file      <- read.fasta( filedir )
  seq_name0 <- attributes( file )$names
  seq0      <- getSequence( file )
  filename  <- str_match( filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)")[,2]
  
  for (k in 1: s_times )
  {
    set.seed( seed*k )     
    y            <- sample( 1: length(seq0), n )  
    y            <- sort( unique(y) )
    out_seq      <- seq0[y]
    out_seq_name <- seq_name0[y]
    
    write.fasta( out_seq, 
                 names = out_seq_name, 
                 file.out = paste0("./", filename, "_r", k, ".fasta") )
  }
  
  print("DONE") 
  
  #v20170815
}



### Extract epi info from sequence name --------------------------------

idInfo <- function( rawid, 
                    datasource = "n",
                    g.csv      = "")
{
  # format: 
  # N:  >{accession}_{strain}_{serotype}_|{country}|_{year}-{month}-{day}
  # G:  Isolate name Type Collection date Isolate ID
  # Both need replace the blank with underline 
  
  library(seqinr)
  library(stringr)
  
  # g
  a.string.g  <- "EPI_ISL_([0-9]+)"
  s.string.g  <- "_A_/_(H5N[0-9xX]{1,2})_"
  y.string.g  <- "_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)" 
  n.Nstring.g <- "_EPI_ISL_([0-9]+)|_A_/_H5N[0-9xX]{1,2}|_[0-9]{4}[-0-9]{6}|_[0-9]{4}-[0-9]{2}_\\(Day_unknown\\)|_[0-9]{4}_\\(Month_and_day_unknown\\)"
  
  # n 
  a.string.n   <- "[A-Z]{1,2}[0-9]{5,6}"
  s.g.string.n <- "_(H5[N0-9xX]{0,2})_\\|([a-zA-Z_\\']+)\\|"
  y.string.n   <- "_[0-9]{4}-[0-9]{2}-[0-9]{2}|_[0-9]{4}-[0-9]{2}-|_[0-9]{4}--|_--"
  n.Nstring.n  <- "[A-Z]{1,2}[0-9]{5,6}_|_(H5[NxX0-9]{0,2})_\\|([a-zA-Z_\\']+)\\||_([0-9]{4}[-0-9]{2,6})"
  
  if( datasource == "g")
  {
    id.a <- gsub("_ISL_", "", str_match( rawid, a.string.g )[, 1] )
    
    id.s <- str_match( rawid, s.string.g )[,2] 
    id.s[ which( id.s == "H5"|
                 id.s == "H5Nx"|
                 id.s == "H5NX")  ] = "H5N0"
    
    id.y <- str_match( rawid, y.string.g )
    id.y <- gsub( "^_", "", x = id.y)[,1]
    
    id.n <- gsub( n.Nstring.g, rawid, replacement = "")
    
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- gsub( "_A/", "A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ] )
    
    id.n  <- gsub("\\(|\\)|\\[|\\]|\\.|:|-|/", "_", id.n)
    id.n  <- gsub("__", "_", id.n)
    id.n  <- gsub("\\'|\\?|>", "", id.n)
    id.n  <- gsub("A_", "", id.n)
    id.n  <- gsub("_$", "", id.n)
    
    g <- gsub( " ", "_", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Location )
    g <- gsub("_$", "",  str_match( g, "([A-Za-z_]+)_/_([A-Za-z_]+)" )[,3] ) 
    
    g[ which( is.na(g) == TRUE ) ] = "Unknown"
    g[ which(g == "Russian_Federation") ] = "Russia"
    g[ which(g == "United_States") ] = "USA"
    g[ which(g == "Korea") ] = "South_Korea"
    
    id.g <- g[ match( id.a, gsub("_ISL_", "", read.csv( g.csv, header = TRUE, stringsAsFactors = FALSE)$Isolate_Id ) ) ]
    
    
  }else
  {
    id.a <- str_match( rawid, a.string.n)[,1]
    
    id.s <- str_match( rawid, s.g.string.n)[,2]
    id.s[ which( id.s == "H5"|
                 id.s == "H5Nx"|
                 id.s == "H5NX")  ] = "H5N0"
    
    id.g <- str_match( rawid, s.g.string.n)[,3]
    id.g[ which( id.g == "Viet_Nam") ] = "Vietnam"
    id.g[ which( id.g == "Cote_d'Ivoire") ] = "Cote_dIvoire"
    
    id.y <- str_match( string = rawid, y.string.n)
    id.y <- gsub( "_--", "1900-01-01", id.y)
    id.y <- gsub( "^_", "", id.y)
    
    id.n <- gsub( n.Nstring.n, "", rawid)
    
    id.n[ which( startsWith(id.n, "A/") == FALSE) ] <- 
      paste0("A/", id.n[ which( startsWith(id.n, "A/") == FALSE) ])
    
    id.n  <- gsub("\\(|\\)|\\[|\\]|\\.|:|-|/|__", "_", id.n)
    id.n  <- gsub("__", "_", id.n)
    id.n  <- gsub("\\'|\\?|>", "", id.n)
    id.n  <- gsub("A_", "", id.n)
    id.n  <- gsub("_$", "", id.n)
    
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
  
  #v20171003
}


### remove duplicated strain --------------------------------

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


### time-extraction from sequence id  --------------------------------


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
  
  daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d"),
                                         strptime( yr.0, "%Y-%m-%d"), 
                                         units = "days") 
  )/yday(yr.e)
  
  # bug?
  if ( TRUE %in% is.na(daydifference) )
  {
    
    rawdata.2[ which(is.na(daydifference)) ] <- 
      sub("01$", "02", rawdata.2[ which(is.na(daydifference)) ] )
    
    
    daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d"),
                                           strptime( yr.0, "%Y-%m-%d"), 
                                           units = "days") 
    )/yday(yr.e)
  }
  
  yr.daydifference <- yr + daydifference
  yr.daydifference <- format( round( yr.daydifference, 3 ), nsmall = 3)
  
  return(yr.daydifference)
  
  #v20180125
}


### curate the sequence  --------------------------------

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


### extract of labeled tip --------------------------------

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


### remove duplicate and filter seq based on arranged name --------------------------------

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
                  endsWith( str_match( id, "_([0-9]{4}\\.[0-9]+$)")[,2], ".499" ) )
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
               file.out = sub(".fasta", "_cr.fasta", fasfile) )
  
  print( length( remain ) )
  
  #v20171019
}

### extract id info --------------------------------

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
    id.s <- str_match( id, "_(H5N[0-9]{1,2})_")[,2]
    id.y <- as.numeric( str_match( id, "_([0-9]{4}.[0-9]{3})$")[,2] )  
    id.n <- gsub("^[A-Z]{1,2}[0-9]{5,6}_|^EPI[0-9]+_|_\\|[A-Za-z_]+\\|_|H5N[0-9]{1,2}_[0-9]{4}.[0-9]{3}$", "", id)
    
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
    id.s <- str_match( id, "_(H5N[0-9]{1,2})_")[,2]
    id.y <- as.numeric( str_match( id, "_([0-9]{4}.[0-9]{3})$")[,2] )  
    id.n <- gsub("^[A-Z]{1,2}[0-9]{5,6}_|^EPI[0-9]+_|_\\|[A-Za-z_]+\\|_|H5N[0-9]{1,2}_[0-9]{4}.[0-9]{3}$", "", id)
    
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
  
  #v20171108
}

### cladeSampling --------------------------------

cladeSampling <- function( trefile      = file.choose(),
                           fasfile      = file.choose(),
                           listinput    = list(),
                           seed         = 666,
                           grid         = 1,
                           minBranchlth = TRUE, 
                           showTree     = FALSE, 
                           saveFasta    = FALSE,  
                           suppList     = FALSE,
                           list.x       = c("id", "y", "geo") )
{
  require( ape )
  require( seqinr )
  require( ggtree )
  require( stringr )
  
  # getdescendants
  getDes <- function( node, curr = NULL )
  {
    if( is.null(curr) ){ curr <- vector() }
    
    edgemax   <- tre.d[ c(2,1) ]
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
  
  tre.d   <- fortify( nex )
  N.tip   <- length( which( tre.d$isTip ) )
  N.node  <- nex$Nnode
  edgemax <- tre.d[ c(2,1) ]
  
  
  t.id <- gsub("'", "", tre.d$label)
  
  if( suppList )
  {
    
    tem.m  <- match( t.id[ 1:  N.tip], listinput[[ list.x[1] ]])
    
    if( TRUE %in% is.na(tem.m) ){stop()}
    
    i.id.y <- listinput[[ list.x[2] ]][ tem.m ]
    i.id.g <- listinput[[ list.x[3] ]][ tem.m ]
    
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
  
  remain <- c( nogroup_tip, selected_tip )
  seq.o  <- seq[ match( t.id[ sort(remain) ], id) ]
  id.o   <- id[ match( t.id[ sort(remain) ], id) ]
  
  if( saveFasta )
  {
    write.fasta( seq.o, id.o, 
                 file.out = gsub( ".fasta", "_s.fasta", fasfile) )
  }
  
  print( paste0("sampled n = ", length(remain), " from ", length(id) ) )
  print( table( floor( i.id.y[remain] ), i.id.g[remain] ) )
  
  
  #v20171111b
}



### geoID --------------------------------


geoID <- function( strings, 
                   SEA  = SEA,
                   nA   = nA,
                   E    = E, 
                   host = FALSE)
{
  cn.NE <- "Jilin|Sheny|Liaoning|Heilongjiang"
  cn.BH <- "Beijing|Hebei|Shandong|Zhaozhuang"
  cn.YZ <- "Jiangsu|Zhejiang|Shanghai|Gaoyou|Dongtai|Xuzhou|Danyang|Yangzhou|Yuhang|Hangzhou|Wuxi|Wenzhou|Zhenjiang|Taizhou|Taishun"
  
  cn.C  <- "Hunan|Hubei|Henan|Changsha|Chang_Sha|Anhui|Jiangxi|Shanxi|Nanchang|Wuhan|Ganzhou"
  
  cn.SW <- "Guizhou|Guangxi|Yunnan|Guiyang|Tibet|Sichuan|Chongqing|Anning|Dali|TongHai"
  cn.NW <- "Ningxia|Xinjiang|Gansu|Qinghai|Shaanxi"
  cn.S  <- "Hong_Kong|Shantou|Guangdong|Shenzhen|Guangzhou|Fujian|Dongguan|Zhanjiang|Zhongshan|Huizhou|Qingyuan"
  
  SEA   <- "Bangladesh|Laos|Malaysia|Myanmar|Thailand|Vietnam|India|Nepal"
  nA    <- "Japan|Mongolia|Russia|South_Korea"
  E     <- "Bulgaria"
  
  nonML <- "avian|bird|wildbird|poultry|chicken|duck|dove|pigeon|mallard|goose|environment|water|teal|hawk|eagle|magpie|munia|myna|kestrel|peregrine|crow|sparrow|robin|mesia|gull|egret|swan|shrike|buzzard|heron|quail|pheasant|grebe|starling|swallow|white_eye|cormorant|goldeneye|fowl|pochard|crane|peacock|turkey|falcon|swiftlet|rook|pintail|curlew"
  
  
  geo.key <- c( cn.NE, cn.BH, cn.YZ, cn.C, cn.SW, cn.NW, cn.S, SEA, nA, E)
  geo     <- c("cnNE", "cnBH", "cnYZ", "cnC", "cnSW", "cnNW", "cnS", "SEA", "nA", "E", "Unknown")
  
  if( host )
  {
    geo.key <- c(nonML) 
    geo     <- c("nonML", "Unknown")
  }  
  
  y =
  geo[
    sapply( as.list( strings ), 
            function(x)
            {
              g = c()
              for(i in 1: length(geo.key) )
              {
                g[ length(g) + 1] <- ifelse( grepl( geo.key[i], x, ignore.case = TRUE), i, NA) 
              }
              
              if( length( which( !is.na(g) ) ) == 1 ){ return( which( !is.na(g) ) ) }else{ return( length(geo.key) + 1 ) }
            }
    ) ]
  
  print( strings[ which(y == "Unknown")  ]  )
  return( y )
  
  #v20180110
}


### hyphy_slac -------------------------------

hyphy_slac <- function( folder = getwd() )
{
  require(stringr)
  require(readr)
  
  # prep space 
  gene  <- c()
  clade <- c()
  year  <- c()
  geo   <- c()
  
  w          <- c()
  w.u        <- c()
  w.l        <- c()
  dn.mean    <- c() 
  ds.mean    <- c() 
  dn_ds.mean <- c()
  dn_ds.sd   <- c()
  
  po.site <- c()
  ne.site <- c()
  
  file.log <- paste0( folder, list.files( folder )[ grep( "t.log", list.files( folder ) ) ] )
  file.txt <- paste0( folder, list.files( folder )[ grep( "t.txt", list.files( folder ) ) ] )
  
  if( !length(file.log) == length(file.txt) ){ stop("inconsistent .log and .txt number") }else
  {
    
    for(i in 1: length( file.log ) )
    {
      
      str_i <- str_match( list.files( folder )[ grep( "t.log", list.files( folder ) ) ][i], "([0-9]+[a-z])_([a-zA-Z]+)_([0-9T]+)-([a-zA-Z0-9]{2,3})" )
      
      log.i <- read_file( file.log[i] )
      txt.i <- read.table( file.txt[i], sep = "\t", header = T)
      
      gene  <- c( gene, str_i[1,5] ) 
      clade <- c( clade, str_i[1,2] )
      year  <- c( year, str_i[1,4] )
      geo   <- c( geo, str_i[1,3] )
      
      w          <- c( w, as.numeric( str_match(log.i, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2] ) )
      w.u        <- c( w.u, as.numeric( str_match(log.i, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4] ) )
      w.l        <- c( w.l, as.numeric( str_match(log.i, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3] ) )
      
      dn.mean    <- c( dn.mean, mean( txt.i[,8] ) ) 
      ds.mean    <- c( ds.mean, mean( txt.i[,7] )  ) 
      dn_ds.mean <- c( dn_ds.mean, mean( txt.i[,9] ) )
      dn_ds.sd   <- c( dn_ds.sd, sd( txt.i[,9] ) )  
      
      po.site <- c( sum( txt.i[,10][ which( txt.i[,10] != 0 ) ] < 0.1 ), po.site )
      ne.site <- c( sum( txt.i[,11][ which( txt.i[,11] != 0 ) ] < 0.1 ), ne.site )
      
      
    }
  }
  
  df.slac <- data.frame( gene, clade, year, geo, w, w.u, w.l, dn.mean, ds.mean, dn_ds.mean, dn_ds.sd, po.site, ne.site, 
                         stringsAsFactors = FALSE)
  return(df.slac)
  
  #v20171023
}



### temSample --------------------------------

temSample <- function( fasfile    = file.choose(), 
                       faslist    = list(), 
                       list.x     = c("id", "time", "chr"), 
                       samplelist = list( chr1 = c("t1", "n-to-keep") ), 
                       seed       = 666 )
{
  fas.n <- attributes( read.fasta( fasfile ) )$names
  fas.s <- getSequence( read.fasta( fasfile ) )
  
  tem.m <- match( fas.n, faslist[[ list.x[1] ]] )
  
  if( TRUE %in% is.na(tem.m) ){stop()}
  
  x.y <- floor( faslist[[ list.x[2] ]][ tem.m ] )
  x.c <- faslist[[ list.x[3] ]][ tem.m ]
  
  set.seed( seed )
  tokill <- c()
  for( k in 1: length(samplelist) )
  {
    for( j in 1: ( length( samplelist[[k]] )/2 ) )
    {
      tem.i  <- which(  x.c ==  attributes(samplelist)$names[k] &  x.y == samplelist[[k]][ (j-1)*2 + 1 ] )
      tokill <- c( tokill, sample( fas.n[tem.i], ( length(tem.i) - samplelist[[k]][ j*2 ] ) ) )
    }
  }
  
  k.m <- match( tokill, fas.n )
  r.i <- seq(1, length(fas.n) )[ -k.m ]
  
  print( table( x.c[ r.i ] ) )
  
  write.fasta( fas.s[ r.i ], fas.n[ r.i ], file.out = gsub( ".fasta", ".ts.fasta", fasfile  ) )
  
  
  #v20171110b
}

### alternative way to do sequence extraction from tree --------------------------------

leafEx <- function( filedir   = file.choose(), 
                    leaflist  = c() )
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
  #v20171130
}


### site-level unique sequences --------------------------------

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
  
  write.fasta( sequences = fas.s1[remain], names = fas.i1[remain], file.out = gsub( ".fasta", "_rmdP.fasta", fasdir) )
  
  #v20171130
}

### search for ac in the faslist --------------------------------


acSearch <- function( faslist     = list(),
                      ac.dir      = 1,
                      keyword     = NULL, 
                      keyword.dir = c(9),
                      range       = NULL, 
                      range.dir   = 4 
)
{
  a.1 <- c()
  a.2 <- c()
  
  if( length( keyword.dir ) > 0 )
  {
    a.1 = seq(1, length( faslist[[ ac.dir ]] ) )
      
    for( i in 1: length( keyword.dir ) )
    {
      a.t <- which( faslist[[ keyword.dir[i] ]] == keyword[i] )
      a.1 <- intersect( a.1, a.t )
    }
  }
  
  
  if( length(range) >1 )
  {
    a.2 <- c( a.2, which( floor( faslist[[ range.dir ]] ) <= range[2] & floor( faslist[[ range.dir ]] ) >= range[1] ) )
  }
  
  if( ( length(keyword) > 0 ) & ( length(range) > 0) )
  {
    return( faslist[[ ac.dir ]][ sort( unique( intersect( a.1, a.2 ) ) ) ] ) 
    
  }else
  {
    return( faslist[[ ac.dir ]][ sort( unique( c( a.1, a.2 ) ) ) ] ) 
    
  }
  
  #v20180108
}


### python color --------------------------------  

pyCol <- function( name = c( "red", "blue", "green" ) )
{
  Col.code <- c( "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                 "#8c564b", "e3777c2", "#7f7f7f", "#bcbd22", "#17becf" )
  Col.name <- c( "blue", "orange", "green", "red", "purple",
                 "brown", "pink", "gray", "yellow", "cyan")
  
  
  return( Col.code[ match( name, Col.name ) ] )
  
  #v20180110
  
}

### jumpMx --------------------------------  

jumpMx <- function( states = c("cnC", "cnE", "cnN", "cnS", "cnSW") )
{
  lth <- length( states )
  mx  <- matrix( rep(0,lth^2), nrow = lth)
  ord <- combn( lth, 2)
  
  dirname <- c() 
  string  <- c()
  for( x in 1: ( dim(ord)[2] ) )
  {
    dirname = c( dirname, paste0( states[ ord[,x][1] ], "_to_", states[ ord[,x][2] ] ) )
    dirname = c( dirname, paste0( states[ rev(ord[,x])[1] ], "_to_", states[ rev(ord[,x])[2] ] ) )
  }
  
  for( x in 1: ( dim(ord)[2] ) )
  {
    mx1  <- matrix( rep(0,lth^2), nrow = lth)
    mx1[ ord[,x][1], ord[,x][2]  ]          <- 1
    string  <- c( string, paste( paste( as.vector( t(mx1) ), ".0", sep = ""), collapse = " ") )
    
    mx2  <- matrix( rep(0,lth^2), nrow = lth)
    mx2[ ord[,x][2], ord[,x][1] ] <- 1
    string  <- c( string, paste( paste( as.vector( t(mx2) ), ".0", sep = ""), collapse = " ") )
  }
    
  out <- c()
  for( y in 1: length(dirname) )
  {
    out <- c(out, paste("<parameter id=", dirname[y], " value=", string[y], "/>",  sep = 'b' ) )
  }
  
  print("Remember to replace b to quotation marks!!!")
  write.table( out, file = "trans.mx", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #v20180112
}


### improved findtaxa --------------------------------

findtaxa2 <- function( type       = 2, 
                       tree       = NULL, 
                       targetid   = "chiken", 
                       target     = "red",
                       Refdf      = data.frame(name = NA, state = NA),
                       useRefdf   = TRUE,
                       background = "black")
{
  
  # type 1 = branch coloring 
  # type 0 = tip shape
  # default branch color = black
  
  library(ape)
  library(ggtree)
  
  # extract tree data
  
  tree.d                         <- fortify(tree)
  tree.d[, ncol(tree.d) + 1]     <- gsub("'", "", tree.d$label)
  colnames(tree.d)[ncol(tree.d)] <- "taxaid"
  
  # for tip shape
  
  if ( useRefdf )
  {
    tem.m    <- match( gsub( "'", "", tree$tip.label ), Refdf$name )
    if( TRUE %in% is.na(tem.m) ){ stop() }
    
    taxainfo <- Refdf$state[ tem.m ]
    
  }else
  {
    taxainfo <- tree.d$taxaid
  }
  
  
  if (type == 0)
  {
    
    shapetaxa <- data.frame( node   = c( 1:length( tree.d$isTip ) ), 
                             shapee = NA)
    
    for (i in 1: length(targetid) )
    {
      shapetaxa$shapee[  grep(  tolower( targetid[i] ), tolower( taxainfo ) ) ] <- target[i]
    }
    
    return(shapetaxa)
    
  }else {
    
    # for branch colorring  
    
    # new column
    tree.d[, ncol(tree.d) + 1]     <- background
    colnames(tree.d)[ncol(tree.d)] <- "colorr"
    
    # for branch extension
    edgematrix <- as.matrix( tree.d[,c(2,1)] )
    
    # color grouping 
    group_color <- unique( target )
    
    for ( i in 1: length( group_color ) )
    {
      
      # color as group to combine key word to targetno
      sub_color <- which( target == group_color[i] )
      targetno  <- c()
      
      for ( t in 1: length( sub_color ) )
      {
        targetno <- 
          unique( c( targetno, grep( tolower( targetid[ sub_color[t] ] ), 
                                     tolower( taxainfo ) )
          ) )
      }
      
      pre_targetno  <- length( targetno )
      post_targetno = 0
      
      # while loop 
      
      while( pre_targetno != post_targetno )
      {
        
        pre_targetno = length( targetno )
        
        for( k in 1: length(targetno) )
        {
          # all sibiling 
          sibs <- as.integer( edgematrix[ 
            which( edgematrix[,1] == 
                     edgematrix[ which( edgematrix[,2] == targetno[k] ), ][1] ),][,2] )
          
          if ( length( sibs ) == 1 )
          {
            targetno = c( targetno, as.integer( edgematrix[which( edgematrix[,2] == targetno[k]), ][1] ) )
            
          }else{
            
            if ( length( which(! sibs %in% targetno ) ) == 0 )
            {
              # add parent no
              targetno = c( targetno, as.integer( edgematrix[which( edgematrix[,2] == targetno[k]),][1] ) )
            }
          }
          
          targetno = unique( targetno )
        }
        
        post_targetno = length(targetno)
        
      }
      
      # coloring
      
      tree.d$colorr[targetno] <- group_color[i]
      
    }
    return(tree.d)    
    
  }
  #v20180201
}








