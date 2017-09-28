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
    
    shapetaxa <- data.frame(node = c(1:length(tree.d$isTip)), shapee = NA)
    
    for (i in 1: length(targetid))
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
    
    edgematrix <- as.matrix(tree.d[,c(2,1)])
    
    # color grouping 
    
    group_color <- unique(target)
    
    for (i in 1: length(group_color) )
      {
      
      # color as group to combine key word to targetno
      
      sub_color <- which(target == group_color[i] )
      targetno  <- c()
      
      for (t in 1: length(sub_color) )
        {
        
        targetno <- 
          unique( c(targetno, grep( tolower( targetid[ sub_color[t] ] ), tolower(tree.d$taxaid) )) )
      
        }
      
      tobecolor     <- c()
      pre_targetno  <- length(targetno)
      post_targetno = 0
      
      # while loop 
      
      while( pre_targetno != post_targetno )
        {
        
        pre_targetno = length(targetno)
        
        for(k in 1:length(targetno))
          {
          
          # all sibiling 
          sibs <- edgematrix[
            which(edgematrix[,1] == 
                    edgematrix[which(edgematrix[,2] == targetno[k]),][1]),][,2]
          
          if (length(sibs) == 1)
            {
            
            targetno = c(targetno, edgematrix[which(edgematrix[,2] == targetno[k]),][1])
            
            }else{
            
            if (length(which(sibs %in% targetno == "FALSE")) == 0){
              
              tobecolor = c(edgematrix[which(edgematrix[,2] == targetno[k]),][1], tobecolor)
              targetno  = c(targetno, edgematrix[which(edgematrix[,2] == targetno[k]),][1])
            }
            
          }
          targetno  = unique(targetno)
          tobecolor = unique(c(targetno, tobecolor))
        
          }
        
        post_targetno = length(targetno)
        
        }
      
      # coloring
      
      tree.d$colorr[tobecolor] <- group_color[i]
      
      }
    return(tree.d)    
    
  }
  #v201706
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
          
          y = c2s(x)
          z = gsub("-", "", y)
          z = gsub("~", "", z)
          
          l = length( s2c(z) )
          
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
  
  #v20170907
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
    "Shandong|Jiangsu|Huadong|Eastern_China|Fujian|Anhui|Shanghai|Zhejiang|Jiangxi|Nanchang"
  
  C    <- 
    "Hunan|Hubei|Henan"
  
  S    <- 
    "Hong_Kong|Shantou|Guangdong|GD|ST|Shenzhen|Guangzhou"
  
  SW   <- 
    "Guizhou|Guangxi|Yunnan|Guiyang|Tibet|Sichuan|Chongqing"
  
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
    ac_code  <- "EPI[0-9]+|[A-Z]{1,2}[0-9]{5,6}"
    ac.file  <- str_match( seq_name0 , ac_code )
    
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
  
  
  #v20170926
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
    id.s[ which(id.s == "H5|H5Nx|H5NX")  ] = "H5N0"
    
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
    id.s[ which(id.s == "H5|H5Nx|H5NX")  ] = "H5N0"
    
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
  
  #v20170920b
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
                           y = c2s(x)
                           z = gsub("-|~", "", y)
                           z = grep( pattern = "a|t|c|g", 
                                     x = y, 
                                     ignore.case = TRUE, value = TRUE )
                           
                           l = length( s2c(z) )
                           
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
  
  #v20170920b
}


### time-extraction from sequence id  --------------------------------


seqDate <- function( rawdata )
{
  library(stringr)
  
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
  
  daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d"),
                                         strptime( yr.0, "%Y-%m-%d"), 
                                         units = "days") 
  )/365
  
  # bug?
  if ( TRUE %in% is.na(daydifference) )
  {
    
    rawdata.2[ which(is.na(daydifference)) ] <- 
      sub("01$", "02", rawdata.2[ which(is.na(daydifference)) ] )
    
    
    daydifference <- as.numeric( difftime( strptime( rawdata.2, "%Y-%m-%d"),
                                           strptime( yr.0, "%Y-%m-%d"), 
                                           units = "days") 
    )/365
  }
  
  yr.daydifference <- yr + daydifference
  yr.daydifference <- format( round( yr.daydifference, 3 ), nsmall = 3)
  
  return(yr.daydifference)
  
  #v20170921b
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
  anno.tre <- read.csv( filedir, stringsAsFactors = FALSE)
  taxa.s   <- grep( x = anno.tre[,1], pattern = "taxlabels" ) + 1
  ntax     <- 
    as.numeric( str_match( grep( x       = anno.tre[,1], 
                                 pattern = "ntax", value = TRUE) , "(ntax=)([0-9]+)")[,3] 
                )
  
  taxa.e <- taxa.s + ntax - 1
  
  GsGDlike_name <- str_match( string = anno.tre[, 1][taxa.s: taxa.e], 
                              pattern = "\'([0-9A-Za-z_\\|.]+)\'" )[,2]
  
  GsGDlike_tag  <- str_match( string = anno.tre[, 1][taxa.s: taxa.e], 
                              pattern = "color=#([a-z0-9]{6})")[, 2]
  
  return(df = data.frame(id  = GsGDlike_name, 
                         tag = GsGDlike_tag, stringsAsFactors = FALSE))
  
}



### remove duplicate and filter seq based on arranged name --------------------------------

rmDup <- function( fasfile = file.choose(), 
                   year    = c(1000,3000),
                   geo     = c(),
                   sero    = "",
                   rmdup   = TRUE)
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
  
  s     <- grep( sero, str_match( id, "_(H5N[0-9]{1,2})_" )[,2] )
  
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
  #v20170927v
}

### extract id info --------------------------------

taxaInfo <- function( file    = file.choose(), 
                      useTree = FALSE, 
                      makecsv = FALSE )
{
  # input: 
  # 1 colored .tre file
  # 2 .fas file with clean id 
  
  require(seqinr)
  require(stringr)
  
  if ( useTree )
  {
    anno.tre <- read.csv( file, stringsAsFactors = FALSE)
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
    
    
    if( TRUE %in% is.na(unlist( ls[c(1:6)] ) ) ){ Stop() }
    
    
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
                       y = c2s(x)
                       z = gsub("-|~", "", y)
                       z = grep( pattern = "a|t|c|g", 
                                 x = y, 
                                 ignore.case = TRUE, value = TRUE )
                       l = length( s2c(z) )
                       
                       return(l)
                     } )
    
    ls <- list( id.a, id.g, id.s, id.y, id.n, id, seq.l)
    df <- data.frame( id.a, id.g, id.s, id.y, id.n, id, seq.l )
    
    if( TRUE %in% is.na(unlist( ls[c(1:6)] ) ) ){ Stop() }
    
  }
  
  if ( makecsv ){ write.csv(df, file = sub( ".fasta", "_info.csv", file) , row.names = FALSE) }
  
  return( ls )
  
  #v20170928v
}
