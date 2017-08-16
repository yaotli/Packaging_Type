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
  if ( is.na(daydifference) )
  {
    x   <- sub(pattern = "01", replacement = "02", x)
    
    daydifference <- as.numeric( difftime( strptime(x, format = "%Y-%m-%d"),
                                           strptime(yr.0, format = "%Y-%m-%d"), 
                                           units = "days"
                                           )
                                 )/365
  }
    
  yr.daydifference <- yr + daydifference
  yr.daydifference <- format( round( yr.daydifference, 3 ), nsmall = 3)
  
  return(yr.daydifference)
  
  #v20170614
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

trimtool <- function(propblank = 0.8, filedir = file.choose()){
  
  library(stringr)
  library(seqinr)
  
  file = read.fasta(filedir)
  seq_name0 = attributes(file)$names
  seq0 = getSequence(file)
  
  seq_matrix = do.call(rbind, seq0)
  
  coltoberemove = apply(seq_matrix, 2, function(x){
    
    blank = ( length( which( x == "-") ) + length( which( x == "~") ) ) 
    fl = length(x)
    
    if (fl*propblank < blank){
      
      return(1)
      
    }else{
      
      return(0)
    }   
  }
  )
  
  cut_matrix = seq_matrix[,-which(coltoberemove == 1)]
  
  seq_cut = as.list( data.frame(t(cut_matrix), stringsAsFactors = FALSE) )
  
  filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fas)" )[,2]
  
  write.fasta(seq_cut, 
              file.out = paste0("./trim_", filename, ".fasta"),
              names = seq_name0)
  print("DONE")

  #20170729  
}


## keepLongSeq ----------------


keepLongSeq <- function(seq_0, 
                        id_0, 
                        showRemain = FALSE)
{
  library(seqinr)
  
  if( length( which( duplicated(id_0) == TRUE) ) )
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
    
    remain <- seq(1, length(seq_0))[-toberemove]
    
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
  
  #v201706
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

subfastaSeq <- function(subtype = "H5N1", 
                        time_s  =  1000,
                        time_e  =  3000,
                        filedir = file.choose(),
                        invertedsubtype = FALSE, 
                        AC      = FALSE, 
                        AC_list = NA, 
                        no = "")
{
  require(seqinr)
  require(stringr)
  
  file  = read.fasta(filedir)
  
  seq_name0 <- attributes(file)$names
  seq0      <- getSequence(file)
  
  # subtype 
  
  if (AC == TRUE)
  {
    ac_code  <- "EPI[0-9]+|[A-Z]{1,2}[0-9]{5,6}"
    ac.file  <- str_match( seq_name0 , ac_code )
    
    ac.i     <- match( AC_list, ac.file)
    
    if ( TRUE %in% is.na(ac.i) )
    {
      print("ERROR")
      
    }else
    {
      filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)")[,2]
      
      write.fasta(sequences = seq0[ac.i], 
                  names     = seq_name0[ac.i],
                  file.out  = paste0("./", "ac", no, "_", filename, ".fasta") )
      
      print( seq_name0[ac.i] )
      
    }
    
  }else
  {
    
    if( invertedsubtype == TRUE )
    {
      subtype_i <- grep(pattern = paste0("_", subtype, "_"), 
                        x       = seq_name0, 
                        invert  = TRUE)   
    }else
    {
      subtype_i <- grep(pattern = paste0("_", subtype, "_"), 
                        x       = seq_name0)        
    }
    
    
    # time of isolation
    
    T_iso  <- as.numeric( gsub("_", "", str_match(seq_name0, "_([0-9]{4})\\.([0-9]+)")[,1] )  )
    
    time_i <- which(T_iso < time_e & T_iso > time_s)
    
    # output 
    
    remain   <- sort( intersect( subtype_i, time_i ) )
    
    filename <- str_match(filedir, "([a-zA-Z0-9_-]+)(\\.)(fasta)")[,2]
    
    write.fasta(sequences = seq0[remain], 
                names     = seq_name0[remain],
                file.out  = paste0("./", filename, "_", subtype, "-", time_s, "-", time_e, ".fasta") )
    
    print( seq_name0[remain] )
    
    
  }
  
  
  #v20170729
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





