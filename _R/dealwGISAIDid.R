# manually replace ">_" with ">" in files from GISAID
# manually replace blank with "_" in files from NCBI

### read-in --------------------------------

library(seqinr)
library(stringr)
source("~/Packaging_Type/_R/Function.R")

setwd("~/Desktop/Geo/sourceSeq/")

# ncbi

fasta_n = "./H5_ncbi_7024_20170601.fasta"

data_n <- read.fasta(fasta_n)
id_n   <- attributes(data_n)$names
seq_n  <- getSequence(data_n)



# gisaid

fasta_1 = "./H5_gisaid1_1181_20170602.fasta"
fasta_2 = "./H5_gisaid2_1317_20170602.fasta"
csv_1   = "./H5_gisaid1_1181_20170602.csv"
csv_2   = "./H5_gisaid2_1315_20170602.csv"

data_1  <- read.fasta(fasta_1)
id_1    <- attributes(data_1)$names
seq_1   <- getSequence(data_1)
table_1 <- read.csv(csv_1, header = TRUE, stringsAsFactors = FALSE)

data_2  <- read.fasta(fasta_2)
id_2    <- attributes(data_2)$names
seq_2   <- getSequence(data_2)
table_2 <- read.csv(csv_2,header = TRUE, stringsAsFactors = FALSE)


# combine files from gisaid 

id    <- c(id_1, id_2)
seq   <- c(seq_1, seq_2)
sheet <- rbind(table_1, table_2)


### remove duplicated id (same isolateID) -------------------------------- 

# keep the one with longest sequence

if (  length(which(duplicated(id) == TRUE) ) > 0  )
{
    toberemove <- c()
    dup        <- which(duplicated(id) == TRUE)
    
    print(id[dup])
  
    for(k in 1: length(dup) )
    {
      
      id_dup_k <- which( id %in% id[ dup[k] ] == TRUE)
        
      seqL     <- which.max(
        
        sapply(seq[id_dup_k], function(x)
        {
          
          y = c2s(x)
          z = gsub("-", "", y)
          z = gsub("~", "", z)
          
          l = length(s2c(z))
          
          return(l)
        }
        )

      )
      
      toberemove <- c(toberemove, id_dup_k[-seqL])
      
    }
    
    remain <- seq(1, length(id))[-toberemove]
    id     <- id[remain]
    seq    <- seq[remain]
    
}
  

### match .fasta id with .csv --------------------------------

isolateID    <- sheet$Isolate_Id
seqIsolateID <- str_match(id, "EPI_ISL_([0-9]+)")[,1]


# _m for matched id and seq; the same order as data in sheet

id_m  <- gsub("_EPI_ISL_([0-9]+)", "",id)[match(isolateID,seqIsolateID)]
seq_m <- seq[match(isolateID,seqIsolateID)]


### basic id cleaning --------------------------------

    id_m = gsub(" ", "_", id_m)
    id_m = gsub("\\(|\\)|\\[|\\[|\\.", "-", id_m)
    id_m = gsub("\\'|\\?|\\>", "", id_m)
    id_m = gsub("_A_/", "", id_m)
    
# deal with date
    
    id_m = gsub("_-Month_and_day_unknown-", "-07-01", id_m)
    id_m = gsub("_-Day_unknown-", "-15", id_m)

# generate temporal info
    
    d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
    
    for (i in 1: length(id_m))
    {
      id_m[i] <- gsub(x = id_m[i], pattern = d, replacement = 
                     
                     phylo_date(
                       
                       str_match( string = id_m[i], pattern = d)[,1]
                       
                     ))
    }


# check " TRUE %in% endsWith(id_m, "NA") "


### add Geo info. (country level) --------------------------------
    
# extract from sheet
    
rawLocation   <- gsub(pattern = " ", "_", sheet$Location)
sheetLocation <- gsub("_$", "", 
                      str_match(rawLocation, "([A-Za-z_]+)_/_([A-Za-z_]+)")[,3])

troubleL      <- which(is.na(sheetLocation) == TRUE)

if( length(troubleL) > 0 )
{
  for (i in 1: length(troubleL))
  {
    if( !is.na( str_match(sheet$Location[troubleL[i]], "[A-Za-z]+")[1] ) )
    {
      
      sheetLocation[ troubleL[i] ] = str_match(sheet$Location[troubleL[i]], "[A-Za-z]+")[1]
      
    }else
    {
      sheetLocation[ troubleL[i] ]  = "NA"
    } 
  }
  
  if ( TRUE %in% is.na(troubleL) )
  {
    print("ERROR") 
  }
}
    

# insert geo info  
# _g stands for id with geo info
    
dd   = "([0-9]{4})\\.([0-9]{2})"
id_g = c()

for(i in 1: length(seq_m))
{
  id_g[i] <- sub(pattern = str_match(id_m[i], dd)[,1],
                 replacement = paste0("|", sheetLocation[i], "|_", str_match(id_m[i], dd)[,1]),
                 x = id_m[i]
                 )
}

    
### deal with identical ID and mixed strains --------------------------------

## identical ID ----------------

# keep the one with longest sequence

if ( length( which(duplicated(id_g) == TRUE) ) > 0 )
{
  toberemove_2 <- c()
  dup_id       <- which(duplicated(id_g) == TRUE)
  
  for(i in 1: length(dup_id) )
  {
    
    dup_id_i <- which( id_g %in% id_g[dup_all[i]] == TRUE)
    
    SeqL     <- which.max(
      sapply(seq_m[dup_id_i], function(x)
      {
        
        y = c2s(x)
        z = gsub("-", "", y)
        z = gsub("~", "", z)
        
        l = length(s2c(z))
        
        return(l)
      }
      )
    )
    
    toberemove_2 <- c(toberemove_2, dup_id_i[-SeqL] ) 
    
  }
}


## mixed sample ----------------

if ( length( grep("Mixed", id_g) ) > 0 )
{
  mixedid <- grep("Mixed", id_g)
  print(id_g[mixedid])

} 

toberemove_2 <- sort( unique(c(toberemove_2, mixedid)) )

# remove and make final files with _gisaid

remain_2   <- seq(1, length(seq_m))[-toberemove_2]
id_gisaid  <- id_g[remain_2]
seq_gisaid <- seq_m[remain_2]


### prepare ncbi data and merge two database --------------------------------


## basic ID cleaing ----------------

    id_n_ed = gsub(" ", "_", id_n)
    id_n_ed = gsub("\\(|\\)|\\[|\\[|\\.", "-", id_n_ed)
    id_n_ed = gsub("\\'|\\?|\\>", "", id_n_ed)
    id_n_ed = gsub("_A_/", "", id_n_ed)

# deal with date

    id_n_ed[which( endsWith(id_n_ed, "_--") == "TRUE")] <- 
      gsub("_--", "_9999-99-99", id_n_ed[which( endsWith(id_n_ed, "_--") == "TRUE")])
    
    id_n_ed[which( endsWith(id_n_ed, "--") == "TRUE")] <- 
      gsub("--", "-07-01", id_n_ed[which( endsWith(id_n_ed, "--") == "TRUE")])
    
    id_n_ed[which( endsWith(id_n_ed, "-") == "TRUE" )] <- 
      paste0(id_n_ed[which( endsWith(id_n_ed, "-") == "TRUE" )], "15")
    
    
# generate temporal info

    d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
    
    for (i in 1: length(id_n_ed))
    {
      id_n_ed[i] <- gsub(x = id_n_ed[i], pattern = d, replacement = 
                        
                        phylo_date(
                          
                          str_match( string = id_n_ed[i], pattern = d)[,1]
                          
                        ))
    }

## deal with idnetical id --------------------------------


if (  length(which(duplicated(id_n_ed) == TRUE) ) > 0  )
{
  toberemove <- c()
  dup_id_n   <- which(duplicated(id_n_ed) == TRUE)
  
  for( k in 1: length(dup_id_n) )
  {
    
    id_dup_k <- which( id_n_ed %in% id_n_ed[ dup_id_n[k] ])
    
    seqL     <- which.max(
      
      sapply(seq_n[id_dup_k], function(x)
      {
        
        y = c2s(x)
        z = gsub("-", "", y)
        z = gsub("~", "", z)
        
        l = length(s2c(z))
        
        return(l)
      }
      )
    )
    print(id_n_ed[id_dup_k])
    
    toberemove <- c(toberemove, id_dup_k[-seqL])
    
  }
  
  remain  <- seq(1, length(seq_n))[-toberemove]
  id_n_ed <- id_n_ed[remain]
  seq_n   <- seq_n[remain]
  
}


## combine ncbi and gisaid ----------------   
   
# n = 9246
    
 id_final <- c(id_gisaid, id_n_ed)
seq_final <- c(seq_gisaid, seq_n)

# sort based on time  

    yeari <- as.numeric( gsub("_", "", str_match(id_final, "_([0-9]{4})\\.([0-9]+)")[,1]) )
sortindex <- sort(yeari, na.last = FALSE, index.return = TRUE)$ix

write.fasta(seq_final[sortindex], 
            id_final[sortindex], 
            file.out = "~/Desktop/Geo/H5_merged_9243.fasta")  
  
  
  
    

