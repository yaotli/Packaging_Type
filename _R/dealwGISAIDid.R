# manually replace ">_" with ">"

# read-in ####

library(seqinr)
library(stringr)
source("~/Packaging_Type/_R/Function.R")

setwd("~/Desktop/Geo/sourceSeq/")

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

# combine files

id    <- c(id_1, id_2)
seq   <- c(seq_1, seq_2)
sheet <- rbind(table_1, table_2)

# remove duplicated id ####

if (  length(which(duplicated(id) == TRUE) ) > 0  )
{
    toberemove <- c()
    dup        <- which(duplicated(id) == TRUE)
    
    print(id[dup])
  
    for(k in 1: length(dup))
    {
      
      id_dup_k <- grep(id[dup[k]], id)  
      
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
      remain     <- seq(1, length(id))[-toberemove]
      
    }
    
    id  <- id[remain]
    seq <- seq[remain]
    
}
  

# match .fasta id with .csv  ####

isolateID    <- sheet$Isolate_Id
seqIsolateID <- str_match(id, "EPI_ISL_([0-9]+)")[,1]

# _m for matched id and seq; the same order as data in sheet 

id_m         <- gsub("_EPI_ISL_([0-9]+)", "",id)[match(isolateID,seqIsolateID)]
seq_m        <- seq[match(isolateID,seqIsolateID)]


# extract Geo info from sheet ####

sheetLocation <- str_match(sheet$Location, "([A-Za-z]+) / ([A-Za-z]+)")[,3]
troubleL      <- which(is.na(str_match(sheet$Location, "([A-Za-z]+) / ([A-Za-z]+)")[,3]) == TRUE)

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

# basic id cleaning #### 

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










# deal with duplicated ####


dup_id  <- which(duplicated(id_m) == TRUE)
dup_all <- c()

for(i in 1: length(dup_id))
{
  dup_all = c(dup_all, which(match(id_m, id_m[dup_id[i]]) != "NA") )
  
}

dup_all <- unique(dup_all)

  print(id_m[dup_all])
  print(sheet$Passage_History[dup_all])
  
# keep: 
# print(sheet$Passage_History[dup_all[c(1,3,6,7,9,11,14)]])

toberemoved <- dup_all[-c(1,3,6,7,9,11,14)]








# deal with time ####





id[which(duplicated(id_rmEPI) == TRUE)]












# remove Mixed and Lab-derived sample ####

# remove 


head(id_rmEPI)
head(sheet$Isolate_Name[sheet2seq])


if ( length( grep("Mixed", id) ) > 0 )
{
  mixedid <- grep("Mixed", id)
  
  print(id[mixedid])
  
  remain <- seq(1, length(seq))[-mixedid]
  id     <- id[remain]
  seq    <- seq[remain]
  
} 


# add Geo info. (country level)




# match isolate id 

isolateID    <- sheet$Isolate_Id
seqIsolateID <- str_match(id, "EPI_ISL_([0-9]+)")[,1]

# remove isolateID

id_rmEPI <- gsub("_EPI_ISL_([0-9]+)", "",id)



id[which(startsWith(id, "A/") == FALSE)]



sheet$Isolate_Name[match(seqIsolateID, isolateID)]







