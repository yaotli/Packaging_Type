# Aim of the project: 
# to understand the subtype dist. in Gs/GD and 
# inspect the essential cis-determinants
# for this section:
# 1. clean the sequence ID
# 2. try to remove apparent duplicated sequences
#
# Data sources:
# 1. gisaid: 6278/ H5
# 2. ncbi: 6677/ H5
# 3. reference strains with clade label from Dr. Gavin Smith


# clean ID ####

library(seqinr)
library(stringr)

  
  file = read.fasta("~/Desktop/ReassortSubtype/pool_ha_13193.fasta")
  
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)

  seq_name = gsub(" ", "_", seq_name0)
  seq_name = gsub("\\(", "-", seq_name)
  seq_name = gsub("\\)", "-", seq_name)
  seq_name = gsub("\\[", "-", seq_name)
  seq_name = gsub("\\]", "-", seq_name)
  seq_name = gsub("\\'", "", seq_name)
  seq_name = gsub("\\.", "-", seq_name)  
  
  seq_name = gsub(">", "", seq_name)  
  
  # for gisaid
  seq_name = gsub("_A_/", "", seq_name)

  
  # deal with date
  
  # for gisaid
  seq_name = gsub("_-Month_and_day_unknown-", "-07-01", seq_name)
  seq_name = gsub("_-Day_unknown-", "-15", seq_name)
  
  # for ncbi
  
  # w/o all
  seq_name[which( endsWith(seq_name, "_--") == "TRUE")] <- 
    gsub("_--", "-9999-99-99", seq_name[which( endsWith(seq_name, "_--") == "TRUE")])
  
  # w/o month and day
  seq_name[which( endsWith(seq_name, "--") == "TRUE")] <- 
    gsub("--", "-07-01", seq_name[which( endsWith(seq_name, "--") == "TRUE")])
  
  # w/o day
  
  sub_seq <- which( endsWith(seq_name, "-") == "TRUE" )[ 
    !which( endsWith(seq_name, "-") == "TRUE" ) %in% seq(12956, 13193) ]
   
  seq_name[sub_seq] <- paste0(seq_name[sub_seq], "15")
    
# deal with duplicated ID ####
  
 # D (n = 4354)
  
  duplicated_ID <- which(duplicated(seq_name) == "TRUE")
  
 # D2 (n = 5003)
  
  n <- "A/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)/([A-Za-z0-9-_]+)"
  
  duplicated_ID2 <- which(duplicated(str_match(seq_name, n)[,1]) == TRUE)
  
 # S (n = 5687)
  
 duplicated_seq <-  which(
   
   duplicated(sapply(seq0, function(x){
     
     y = c2s(x)
     z = gsub("-", "", y)
  
   return(z) 
     }  
   )) == "TRUE")

 
 # D2 & S (n = 4461)
 
   duplicated_rm <- intersect(duplicated_ID2, duplicated_seq)
 
 # to be edit: (n = 205)
 
 duplicated_edit <- duplicated_ID[!duplicated_ID %in% duplicated_rm]
  
 
 # loop for each replicated case
 
 if( length(duplicated_edit) > 0 ){
   
  # for all the duplicated id
   
   for (i in 1: length(duplicated_edit)){
     
     dup0 = which(match(seq_name, seq_name[duplicated_edit[i]]) != "NA")
     
     # but not in remove group
     
     dup0 = dup0[!dup0 %in% duplicated_rm]
     
     # create a null vector to appendex  
     
     app = c("", paste0(letters, "_"), toupper(paste0(letters, "_") ) )
     
     ap.id = seq_name[dup0]
     app.id = c()
     
     # loop to deal with multiple replicated
     # find the date info, insert labeling in the middle
     
     d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
     
     time_rep <- str_match(ap.id, d)[,1]
     
     for (k in 1: length(ap.id)){
       
       app.id[k] = sub(time_rep[k], paste0(app[k], time_rep[k]), ap.id[k])
     }
     
     # back to seq_name0  
     
     seq_name[dup0] = app.id
     
   }
 }     
  
 
# check and output ####
 
 which(duplicated(seq_name[seq(1, length(seq_name))[-duplicated_rm]]) == TRUE)
 

# sequence curation ####
 
 # for ambiguous nucleotide and seq length
 
 tobedelect1 <- c()

 for(i in 1: length(seq0)){
   
   ATCG = c("a", "t", "c", "g")
   
   # convert character to string
   
   seq_0 = c2s(seq0[[i]])
   seq_i = gsub("-", "", seq_0)
   seq_i = gsub("~", "", seq_i)
   
   # seq length and number of ambiguous nucleotide
   
   seqlth = length( s2c(seq_i) )
      amb = length( which(! s2c(seq_i) %in% ATCG ) )
   
   if( ( seqlth < 1500 ) | ( amb > 5 ) ){
     
     tobedelect1[length(tobedelect1) + 1 ] = i
     
   }
   
 }
 
 # sum duplicated_rm & tobedelect1 PLUS A/unknown/unknown/G444/unknown
 tobedelect <- sort (unique( c(1, duplicated_rm, tobedelect1) ))
 
 # generate temporal info
 
 d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"

 for (i in 1: (length(seq0) - 238)) {
   
   seq_name[i] <- gsub(x = seq_name[i], pattern = d, replacement = 
                         phylo_date( str_match(string = seq_name[i], pattern = d)[,1] ))
   
 }
   
 
 # first output for alignment
 
 remain <- seq(1: length(seq0))[-tobedelect]
 
 seq_name_out <- seq_name[remain]
      seq_out <- seq0[remain]
      
 write.fasta(seq_out, 
             file.out = "~/Desktop/outputSeq1.fasta",
             names = seq_name_out)
 
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
   
  
