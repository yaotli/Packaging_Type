# Aim of the section: 
# to understand the subtype dist. in Gs/GD and 
# inspect the essential cis-determinants


# clean ID and extract info ####

library(seqinr)
library(stringr)

  file = read.fasta(file.choose())
  
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
  seq_name[which( endsWith(seq_name, "-") == "TRUE")] <- 
    gsub("-", "-15", seq_name[which( endsWith(seq_name, "-") == "TRUE")])
  

  
  # deal with duplicate ID / seq
  
  duplicated_ID <- which(duplicated(seq_name) == "TRUE")
  duplicated_note = duplicated(seq_name)
  
  
  
  


