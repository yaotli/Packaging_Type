# Aim of the project: 
# to understand the subtype dist. in Gs/GD and 
# inspect the essential cis-determinants
# for this section:
# 1. extract subtree seq
# 2. call back deleted seq
#
# Data source: 
# 1. gisaid: 6278/ H5
# 2. ncbi: 6677/ H5
#
# download: 20170323
# minlth: 1500
# pool: pool_ha_12955.fasta

# repeat previous work ####

  
source("_R/Function.R")

  cleanID()
  
  # source: pool_ha_12955.fasta
  # result: cleanID_pool_ha_12955.fasta
  
curateSeq(maxamb = 5, minseq = 1600, mode = 5)

  # mode 5: duplicated seq x duplicated id

  # source: cleanID_pool_ha_12955.fasta
  # result: curateSeq-5_pool_ha_12955.fasta
  # add reference file: smallref.fasta
  # result: addref_curateSeq-5_pool_ha_12955.fasta


  # MAFFT
  # source: addref_curateSeq-5_pool_ha_12955.fasta
  # result: align_addref_pool_ha_12955.fasta

  # trim: lth = 1581

curateSeq(maxamb = 1500, minseq = 1, mode = 2, vip = 238)

  # source: trim_addref_pool_ha_12955.fasta
  # result: curateSeq-2_trim_pool_ha_12955.fasta (n = 6582)


  # FastTree
  # source: curateSeq-2_trim_pool_ha_12955.fasta 
  # ./FastTree -nt -spr 4 -nni 10 -gtr -cat 20 -gamma  <tobetree.fasta> tree
  # result: 

# extract subclade seq ####

  # GsGD / non-GsGD

subtreseq(findrep = 1, outlier = 1, originfile = 1)

  # source: 1. trim_addref_pool_ha_12955.fasta, 
  #         2. id_sub5496, 
  #         3. align_addref_pool_ha_12955.fasta

  # result: 1. align_addref_GsGD.fasta, 2 align_addref_GsGDc.fasta 
  
  # trim: align_addref_GsGD.fasta to the length of 1581


  # 2344 / non-2344

subtreseq(findrep = 1, outlier = 1, originfile = 1)

  # source: 1. trim_align_addref_GsGD.fasta, 
  #         2. id_sub888,
  #         3. align_addref_GsGD.fasta
          
  # result: 1. align_addref_GsGD_2322.fasta, 2. align_addref_GsGD_2322c.fasta

  
  # nonN1 in non2344 clade

subtreseq(findrep = 1, outlier = 0, originfile = 1)

  # souce: 1. trim_align_addref_GsGD.fasta,   
  #        2. non2344_nonN1,
  #        3. align_addref_GsGD.fasta

  
# seperate N1 and nonN1 ####

  library(seqinr)
  library(stringr)

  file <- read.fasta(file.choose())

  
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
       
  # eliminate unsubtyped seq
       
  nonsubtype <- which(is.na(str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name0)[,1]) == TRUE)
  nonsubtype <- unique( sort( c(nonsubtype, grep("H5N0", x = seq_name0)) ) )
  
  if ( length(nonsubtype) > 1 ){
  
    seq_name <- seq_name0[-nonsubtype]
    seq <- seq0[-nonsubtype]       
    
  }else{
    
    seq_name = seq_name0
         seq = seq0
    
  }
  

     n1_id <- which( str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name)[,2] == "H5N1" )
  nonN1_id <- seq(1:length(seq_name))[-n1_id]
  
  
  write.fasta(seq[n1_id], 
              file.out = "~/Desktop/N1.fasta", 
              names = seq_name[n1_id])
  
  write.fasta(seq[nonN1_id], 
              file.out = "~/Desktop/Nx.fasta", 
              names = seq_name[nonN1_id])
  

# Difference ####

  # source: 5pRNA_2344.fas and 5pRNA_GsGDnon2344.fas

setwd("~/Desktop/3groups/UTR/")  
filenames <- list.files(getwd())  

library(RWebLogo)

for (k in 1: length(filenames)){

  weblogo(file.in = filenames[k], open = FALSE, 
          file.out = paste0("~/Desktop/", gsub(".fasta", "", filenames[k]), ".pdf"), 
          color.scheme = 'classic', stacks.per.line = 150, units = 'probability')

}

 


 
  
  
  