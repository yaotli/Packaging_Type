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
# minlth: 1600
# pool: pool_ha_12955.fasta

# repeat previous work 

  # source: pool_ha_12955.fasta

source("_R/Function.R")

  cleanID()

  # curateSeq()
  # source: cleanID_pool_ha_12955.fasta
  # mode 5: duplicated seq x duplicated id
  
curateSeq(maxamb = 5, minseq = 1600, mode = 5)

  # result: curateSeq-5_pool_ha_12955.fasta
  # add reference

  # MAFFT
  # source: addref_curateSeq-5_pool_ha_12955
  # result: align_addref_pool_ha_12955.txt

  # trim: lth = 1686

curateSeq(maxamb = 1500, minseq = 1, mode = 2, vip = 238)

  # FastTree
  # source: curateSeq-2_trim_pool_ha_12955.fasta (n = 6582)
  # result: tree.tre

# extract subclade seq


subtreseq(rmrep = 1, outlier = 0)

  # source: 1. addref_curateSeq-5_pool_ha_12955.fasta 2. id_sub5496
  # result: Seq_sub5496_addrep.fasta

subtreseq(rmrep = 1, outlier = 1)

  # source: 1. Seq_sub5496_addrep.fasta 2. id_sub888
  # result: 1. Seq_sub888_addrep.fasta 2. Seq_sub888c_addref.fasta



# Difference ####

library(RWebLogo)

  weblogo(file.in = file.choose(), color.scheme = 'classic', stacks.per.line = 150)

  # source: 5_non2344.fasta and 3_non2344.fasta

library(seqinr)

  non2344 <- read.fasta(file.choose())
  
  seq_name0 = attributes(non2344)$names
       seq0 = getSequence(non2344)

  
  non2344_nonN1 <- c( grep("_H5N2_", seq_name0), 
                      grep("_H5N5_", seq_name0), 
                      grep("_H5N6_", seq_name0), 
                      grep("_H5N9_", seq_name0) )

  write.fasta(seq0[non2344_nonN1], file.out = "~/Desktop/non2344_nonN1.fasta", 
              names = seq_name0[non2344_nonN1])
  


