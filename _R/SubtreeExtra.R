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

  # cleanID
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

  # source: 1. addref_curateSeq-5_pool_ha_12955.fasta 2. id_sub5496

subtreseq(rmrep = 1, outlier = 0)

  # result: Seq_sub888_addrep.fasta


  # source: 1. Seq_sub888_addrep.fasta 2. id_sub888

subtreseq(rmrep = 1, outlier = 1)

  # result: 1. Seq_sub888_addrep.fasta 2. Seq_sub888c_addref.fasta

















