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
# 3. reference strains with clade label from Dr. Gavin Smith (238)


# Prepare for tree ####


  # clean ID: pool_ha_12995 (gisaid + ncbi)

cleanID()

  # curateSeq

curateSeq(mode = 5)

  # add reference strains
  # align by MAFFT, and trim in BioEdit

  # curateSeq

curateSeq(maxamb = 1500, minseq = 1, mode = 2, vip = 238)


  # FastTree
  # ./FastTree -nt -spr 4 -nni 10 -gtr -cat 20 -gamma <tobetree.fasta> out


