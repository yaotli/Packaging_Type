# aim: the inspect N1 (avian) UTR 
#

source("_R/Function.R")

  cleanID()
  
  # source: pool_n1_avian_10821. fasta
  # result: cleanID_pool_n1_avian_10821.fasta
  
  
  curateSeq(maxamb = 5, minseq = 1300, mode = 5)
  
  # source: cleanID_pool_n1_avian_10821.fasta
  # result: curateSeq-5_cleanID_pool_n1_avian_10821.fasta
  
  
  # MAFFT
  
  # source: curateSeq-5_cleanID_pool_n1_avian_10821.fasta
  # result: align_curateSeq-5_pool_n1_avian.fasta
  
  
# split the file accorging to time
  
  library(seqinr)
  
  file <- read.fasta(file.choose())
  
  seq_name0 = attributes(file)$names
       seq0 = getSequence(file)
  
  year <- as.numeric( 
    gsub("_", "", str_match(pattern = "_([0-9]{4})\\.([0-9]+)", seq_name0)[,1]) 
    )
       
  y1 <- which(year < 2008)
  
   seq_name = seq_name0[y1]
  seq_name2 = seq_name0[ seq(1, length(seq_name0))[-y1] ]
  
   seq = seq0[y1]
  seq2 = seq0[ seq(1, length(seq_name0))[-y1] ]
  
  
  write.fasta(seq, 
              file.out = "~/Desktop/outputY1.fasta", 
              names = seq_name)
  
  write.fasta(seq2, 
              file.out = "~/Desktop/outputY2.fasta", 
              names = seq_name2)
  
# Logo
  
  library(RWebLogo)
  
  weblogo(file.in = file.choose(), 
          color.scheme = 'classic', stacks.per.line = 100, units = 'probability')
  
  
  
       
       
       
       
  
  
  
  
  
  
  
  
  
  
  
  