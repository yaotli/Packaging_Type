library(seqinr)
library(stringr)

ha_file = read.fasta("~/Desktop/HA_NA/g_HA_h5_20170503.fasta")
nu_file = read.fasta("~/Desktop/HA_NA/g_NA_h5_20170503.fasta")


  HAseq_name0 = attributes(ha_file)$names
       HAseq0 = getSequence(ha_file)
     
  NUseq_name0 = attributes(nu_file)$names
  NUseq0 = getSequence(nu_file)
       
# make paired seq files ####
    
  for(n in 1:2){
  
    hana = c("HA", "NU")
    
    seq_name = get(paste0(hana[n], "seq_name0")) 
         seq = get(paste0(hana[n], "seq0")) 
    
    dupid <- which(duplicated(seq_name) == TRUE)

    toberemove = c()
        
    for (k in 1: length(dupid)){
      
        dup_i <- which(seq_name %in% seq_name[dupid[k]] == TRUE)
      
      keep_id <- which.max(sapply(seq[dup_i], length))
      
      toberemove = c(toberemove, dup_i[-keep_id])    
      
    }
    
    tobekeep = seq(1: length(seq))[-toberemove]
    
    write.fasta(seq[tobekeep], 
                names = seq_name[tobekeep], 
                file.out = paste0("~/Desktop/HA_NA/", hana[n],"p.fasta") )
      
    }
  
source("./_R/Function.R")    
  
  cleanID("~/Desktop/HA_NA/HAp.fasta")
  cleanID("~/Desktop/HA_NA/NUp.fasta")
  

# MAFFT 
  
  # try to sort and inspect the alignment
  
   align_NU <- read.fasta("~/Desktop/HA_NA/cleanID_NUp.fasta")
  seq.name0 <- attributes(align_NU)$names  
       seq0 <- getSequence(align_NU)
  
     sortid <- sort(
       str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq.name0)[,2], 
       index.return =  TRUE)$ix
  
  write.fasta(seq0[sortid], 
              names = seq.name0[sortid], 
              file.out = "~/Desktop/sorted_cleanID_NU.fasta")
  
  # manually split according to NA subtype (align_NU.fasta)
  # try to auto-trim first 
  
  for (y in 1:9){
    
    dirfas = paste0("~/Desktop/HA_NA/aligntrim/n", y, ".fas")
    
    trimtool(filedir = dirfas)  
      
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
  
      
  