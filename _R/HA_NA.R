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
  
  # in Cmd line #
  # mafft ~/Desktop/HA_NA/cleanID_HAp.fasta > ~/Desktop/HA_NA/align_cleanID_HAp.fasta
  
# trim
  
  trimtool(propblank = 0.8, filedir = "~/Desktop/HA_NA/align_cleanID_HAp.fasta")
  
  # manualy checked by BioEdit (lth = 1686)
  
  curateSeq(maxamb = 15, minseq = 1500, mode = 2, 
            filedir = "~/Desktop/HA_NA/align_trim/trim_align_cleanID_HAp.fasta")
  
# FastTree
  
  # ./FastTree -nt -spr 4 -nni 10 -gtr -cat 20 -gamma <~/Desktop/HA_NA/align_trim/curateSeq-2_trim_align_cleanID_HAp.fasta > tree_curateSeq-2_trim_align_cleanID_HAp.tre
  
  
# use HA tree to fish out NA seq
  
  # inputfile: cleanID_NUp.fasta
  
  subtreseq(findrep = 0, outlier = 0, originfile = 0)
  
  curateSeq(maxamb = 13, minseq = 1350, mode = 1, 
            filedir = "~/Desktop/HA_NA/tree/GsGD_NU_4575.fasta")
  

  
# Nx splitting
  
  curateSeq(maxamb = 13, minseq = 1350, mode = 1, 
            filedir = "~/Desktop/HA_NA/cleanID_NUp.fasta")
  
  # mafft
  # mafft mafft ~/Desktop/HA_NA/curateSeq-1_cleanID_NUp.fasta > ~/Desktop/HA_NA/align_curateSeq-1_cleanID_NUp.fasta

  
  # split NA according to their subtypes
  
   align_NU <- read.fasta("~/Desktop/HA_NA/align_trim/align_curateSeq-1_cleanID_NUp.fasta")
  seq.name0 <- attributes(align_NU)$names  
       seq0 <- getSequence(align_NU)
  
     NAtype <- str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq.name0)[,2]
  
  
  for (k in 1: length(unique(NAtype)) ){
    
     indexF <- c()
     
    for(n in 1: 10){
      indexF <- c(indexF, rep(n, as.vector(table(NAtype)[n] ) ))
    }
    
    Nxindex <- split(sort(NAtype, index.return = TRUE)$ix, indexF)[[k]]
    
    write.fasta(seq0[Nxindex], names = seq.name0[Nxindex],
                
                file.out = paste0("~/Desktop/align_NU_", sort(unique(NAtype))[k], ".fasta") )
  
  }
  
  # try to auto-trim first 
  # N1 did not apply trimtool becasue of most of the seq with truncation
  # combine N0 and N1
     
  # create concensus seq from N0 - N9 
     
  lsConSeq = list()   
  for (d in 1: 9){
    
    trimN <- paste0("~/Desktop/HA_NA/Nx_splitting/trim/trim_align_NU_H5N", d, ".fasta")
    
    fileTrimN <- read.fasta(trimN)
         seq0 <- getSequence(fileTrimN)
    
    seq_matrix <- do.call(rbind, seq0)
    
    lsConSeq[[d]] =
    apply(seq_matrix, 2, function(x){
      
      most <- which.max(
        as.data.frame( 
          table(x[which( x !="~" & x !="-")]) , stringsAsFactors = FALSE)[,2])
      
         y <- as.data.frame( 
           table(x[which( x !="~" & x !="-")]) , stringsAsFactors = FALSE)[,1][most]
      
      return(y)
      
    })
    
  }
     
  write.fasta(lsConSeq, names = seq(1,9), file.out = "~/Desktop/ConSeqNx.fasta")   
     
# translatorx
  
  # result: ConSeqNx.nt_ali.fasta & ConSeqNx.aa_ali.fasta
  
   Conseq <- read.fasta("~/Desktop/HA_NA/Nx_splitting/ConSeqNx.nt_ali.fasta")
  seq_con <- getSequence(Conseq)
  
  for(z in 1:9){
  
    file <- read.fasta( paste0("~/Desktop/HA_NA/Nx_splitting/trim/trim_align_NU_H5N",z, ".fasta") )
    
         Seq0 <- getSequence(file)
    seq_name0 <- attributes(file)$names
    
    seq_matrix <- do.call(rbind, Seq0)
    
        blankv <- which(seq_con[[z]] == "-")
           ntv <- which(seq_con[[z]] != "-")
      
    mx = c()
    
    for(w in 1: length(seq_con[[z]]) ){
      
      if (w %in% blankv){
        
        mx = cbind(mx, "-")
        
      }else{
        
        mx = cbind(mx, seq_matrix[, which(ntv == w)])
        
      }
    }
    
   mx = gsub("~", "-", mx)
   mx = as.list( data.frame(t(mx), stringsAsFactors = FALSE) )
   
   write.fasta(mx, 
               file.out = paste0("~/Desktop/rearrange_H5N", z, ".fasta"),
               names = seq_name0)
      
  }
  
# fish out GsGD
  
  subtreseq()
  
  gsgd_fas <- read.fasta("~/Desktop/HA_NA/tree/GsGD_HAp_4575.fasta")
  gsgd_Seq <- getSequence(gsgd_fas)
   gsgd_ID <- attributes(gsgd_fas)$names
  
  com_NU_fas <- read.fasta("~/Desktop/HA_NA/Nx_splitting/rearrang/comb_rearrange_H5_NU")
  com_NU_Seq <- getSequence(com_NU_fas)
   com_NU_ID <- attributes(com_NU_fas)$names
  
  inHAbutnotinNA <- which( is.na(match(gsgd_ID, com_NU_ID) ) == TRUE)
  
  # write NA
  
  fishoutfromHA <- match(gsgd_ID[-inHAbutnotinNA], com_NU_ID)
  
  write.fasta(com_NU_Seq[fishoutfromHA], 
              names = com_NU_ID[fishoutfromHA], 
              file.out = "~/Desktop/GsGD_NU_p.fasta")
  

  # trim NA 
  # result: Tree_trim_GsGD_NU_p (n = 3613)
  trimtool(propblank = 0.5, filedir = "~/Desktop/GsGD_NU_p.fasta")
  
  
# prepare the cophy (tree file for tanglegram)

library(ape)    
  
  # nwk from GsGD_HAp_4575
  HAtree <- read.tree("~/Desktop/HA_NA/tanglegram/GsGD_HA_nwk")
  
  # nwk from Tree_trim_GsGD_NU_p
  NUtree <- read.tree("~/Desktop/HA_NA/tanglegram/GsGD_NU_nwk")
  
  
  NUtree$edge.length = NULL
  HAtree$edge.length = NULL
  
  hamatch <- match(NUtree$tip.label, HAtree$tip.label)
  
  NUtree$tip.label = seq(1, length(NUtree$tip.label))
  
  HAtree$tip.label[hamatch] = seq(1, length(NUtree$tip.label))
  HAtree$tip.label[-hamatch] = seq((length(NUtree$tip.label) + 1), length(HAtree$tip.label))
  
  write.tree(HAtree, file = "~/Desktop/HA_NA/tanglegram/GsGD_HA_R")
  write.tree(NUtree, file = "~/Desktop/HA_NA/tanglegram/GsGD_NU_R")
  
  # association block in NEXUS
  
  col1 <- paste0("'", seq(1, length(NUtree$tip.label)), ":")
  col2 <- paste0(seq(1, length(NUtree$tip.label)), ",")
  
  write.csv(data.frame(col1,col2, stringsAsFactors = FALSE), 
            "~/Desktop/HA_NA/tanglegram/IDlink.csv", row.names = FALSE, col.names = FALSE)
  
  
  