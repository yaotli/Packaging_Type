# The aim of this section is to examine the nucleotide change distribution on tree
#
# 

  # clean and curate

  source("./_R/Function.R")

  cleanID("~/Desktop/pool_ha_18100.fasta")

  curateSeq(maxamb = 1000, minseq = 1, mode = 5, 
            filedir = "~/Desktop/cleanID_pool_ha_18100.fasta")  
  
  # seq number = 11239
  
  
  # mafft: 
  # source: cuated_5pRNA_pool_ha_18100.fasta
  # result: align_curateSeq-5_pool_ha_18100.fasta
  
  # reverse-complement the aligned file and keep the 5 prime region
  # source: align_curateSeq-5_pool_ha_18100.fasta
  # : 5pRNA_pool_ha_18100.fasta
  
  
# extract seq with intact 5 prime end ####
  
  
  library(seqinr)
  library(stringr)
  
  
  file <- read.fasta("~/Desktop/gheatmap/5pRNA_pool_ha_18100.fasta")
  seq.name0 = attributes(file)$names
       seq0 = getSequence(file)
  
  # no subtype info     
  nonsubtype <- which(is.na(str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq.name0)[,1]) == TRUE)
  nonsubtype <- unique( sort( c(nonsubtype, grep("H5N0", x = seq.name0)) ) )
       
  seq_matrix <- do.call(rbind, seq0)
  
  # no both positions we are interested
        nont <- intersect(which(seq_matrix[,24] == "-"), which(seq_matrix[,35] == "-"))
  
  toberemove <- unique(sort(c(nonsubtype, nont)))
     
    toremain <- seq(1, length(seq0))[-toberemove]     
  
  # remove no subtype info and no position seq at 24 and 35       
    
  file2 <- read.fasta("~/Desktop/align_curateSeq-5_pool_ha_18100.fasta")    
  seq.name2 = attributes(file2)$names
       seq2 = getSequence(file2)
  
         
  write.fasta(seq2[toremain], names = seq.name2[toremain], 
              file.out = "~/Desktop/partial_align_pool_ha_18100.fasta")  
  
  # n = 3452
  
  # curate the file for tree
  
  curateSeq(maxamb = 15, minseq = 1500, mode = 7, 
            filedir = "~/Desktop/partial_align_pool_ha_18100.fasta")

  
  # trim by BioEdit
  # results: trim_curateSeq-7_partial_pool_ha_18100
  
  
  # FastTree
  # source: trim_curateSeq-7_partial_pool_ha_18100
  # result: tree_trim_curateSeq-7_partial_align_pool_ha_18100
  
  
  
# prepare the tree
  
  library(ggtree)
  library(ape)
  
    h5tree <- read.tree("~/Desktop/gheatmap/nwtree_trim_curateSeq-7_partial_align_pool_ha_18100.tre")
  T_h5tree <-  ggtree(h5tree, size = 0.3)
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n] }
  
  
  serotype <- c("H5N1", "H5N2", "H5N3", "H5N4", "H5N5", 
                "H5N6", "H5N7", "H5N8", "H5N9")
  
  sero_color <- c("black", gg_color_hue(8))
  
  sero_color_h5 <- findtaxa(type = 1, h5tree, targetid = serotype, target = sero_color)

  
  # tree with color on branch
  
  T_h5tree_note = T_h5tree %<+% sero_color_h5 + aes(color = I(colorr)) + 
    geom_tippoint(size = 0.1)
  
  # to inspect the node number:
  # T_h5tree_note + geom_text(aes(label = node), size = 0.5)
  
  T_h5tree_note_r = ggtree::rotate(T_h5tree_note, 2766)
  T_h5tree_note_r = ggtree::rotate(T_h5tree_note_r, 2767)
  T_h5tree_note_r = ggtree::rotate(T_h5tree_note_r, 2992)

  
  sero_color_h5[,12] = "N1"
  
  for (k in 1: length(sero_color_h5[,12])){
    
    sero_color_h5[,12][k] = 
      sub(pattern = "H5", replacement = "", serotype)[
        match(sero_color_h5[,11][k], sero_color)]
  }
  
  colnames(sero_color_h5)[12] = "NAtype"
  
  # tree with annotation
  
  T_h5tree_annotate = T_h5tree %<+% sero_color_h5 + 
    geom_tippoint(aes(color = NAtype ), size = 0.8) + 
    scale_color_manual(values=sero_color) + 
    theme(legend.position = "left", legend.text = element_text(size=30)) +
    guides(colour=guide_legend("NA", override.aes = list(size =20)))
  
# extract info from .fasta and make the matrix heatmap ####
  
  file3 <- read.fasta("~/Desktop/gheatmap//trim_curateSeq-7_partial_align_pool_ha_18100.fasta")    
  seq.name3 = attributes(file3)$names
  seq3 = getSequence(file3)
  
  # seq.name0 from 5pRNA_pool_ha_18100.fasta
  treematch <- match(seq.name3, seq.name0)
  
  write.fasta(seq0[treematch], names = seq.name0[treematch], 
              file.out = "~/Desktop/curated_5pRNA_pool_ha_18100.fasta")
  
  file5 <- read.fasta("~/Desktop/gheatmap//curated_5pRNA_pool_ha_18100.fasta")    
  seq.name5 = attributes(file5)$names
       seq5 = getSequence(file5)
  
  treeid <-        
  match(gsub("'", "", T_h5tree$data$label[which(T_h5tree$data$isTip == TRUE)]), seq.name5)
  
  seq_matrix2 = do.call(rbind, seq5[treeid])   
       
  tips = sero_color_h5$taxaid
  tips = tips[which(is.na(tips ) == FALSE)]
  
   rna_24_5p <- c()
   rna_35_5p <- c()
   rna_72_5p <- c()
  rna_108_5p <- c()

  for(k in 1: length(tips)){
    
    rna_24_5p[length(rna_24_5p) + 1] = seq_matrix2[,24][match(tips[k], seq.name3[treeid] )] 
    rna_35_5p[length(rna_35_5p) + 1] = seq_matrix2[,35][match(tips[k], seq.name3[treeid] )]
    rna_72_5p[length(rna_72_5p) + 1] = seq_matrix2[,72][match(tips[k], seq.name3[treeid] )] 
    rna_108_5p[length(rna_108_5p) + 1] = seq_matrix2[,108][match(tips[k], seq.name3[treeid] )] 
    
    
  }
  
  rna_108_5p = gsub(pattern = "y", "-", rna_108_5p)
  
  rna_matrix = data.frame(p24 = rna_24_5p, p35 = rna_35_5p, 
                          p72 = rna_72_5p, p108 = rna_108_5p,
                          stringsAsFactors = FALSE)  
  
  rownames(rna_matrix) = sero_color_h5$label[ which(sero_color_h5$isTip == TRUE) ]
  
  
  gheatmap(T_h5tree_note_r, rna_matrix, width=0.25) +
    scale_fill_manual(breaks=c( "-", "a", "u", "c", "g"), 
                      values=c( "white", "darkgreen", "steelblue", "orange", "firebrick") ) 


# msaplot 
  
  file4 <- read.fasta("~/Desktop/cuated_5pRNA_pool_ha_18100 2.fasta")    
  
  seq.name4 = attributes(file4)$names
       seq4 = getSequence(file4)
    
  seq.name4 = paste0("'", seq.name4, "'")     
       
  write.fasta(seq4, names = seq.name4, file.out = "~/Desktop/out.fasta")
  
  msaplot(T_h5tree_annotate, "~/Desktop/out.fasta", width = 0.25)
  
  
  
  
  
  
  
  
  
  
  