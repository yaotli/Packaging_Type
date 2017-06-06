# The aim of this section is to examine the RNA 5p nucleotide distribution on tree

source("~/Packaging_Type/_R/Function.R")
setwd("~/Desktop/Geo/")

library(seqinr)
library(stringr)
library(ggtree)
library(ggplot2)


### curate and align --------------------------------

# curate the sequence file

curateSeq(maxamb = 150, minseq = 1, mode = 2, filedir = "./H5_merged_9197.fasta")

## mafft ----------------

# file: mafft_1
  

### extract seq with intact 5p end --------------------------------
  
file_5p124     <- read.fasta("./align_cut/5p124nt_H5_merged_8602.fasta")
file_5p124_id  <- attributes(file_5p124)$names
file_5p124_seq <- getSequence(file_5p124)
  

## remove no subtype and p24/p35 ----------------

nonsubtype <- 
  which(is.na(str_match(pattern = "_(H[0-9]+N[0-9]+)_", file_5p124_id)[,1]) == TRUE)

nonsubtype <- unique( sort( c(nonsubtype, grep("H5N0", x = file_5p124_id)) ) )

# read seq into matrix

seqMx_5p124      <- do.call(rbind, file_5p124_seq)
both2rDash       <- intersect( which(seqMx_5p124[,24] == "-"), which(seqMx_5p124[,35] == "-") )
toberemove_5p124 <- unique( sort( c(nonsubtype, both2rDash) ) )
remain_5p124     <- seq(1, length(file_5p124_seq))[-toberemove_5p124]     
  
# output file

alined_8602 <- read.fasta("./align_cut/align_H5_merged_8602.fasta")    

write.fasta(getSequence(alined_8602)[remain_5p124], 
            attributes(alined_8602)$names[remain_5p124], 
            file.out = "./align_cut/5p124_align_H5_merged_2769.fasta")  


### build the ML tree --------------------------------  
  
curateSeq(maxamb = 150, minseq = 1500, mode = 2, 
            filedir = "./align_cut/5p124_align_H5_merged_2769.fasta")

# trim by trimtool plus manually modify to 1680 nt

trimtool(propblank = 0.99, 
         filedir = "./align_cut/curateSeq-2_5p124_align_H5_merged_2744.fasta")

# last curation to remove identical sequences

curateSeq(maxamb = 150, minseq = 10, mode = 2, 
          filedir = "./align_cut/trim_curateSeq-2_5p124_align_H5_merged_2744.fasta")


## FastTree ----------------

# file: fasttree_1

# extract Eurasian strains (remove one Korea/2004 at the root) and build the tree
# n = 2255

# convert .tre to .nwk

### ggtree -------------------------------- 

sub2068_file <- read.tree("./tree/sub_2068")
sub2068_T0   <- ggtree(sub2068_file, size = 0.3)

serotype_all <- c("h5n1", "h5n2", "h5n3", "h5n4", "h5n5", "h5n6", "h5n7", "h5n8", "h5n9")
sero_color   <- c("black", gg_color_hue(8))
  
sub2068_color_sero <- 
  findtaxa(type = 1, sub2068_file, targetid = serotype_all, target = sero_color)


## basic ggtree with NA subtype ----------------

  
sub2068_sero <- 
  sub2068_T0 %<+%     
  sub2068_color_sero + 
  aes(color = I(colorr)) + 
  geom_tippoint(size = 0.1)

gzoom(sub2068_sero, c(1:571), widths = c(0.5,0.5))


# to inspect the node number:
# sub2068_sero + geom_text(aes(label = node), size = 0.5)
# to rotate the branch:
# ggtree::rotate(T_h5tree_note, 2766)
  
## label by tip ---------------- 
  
sub2068_color_sero[,12] = "N1"
  
for (k in 1: length(sub2068_color_sero[,12]) )
{
  sub2068_color_sero[,12][k] <- 
    sub(pattern = "h5", replacement = "", serotype_all)[
      match( sub2068_color_sero[,11][k], sero_color) ]
  
}
  
colnames(sub2068_color_sero)[12] = "NAtype"
  
# tree with annotation

sub2068_sero_tipanno <- 
  sub2068_T0 %<+% sub2068_color_sero +
  geom_tippoint(aes(color = NAtype)) + 
  scale_color_manual( values =  sero_color ) + 
  theme(legend.position = "left", legend.text = element_text(size= 30 )) +
  guides(colour=guide_legend("NA", override.aes = list(size = 20 )))


## extract sequence info ----------------

# prepare data.frame for heatmap

taxa_sub2068    <- 
  gsub( "'", "", fortify(sub2068_file)$label[which( fortify(sub2068_file)$isTip == TRUE)] )

sub2068_5p_seq  <- file_5p124_seq[ match(taxa_sub2068, file_5p124_id) ]
sub2068_seqMx   <- do.call(rbind, sub2068_5p_seq)

sub2068_heatmap <- data.frame(p24 = sub2068_seqMx[,24], 
                              p35 = sub2068_seqMx[,35], 
                              p72 = sub2068_seqMx[,72],
                              p108 = sub2068_seqMx[,108], stringsAsFactors = FALSE)

sub2068_heatmap$p108[ which(sub2068_heatmap$p108 == "y") ] = "-"

rownames(sub2068_heatmap) <- 
  paste0("'", file_5p124_id[ match(taxa_sub2068, file_5p124_id) ], "'")
  
## heatmaptree ----------------

sub2068_maptree <- 
gheatmap(sub2068_sero, 
         sub2068_heatmap, width=0.25) +
scale_fill_manual(breaks=c( "-", "a", "u", "c", "g"),
                  values=c( "white", "darkgreen", "steelblue", "orange", "firebrick") ) 



### 2.3.4 subtree --------------------------------

sub571_file <- read.tree("./tree/sub_571")
sub571_T0   <- ggtree(sub571_file, size = 0.4)

sub571_color_sero <- 
  findtaxa(type = 1, sub571_file, targetid = serotype_all, target = sero_color)

sub571_sero <- 
  sub571_T0 %<+%     
  sub571_color_sero + 
  aes(color = I(colorr)) 

# highlight N1

sub571_tip_N1 <- 
  findtaxa(type = 0, sub571_file, targetid = "N1", target = "1")

sub571_sero_NA <- 
  sub571_sero %<+% 
  sub571_tip_N1 + 
  geom_tippoint(aes(shape = shapee), size = 1, color = "black")


# heatmap 

taxa_sub571    <- 
  gsub( "'", "", fortify(sub571_file)$label[which( fortify(sub571_file)$isTip == TRUE)] )

sub571_5p_seq  <- file_5p124_seq[ match(taxa_sub571, file_5p124_id) ]
sub571_seqMx   <- do.call(rbind, sub571_5p_seq)

sub571_heatmap <- data.frame(p24 = sub571_seqMx[,24], 
                              p35 = sub571_seqMx[,35], 
                              p72 = sub571_seqMx[,72],
                              p108 = sub571_seqMx[,108], stringsAsFactors = FALSE)

rownames(sub571_heatmap) <- 
  paste0("'", file_5p124_id[ match(taxa_sub571, file_5p124_id) ], "'")


sub571_maptree <- 
  gheatmap(sub571_sero_NA, 
           sub571_heatmap, width=0.25) +
  scale_fill_manual(breaks=c( "-", "a", "u", "c", "g"),
                    values=c( "white", "darkgreen", "steelblue", "orange", "firebrick") ) 


### msaplot --------------------------------
  
  file4 <- read.fasta("~/Desktop/cuated_5pRNA_pool_ha_18100 2.fasta")    
  
  seq.name4 = attributes(file4)$names
       seq4 = getSequence(file4)
    
  seq.name4 = paste0("'", seq.name4, "'")     
       
  write.fasta(seq4, names = seq.name4, file.out = "~/Desktop/out.fasta")
  
  msaplot(T_h5tree_annotate, "~/Desktop/out.fasta", width = 0.25)