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

  
## FastTree ----------------

# file: fasttree_1

# extract Eurasian strains (remove one Korea at the root) and build the tree
# n = 2255

# convert .tre to .nwk

### ggtree -------------------------------- 

sub2255_file <- read.tree("./tree/sub_2255")
sub2255_T0   <- ggtree(sub2255_file, size = 0.3)

serotype_all <- c("h5n1", "h5n2", "h5n3", "h5n4", "h5n5", "h5n6", "h5n7", "h5n8", "h5n9")
sero_color   <- c("black", gg_color_hue(8))
  
sub2255_color_sero <- 
  findtaxa(type = 1, sub2255_file, targetid = serotype_all, target = sero_color)


## basic ggtree with NA subtype ----------------

  
sub2255_sero <- 
  sub2255_T0 %<+%     
  sub2255_color_sero + 
  aes(color = I(colorr)) + 
  geom_tippoint(size = 0.1)

gzoom(sub2255_sero, c(1:641), widths = c(0.5,0.5))


# to inspect the node number:
# sub2255_sero + geom_text(aes(label = node), size = 0.5)
# to rotate the branch:
# ggtree::rotate(T_h5tree_note, 2766)
  
## label by tip ---------------- 
  
sub2255_color_sero[,12] = "N1"
  
for (k in 1: length(sub2255_color_sero[,12]) )
{
  sub2255_color_sero[,12][k] <- 
    sub(pattern = "h5", replacement = "", serotype_all)[
      match( sub2255_color_sero[,11][k], sero_color) ]
  
}
  
colnames(sub2255_color_sero)[12] = "NAtype"
  
# tree with annotation

sub2255_sero_tipanno <- 
  sub2255_T0 %<+% sub2255_color_sero +
  geom_tippoint(aes(color = NAtype)) + 
  scale_color_manual( values =  sero_color ) + 
  theme(legend.position = "left", legend.text = element_text(size= 30 )) +
  guides(colour=guide_legend("NA", override.aes = list(size = 20 )))


## extract sequence info ----------------

# prepare data.frame for heatmap

taxa_sub2255    <- 
  gsub( "'", "", fortify(sub2255_file)$label[which( fortify(sub2255_file)$isTip == TRUE)] )

sub2255_5p_seq  <- file_5p124_seq[ match(taxa_sub2255, file_5p124_id) ]
sub2255_seqMx   <- do.call(rbind, sub2255_5p_seq)

sub2255_heatmap <- data.frame(p24 = sub2255_seqMx[,24], 
                             p35 = sub2255_seqMx[,35], 
                             p72 = sub2255_seqMx[,72],
                             p108 = sub2255_seqMx[,108], stringsAsFactors = FALSE)

sub2255_heatmap$p108[ which(sub2255_heatmap$p108 == "y") ] = "-"

rownames(sub2255_heatmap) <- 
  paste0("'", file_5p124_id[ match(taxa_sub2255, file_5p124_id) ], "'")
  
## heatmaptree ----------------

sub2255_maptree <- 
gheatmap(sub2255_sero, 
         sub2255_heatmap, width=0.25) +
scale_fill_manual(breaks=c( "-", "a", "u", "c", "g"),
                  values=c( "white", "darkgreen", "steelblue", "orange", "firebrick") ) 



### 2.3.4 subtree --------------------------------

sub641_file <- read.tree("./tree/sub_641")
sub641_T0   <- ggtree(sub641_file, size = 0.3)

sub641_color_sero <- 
  findtaxa(type = 1, sub641_file, targetid = serotype_all, target = sero_color)

sub641_sero <- 
  sub641_T0 %<+%     
  sub641_color_sero + 
  aes(color = I(colorr)) + 
  geom_tippoint(size = 0.1)



# heatmap 

taxa_sub641    <- 
  gsub( "'", "", fortify(sub641_file)$label[which( fortify(sub641_file)$isTip == TRUE)] )

sub641_5p_seq  <- file_5p124_seq[ match(taxa_sub641, file_5p124_id) ]
sub641_seqMx   <- do.call(rbind, sub641_5p_seq)

sub641_heatmap <- data.frame(p24 = sub641_seqMx[,24], 
                              p35 = sub641_seqMx[,35], 
                              p72 = sub641_seqMx[,72],
                              p108 = sub641_seqMx[,108], stringsAsFactors = FALSE)


rownames(sub641_heatmap) <- 
  paste0("'", file_5p124_id[ match(taxa_sub641, file_5p124_id) ], "'")


sub641_maptree <- 
  gheatmap(sub641_sero, 
           sub641_heatmap, width=0.25) +
  scale_fill_manual(breaks=c( "-", "a", "u", "c", "g"),
                    values=c( "white", "darkgreen", "steelblue", "orange", "firebrick") ) 


### msaplot --------------------------------
  
  file4 <- read.fasta("~/Desktop/cuated_5pRNA_pool_ha_18100 2.fasta")    
  
  seq.name4 = attributes(file4)$names
       seq4 = getSequence(file4)
    
  seq.name4 = paste0("'", seq.name4, "'")     
       
  write.fasta(seq4, names = seq.name4, file.out = "~/Desktop/out.fasta")
  
  msaplot(T_h5tree_annotate, "~/Desktop/out.fasta", width = 0.25)
  
  
  
  
  
  
  
  
  
  
  