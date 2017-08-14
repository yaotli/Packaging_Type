# data preparation following seqPrep.R

library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

pool_csv          <- "./pool_df.csv"
gsgdtree          <- "./Tree/h5_GsGD"
n1tree            <- "./Tree/N1_pool"

# result files
h5_GsGD_table_dir <- "./h5_GsGD_table.csv"
n1_table_dir      <- "./n1_table.csv"   

### Clade info ----------------------------------

## HA ----------------

# tree file import and data retrive 

pool_table      <- read.csv(pool_csv, header = TRUE, stringsAsFactors = FALSE)

h5_GsGDfile     <- read.tree(gsgdtree)
h5_GsGDtreedata <- fortify(h5_GsGDfile)

rootnode        <- length(h5_GsGDfile$tip.label) + 1

# info for the dataframe

h5_GsGD_R2tip <- dist.nodes(h5_GsGDfile)[rootnode, 1: (rootnode - 1)]

h5_GsGD_name  <- gsub("'", "", h5_GsGDtreedata$label[1: rootnode - 1] )
h5_GsGD_time  <- as.numeric( str_match( h5_GsGD_name, "[0-9]{4}\\.[0-9]{2,3}") )
h5_GsGD_NA    <- str_match( h5_GsGD_name, "(_H5)(N[0-9])")[,3]

ac_code       <- "EPI[0-9]+|[A-Z]{1,2}[0-9]{5,6}"
h5_GsGD_ac    <- str_match( h5_GsGD_name, ac_code )

h5_GsGD_geo   <- pool_table$geo[ match(h5_GsGD_ac, pool_table$ac) ] 

h5_GsGD_table <- data.frame(name    = h5_GsGD_name,  
                            subtype = h5_GsGD_NA,
                            ac      = h5_GsGD_ac,
                            geo     = h5_GsGD_geo,
                            distance = h5_GsGD_R2tip, 
                            time     = h5_GsGD_time, 
                            stringsAsFactors = FALSE)


# clade info
# manually prepare txt file for clade annotation

for(i in 1: length( list.files("./Clade/clade_txt/") ) )
{
  cladename_i   <- read.table(file = 
                                paste0("./Clade/clade_txt/", list.files("./Clade/clade_txt/")[i]),
                              stringsAsFactors = FALSE)
  
  cladematch    <- as.numeric( h5_GsGD_name %in% cladename_i[,1] )
  
  h5_GsGD_table[, ncol(h5_GsGD_table) + 1] =  cladematch
  
  colnames( h5_GsGD_table )[ ncol( h5_GsGD_table ) ] =
    sub(".txt", "", list.files("./Clade/clade_txt/")[i] )
}


## NA ----------------

n1_file     <- read.tree(n1tree)
n1_treedata <- fortify(n1_file)

rootnode_nu <- length(n1_file$tip.label) + 1

# info for the dataframe

n1_R2tip <- dist.nodes(n1_file)[rootnode_nu, 1: (rootnode_nu - 1)]

n1_name  <- gsub("'", "", n1_treedata$label[1: rootnode_nu - 1] )
n1_time  <- as.numeric( str_match( n1_name, "[0-9]{4}\\.[0-9]{2,3}") )

ac_code  <- "EPI[0-9]+|[A-Z]{1,2}[0-9]{5,6}"
n1_ac    <- str_match( n1_name, ac_code)

n1_geo   <- pool_table$geo[ match(n1_ac, pool_table$ac) ] 

n1_table <- data.frame( name     = n1_name,  
                        ac       = n1_ac,
                        geo      = n1_geo,
                        distance = n1_R2tip, 
                        time     = n1_time, 
                        stringsAsFactors = FALSE)

## match NA to HA ----------------

for(k in 7: ncol(h5_GsGD_table) )
{
  
  n1_table_i <- 
    h5_GsGD_table[, k][ match( pool_table$pac[ match( n1_table$ac, pool_table$ac ) ], 
                               h5_GsGD_table$ac ) ]
  
  n1_table_i[ is.na(n1_table_i) ] = 0
  
  n1_table[, ncol(n1_table) + 1] = n1_table_i
  colnames( n1_table )[ ncol(n1_table) ] = colnames(h5_GsGD_table)[k]
  
}


# export 
# write.csv(h5_GsGD_table, "h5_GsGD_table.csv")
# write.csv(n1_table, "n1_table.csv")



### HA and N1 from 3 clades in China --------------------------------

# for further clade-wide sampleing
# geo: China (CN) plus Hong_kong (HK)

# HA 

for(i in 7: 9)
{
  index = which( ( h5_GsGD_table$geo == "China" | h5_GsGD_table$geo == "Hong_Kong" ) &  
                   h5_GsGD_table$subtype == "N1" & h5_GsGD_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_h5_pool_lth1683.fasta", 
              AC      = TRUE,
              AC_list = h5_GsGD_table$ac[ index ], 
              no      = colnames(h5_GsGD_table)[i] )
}

  
# NA

for(i in 6: 8)
{
  index = which( ( n1_table$geo == "China" | n1_table$geo == "Hong_Kong" ) & 
                       n1_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_nu_pool_H5N1_lth1347.fasta", 
              AC      = TRUE,
              AC_list = n1_table$ac[ index ], 
              no      = colnames(n1_table)[i] )
}


# fasttree for 2 genes each clade
#
# cd ~/Desktop/data_souce/Clade/clade_CNHK_fas/
# for f in $(ls *.fasta);
# do ~/fasttree -nt -spr 4 -nni -gtr -cat 20 -gamma -notop <$f> ~/Desktop/data_souce/Clade/sampled_clade_tree_CNHK/$f.tre ;
# done


### updated clade assignment by tree files --------------------------------

# read-in
# h5_GsGD_table <- read.csv( h5_GsGD_table_dir, header = TRUE, stringsAsFactors = FALSE)
# n1_table      <- read.csv( n1_table_dir, header = TRUE, stringsAsFactors = FALSE)


dir_sampledtre  <- "./Clade/sampled_clade_tree_CNHK/sampled_tre/"
list_sampledtre <- list.files( dir_sampledtre )

for( k in 1: length(list_sampledtre) )
{
  tre_k <- read.csv( paste0(dir_sampledtre, list_sampledtre[k]) , stringsAsFactors = FALSE)[, 1]
  
  ntax  <- as.numeric( str_match(tre_k [ grep("ntax", tre_k) ], "(ntax=)([0-9]+)")[,3] )
  tax_s <- grep("ntax", tre_k) + 2
  tax_e <- tax_s + ntax - 1
  
  taxaname <- tre_k[tax_s: tax_e]
  reded    <- grep("#ff0000", taxaname)
  remain   <- taxaname[ - reded ]
  remain   <- gsub("\t|\\'", "", remain)
  remain   <- gsub("\\[\\&\\!color\\=#[A-Za-z0-9]{6}\\]", "", remain)
  
  write.table(remain, 
              paste0("./Clade/sampled_clade_tree_CNHK/sampled_txt/", 
                     gsub("\\.|fasta|tre", "", list_sampledtre[ k ]) ), 
              col.names = F, row.names = F, quote = F)
  
  print( length(remain) )
  
}


# HA

for(i in 1: length( grep("_h5", list.files("./Clade/updated/") ) ) )
{
  
  k = grep("_h5", list.files("./Clade/updated/") )[ i ]
  
  cladename_i   <- read.table( file = 
                                 paste0("./Clade/updated/", list.files("./Clade/updated/")[ k ]),
                               stringsAsFactors = FALSE)
  
  cladematch    <- as.numeric( h5_GsGD_table$name %in% cladename_i[,1] )
  
  updatname     <- list.files("./Clade/updated/")[ k ]
  
  if( updatname %in% colnames(h5_GsGD_table) == TRUE )
  {
    j = which( colnames(h5_GsGD_table) %in% updatname == TRUE )
    
    h5_GsGD_table[ ,j] = cladematch
    
  }else
  {
    h5_GsGD_table[, ncol(h5_GsGD_table) + 1] =  cladematch
    
    colnames( h5_GsGD_table )[ ncol( h5_GsGD_table ) ] = updatname
    }
  
}

# NA 

for(i in 1: length( grep("_n1", list.files("./Clade/updated/") ) ) )
{
  
  k = grep("_n1", list.files("./Clade/updated/") )[ i ]
  
  cladename_i   <- read.table( file = 
                                 paste0("./Clade/updated/", list.files("./Clade/updated/")[ k ]),
                               stringsAsFactors = FALSE)
  
  cladematch    <- as.numeric( n1_table$name %in% cladename_i[,1] )
  
  updatname     <- list.files("./Clade/updated/")[ k ]
  
  if( updatname %in% colnames(n1_table) == TRUE )
  {
    j = which( colnames(n1_table) %in% updatname == TRUE )
    
    n1_table[ ,j] = cladematch
    
  }else
  {
    n1_table[, ncol(n1_table) + 1] =  cladematch
    
    colnames( n1_table )[ ncol( n1_table ) ] = updatname
  }
  
}


## Export table --------------------------------

write.csv(h5_GsGD_table, "h5_GsGD_table.csv")
write.csv(n1_table, "n1_table.csv")

### Extract sampled and filtered seq --------------------------------

# HA 
for(i in 11: ncol(h5_GsGD_table) )
{
  index = which(h5_GsGD_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_h5_pool_lth1683.fasta", 
              AC      = TRUE,
              AC_list = h5_GsGD_table$ac[ index ], 
              no      = colnames(h5_GsGD_table)[i] )
  
}


# NA 
for(i in 10: ncol(n1_table) )
{
  index = which(n1_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_nu_pool_H5N1_lth1347.fasta", 
              AC      = TRUE,
              AC_list = n1_table$ac[ index ], 
              no      = colnames(n1_table)[i] )
  
}

# batch rename 
#
# cd ~/Desktop/data_souce/Clade/sampled_clade_tree_CNHK/sampled_fas/
# for f in *.fasta; do mv "$f" "${f/ac/}"; done
# for f in *.fasta; do mv "$f" "${f/_trim_h5_pool_lth1683/}"; done
# for f in *.fasta; do mv "$f" "${f/_trim_nu_pool_H5N1_lth1347/}"; done

### Sampled seq of clade 232 before 2014 --------------------------------

s_clads232 <- c("./Clade/sampled_clade_tree_CNHK/sampled_fas/s_clade232_h5_CNHK.fasta", 
                "./Clade/sampled_clade_tree_CNHK/sampled_fas/s_clade232_n1_CNHK.fasta")


for(k in 1: length(s_clads232) )
{
  subfastaSeq( time_e = 2014,
               filedir = s_clads232[k] )
}


### stratification by time (2yr) --------------------------------

fas_sample <- paste0("./Clade/sampled_clade_tree_CNHK/sampled_fas/", 
                     list.files("./Clade/sampled_clade_tree_CNHK/sampled_fas/") )

for(i in 1: length(fas_sample))
{
  # per 2 year  
  
  subfastaSeq( time_e = 2006,
               filedir = fas_sample[i] )
  
  subfastaSeq( time_s = 2006, time_e = 2008,
               filedir = fas_sample[i] )
  
  subfastaSeq( time_s = 2008, time_e = 2010,
               filedir = fas_sample[i] )
  
  subfastaSeq( time_s = 2010,
               filedir = fas_sample[i] )
}


# trim name 
# for f in *.fasta; do sed -i "" "s/\\./_/g" $f; done
# for f in *.fasta; do sed -i "" "s/~/-/g" $f; done
# for f in *.fasta; do mv "$f" "${f/s_/}"; done
# for f in *.fasta; do mv "$f" "${f/_CNHK_H5N1/}"; done


## HA1/ HA2 (hyphy) ----------------

fas_sample <- paste0("~/Desktop/data_souce/seq_2yr_CNHK/full/", 
                     list.files("~/Desktop/data_souce/seq_2yr_CNHK/full/") )

fas_sample_HA <- fas_sample[ grep(pattern = "h5", x = fas_sample) ] 

for(k in 1: length(fas_sample_HA) )
{
  ntpartition( position = c(1:1017), filedir = fas_sample_HA[k], no = "-ha1.fasta")
  ntpartition( position = c(1018:1683), filedir = fas_sample_HA[k], no = "-ha2.fasta")
  
}


## NAh/ NAs (hyphy) ----------------

fas_sample <- paste0("~/Desktop/data_souce/seq_2yr_CNHK/full/", 
                     list.files("~/Desktop/data_souce/seq_2yr_CNHK/full/") )

fas_sample_NA <- fas_sample[ grep(pattern = "n1", x = fas_sample) ] 

for(k in 1: length(fas_sample_NA) )
{
  ntpartition( position = c(211:1347), filedir = fas_sample_NA[k], no = "-nah.fasta")
  ntpartition( position = c(1:210), filedir = fas_sample_NA[k], no = "-nas.fasta")
  
}


### stratification by time (1yr) --------------------------------

fas_sample <- paste0("./Clade/sampled_clade_tree_CNHK/sampled_fas/", 
                     list.files("./Clade/sampled_clade_tree_CNHK/sampled_fas/") )

for(i in 1: length(fas_sample))
{
  # per 1 year  
  yr_seq <- seq(2002, 2014)
  
  for (k in 1: (length(yr_seq)-1) )
  {
    subfastaSeq( time_s = yr_seq[k],
                 time_e = yr_seq[k+1],
                 filedir = fas_sample[i] )
  }
}


# trim name 
# for f in *.fasta; do sed -i "" "s/\\./_/g" $f; done
# for f in *.fasta; do sed -i "" "s/~/-/g" $f; done
# for f in *.fasta; do mv "$f" "${f/s_/}"; done
# for f in *.fasta; do mv "$f" "${f/_CNHK_H5N1/}"; done


## remove ambiguous nt -----------------

fas_sample <- paste0("./pool_seq_pyr_CNHK/", 
                     list.files("./pool_seq_pyr_CNHK/") )

for(k in  1: length(fas_sample) ){
  rmAMBnt(fas_sample[k])  }



## HA1/ HA2 ----------------

fas_sample <- paste0("./pool_seq_pyr_CNHK/", 
                     list.files("./pool_seq_pyr_CNHK/") )

fas_sample_HA <- fas_sample[ grep(pattern = "h5", x = fas_sample) ] 

for(k in 1: length(fas_sample_HA) )
{
  ntpartition( position = c(1:1017), filedir = fas_sample_HA[k], no = "-ha1.fasta")
  ntpartition( position = c(1018:1683), filedir = fas_sample_HA[k], no = "-ha2.fasta")
  
}


## NAh/ NAs ----------------


fas_sample_NA <- fas_sample[ grep(pattern = "n1", x = fas_sample) ] 

for(k in 1: length(fas_sample_NA) )
{
  ntpartition( position = c(211:1347), filedir = fas_sample_NA[k], no = "-nah.fasta")
  ntpartition( position = c(1:210), filedir = fas_sample_NA[k], no = "-nas.fasta")
  
}




