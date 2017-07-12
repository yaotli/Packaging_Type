
setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

pool_table    <- read.csv("./pool_df.csv", header = TRUE, stringsAsFactors = FALSE)
h5_GsGD_table <- read.csv("./h5_GsGD_table.csv", header = TRUE, stringsAsFactors = FALSE)
n1_table      <- read.csv("./n1_table.csv", header = TRUE, stringsAsFactors = FALSE)


### HA and N1 from 3 clades in China --------------------------------

# HA 
for(i in 8: 10)
{
  index = which(h5_GsGD_table$geo == "China" &  h5_GsGD_table$subtype == "N1" & h5_GsGD_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_h5_pool_lth1683.fasta", 
              AC      = TRUE,
              AC_list = h5_GsGD_table$ac[ index ], 
              no      = colnames(h5_GsGD_table)[i] )
  
}

# NA
for(i in 7: 9)
{
  index = which(n1_table$geo == "China" & n1_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_nu_pool_H5N1-1000-3000_lth1350.fasta", 
              AC      = TRUE,
              AC_list = n1_table$ac[ index ], 
              no      = colnames(n1_table)[i] )
  
}


### Sampled and filtered seq --------------------------------

# HA 
for(i in 12: ncol(h5_GsGD_table) )
{
  index = which(h5_GsGD_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_h5_pool_lth1683.fasta", 
              AC      = TRUE,
              AC_list = h5_GsGD_table$ac[ index ], 
              no      = colnames(h5_GsGD_table)[i] )
  
}


# NA 
for(i in 11: ncol(n1_table) )
{
  index = which(n1_table[, i] == 1 )
  
  subfastaSeq(filedir = "./align_trim/trim_nu_pool_H5N1-1000-3000_lth1350.fasta", 
              AC      = TRUE,
              AC_list = n1_table$ac[ index ], 
              no      = colnames(n1_table)[i] )
  
}

### Sampled seq of clade 232 before 2014 --------------------------------

s_clads232 <- c("./Clade/sampled_clade_tree_CH/sampled_fas/s_clade232_h5_CH.fasta", 
                "./Clade/sampled_clade_tree_CH/sampled_fas/s_clade232_n1_CH.fasta")


for(k in 1: length(s_clads232) )
{
  subfastaSeq( time_e = 2014,
               filedir = s_clads232[k] )
}










