# data preparation: seqPrep.R

library(seqinr)
library(stringr)
library(ggtree)
library(ape)
library(dplyr)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")

pool_csv <- "./pool_df.csv"
gsgdtree <- "./Tree/h5_GsGD"
n1tree   <- "./Tree/N1_pool"


### HA data ----------------------------------

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
h5_GsGD_ac    <- str_match( h5_GsGD_name, ac_code)

h5_GsGD_geo   <- pool_table$geo[ match(h5_GsGD_ac, pool_table$ac) ] 

h5_GsGD_table <- data.frame(name    = h5_GsGD_name,  
                            subtype = h5_GsGD_NA,
                            ac      = h5_GsGD_ac,
                            geo     = h5_GsGD_geo,
                            distance = h5_GsGD_R2tip, 
                            time     = h5_GsGD_time, 
                            stringsAsFactors = FALSE)


## clade info ----------------

for(i in 1: length( list.files("./Clade/cladetxt/") ) )
{
  cladename_i   <- read.table(file = 
                                paste0("./Clade/cladetxt/", list.files("./Clade/cladetxt/")[i]),
                                stringsAsFactors = FALSE)
  
  cladematch    <- as.numeric( h5_GsGD_name %in% cladename_i[,1] )
  
  h5_GsGD_table[, ncol(h5_GsGD_table) + 1] =  cladematch
  
  colnames( h5_GsGD_table )[ ncol( h5_GsGD_table ) ] =
    sub(".txt", "", list.files("./Clade/cladetxt/")[i] )
}



### NA data ----------------------------------

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

# subtree and sampling 

# write.csv(h5_GsGD_table, "h5_GsGD_table.csv")
# write.csv(n1_table, "n1_table.csv")



### UPDATE THE TABLE USING TREE FILES --------------------------------

# make labeld tree file into txt

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
              paste0("./Clade/sampled_clade_tree_CNHK/sampledtxt/", list_sampledtre[ k ], ".txt"), 
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
  
  h5_GsGD_table[, ncol(h5_GsGD_table) + 1] =  cladematch
  
  colnames( h5_GsGD_table )[ ncol( h5_GsGD_table ) ] =
    sub(".txt", "", list.files("./Clade/updated/")[ k ] )
}

# NA 

for(i in 1: length( grep("_n1", list.files("./Clade/updated/") ) ) )
{
  
  k = grep("_n1", list.files("./Clade/updated/") )[ i ]
  
  cladename_i   <- read.table( file = 
                               paste0("./Clade/updated/", list.files("./Clade/updated/")[ k ]),
                               stringsAsFactors = FALSE)
  
  cladematch    <- as.numeric( n1_table$name %in% cladename_i[,1] )
  
  n1_table[, ncol(n1_table) + 1] =  cladematch
  
  colnames( n1_table )[ ncol( n1_table ) ] =
    sub(".txt", "", list.files("./Clade/updated/")[k] )
}




### Export table --------------------------------

write.csv(h5_GsGD_table, "h5_GsGD_table.csv")
write.csv(n1_table, "n1_table.csv")

# stratified by geo

h5_GsGD_table_cnhk <- h5_GsGD_table[which(h5_GsGD_table$geo == "China" | h5_GsGD_table$geo == "Hong_Kong" ), ]
n1_table_cnhk      <- n1_table[which(n1_table$geo == "China" | n1_table$geo == "Hong_Kong" ), ]



### Figuree --------------------------------

ggplot( data = h5_GsGD_table, 
        aes(x = time, y = distance) ) + 
  geom_point(color = "gray") +
  geom_point( data = h5_GsGD_table[which(h5_GsGD_table$subtype == "N1"), ],
              aes(x = time, y = distance) ) +
  theme_bw() + 
  facet_wrap(~ geo)


## clade234_HA_China ----------------

clade234_HA_cnhk <- 
  
  ggplot( data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade234_h5_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$clade234 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$s_clade234_h5_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  geom_vline(xintercept = 2008, size = 1, color = "black") + 
  ggtitle("clade234_HA_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  xlab("") + ylab("Root-to-tip divergence")



## clade7_HA_China ----------------

clade7_HA_cnhk <- 
  
  ggplot( data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade7_h5_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$clade7 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$s_clade7_h5_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  geom_vline(xintercept = 2006, size = 1, color = "black") + 
  ggtitle("clade7_HA_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  xlab("") + ylab("Root-to-tip divergence")



## clade232_HA_China ----------------

clade232_HA_cnhk <- 
  
  ggplot( data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$s_clade232_h5_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_cnhk[ which(h5_GsGD_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$clade232 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = h5_GsGD_table_cnhk[which(h5_GsGD_table_cnhk$s_clade232_h5_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  geom_vline(xintercept = 2012, size = 1, color = "black") + 
  ggtitle("clade232_HA_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  xlab("") + ylab("Root-to-tip divergence")



## clade234_N1_China ----------------

clade234_N1_cnhk <- 
  
  ggplot( data = n1_table_cnhk[ which(n1_table_cnhk$s_clade234_n1_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = n1_table_cnhk[ which(n1_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$clade234 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$s_clade234_n1_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  geom_vline(xintercept = 2008, size = 1, color = "black") + 
  ggtitle("clade234_N1_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  scale_y_continuous(limits = c(0.05, 0.2) ) + 
  xlab("") + ylab("Root-to-tip divergence")





## clade7_NA_China ----------------

clade7_N1_cnhk <- 
  
  ggplot( data = n1_table_cnhk[ which(n1_table_cnhk$s_clade7_n1_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = n1_table_cnhk[ which(n1_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$clade7 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$s_clade7_n1_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  geom_vline(xintercept = 2006, size = 1, color = "black") + 
  ggtitle("clade7_N1_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  scale_y_continuous(limits = c(0.05, 0.2) ) + 
  xlab("") + ylab("Root-to-tip divergence")



## clade232_NA_China ----------------

clade232_N1_cnhk <- 
  
  ggplot( data = n1_table_cnhk[ which(n1_table_cnhk$s_clade232_n1_CNHK == 1), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) + 
  
  # gray background: all GsGD
  geom_point(data = n1_table_cnhk[ which(n1_table_cnhk$sGsGD == 1), ], 
             aes(x = time, y = distance), 
             color = "gray", alpha = 0.2, stroke = 0, size = 2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$clade232 == 1), ], 
             aes(x = time, y = distance), 
             shape = 1, color = "red", size = 2) + 
  
  geom_point(data = n1_table_cnhk[which(n1_table_cnhk$s_clade232_n1_CNHK == 1), ], 
             aes(x = time, y = distance), 
             shape = 21, color = "black", fill = "red", size = 2, stroke = 0.5 ,alpha = 0.8) +
  
  stat_smooth(method = "glm", colour = "#008B00", fill = "#008B00" ) + 
  geom_vline(xintercept = 2012, size = 1, color = "black") + 
  ggtitle("clade232_N1_CN+HK") + 
  
  scale_x_continuous(limits = c(2000, 2016), breaks = seq(2000, 2016, by = 2), labels = seq(2000, 2016, by = 2) ) +
  scale_y_continuous(limits = c(0.05, 0.2) ) + 
  xlab("") + ylab("Root-to-tip divergence")


## combine ----------------

multiplot(clade234_HA_cnhk, clade7_HA_cnhk, clade232_HA_cnhk, 
          clade234_N1_cnhk, clade7_N1_cnhk, clade232_N1_cnhk, ncol = 2)

