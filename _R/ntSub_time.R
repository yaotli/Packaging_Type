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

### HA substitution over time ----------------------------------

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

# stratified by country 

h5_GsGD_table_CH <- h5_GsGD_table[which(h5_GsGD_table$geo == "China"), ]


## clade234_HA_China ----------------

clade234_HA_China <- 
  
  ggplot( data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade234 == 1 & h5_GsGD_table_CH$subtype == "N1"), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$sGsGD == 1), ], 
             aes(x = time, y = distance), color = "gray", alpha = 0.2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade234 == 1), ], 
             aes(x = time, y = distance), shape = 21, color = "#FF7256", fill = NA, stroke = 1) + 

  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade234 == 1 & h5_GsGD_table_CH$subtype == "N1"), ], 
             aes(x = time, y = distance), color = "red") +

  stat_smooth(method = "auto") + 

  geom_vline(xintercept = 2008, size = 2, color = "black") + 
  ggtitle("clade234_HA_China") + 
  xlab("")



## clade7_HA_China ----------------

clade7_HA_China <- 
  
  ggplot( data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade7 == 1 & h5_GsGD_table_CH$subtype == "N1"), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$sGsGD == 1), ], 
             aes(x = time, y = distance), color = "gray", alpha = 0.2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade7 == 1), ], 
             aes(x = time, y = distance), shape = 21, color = "#7FFF00", fill = NA, stroke = 1) + 
  
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade7 == 1 & h5_GsGD_table_CH$subtype == "N1"), ], 
             aes(x = time, y = distance), color = "green") + 
    
  stat_smooth(method = "auto") + 
    
  geom_vline(xintercept = 2006, size = 2, color = "black") +
  ggtitle("clade7_HA_China") + 
  xlab("")


## clade232_HA_China ----------------

clade232_HA_China <- 
  
  ggplot( data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade232 == 1 & h5_GsGD_table_CH$subtype == "N1"), ], 
          aes(x = time, y = distance) ) + 
  theme_bw() + 
  
  # gray background: all GsGD
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$sGsGD == 1), ], 
             aes(x = time, y = distance), color = "gray", alpha = 0.2) +
  
  # open circle: non-N1; circle: N1
  
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade232 == 1), ], 
             aes(x = time, y = distance), shape = 21, color = "#00BFFF", fill = NA, stroke = 1) + 
  
  geom_point(data = h5_GsGD_table_CH[which(h5_GsGD_table_CH$clade232 == 1 & h5_GsGD_table_CH$subtype == "N1"), ], 
             aes(x = time, y = distance), color = "blue") + 
  
  stat_smooth(method = "auto") + 
  geom_vline(xintercept = 2012, size = 2, color = "black") +
  
  ggtitle("clade232_HA_China") + 
  xlab("")



## combine ----------------

multiplot(clade234_HA_China, clade7_HA_China, clade232_HA_China, ncol = 1)



### NA substitution over time ----------------------------------

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


# stratified by country 

n1_table_CH <- n1_table[which(n1_table$geo == "China"), ]


## clade234_N1_China ----------------

clade234_N1_China <- 
  
  ggplot( data = n1_table_CH[which(n1_table_CH$clade234 == 1), ], 
          aes(x = time, y = distance) ) + 
  
  theme_bw() + 
  
  # gray background: all GsGD
  geom_point( data = n1_table_CH[which(n1_table_CH$sGsGD == 1), ], 
              aes(x = time, y = distance ), 
              color = "gray", 
              alpha = 0.2) +
  
  geom_point(data = n1_table_CH[which(n1_table_CH$clade234 == 1), ], 
             aes(x = time, y = distance), color = "red") +
  
  stat_smooth(method = "auto") + 
  
  geom_vline(xintercept = 2008, size = 2, color = "black") + 
  ggtitle("clade234_N1_China") + 
  xlab("")



## clade7_HA_China ----------------

clade7_N1_China <- 
  
  ggplot( data = n1_table_CH[which(n1_table_CH$clade7 == 1), ], 
          aes(x = time, y = distance) ) + 
  
  theme_bw() + 
  
  # gray background: all GsGD
  geom_point( data = n1_table_CH[which(n1_table_CH$sGsGD == 1), ], 
              aes(x = time, y = distance ), 
              color = "gray", 
              alpha = 0.2) +
  
  geom_point(data = n1_table_CH[which(n1_table_CH$clade7 == 1), ], 
             aes(x = time, y = distance), color = "green") +
  
  stat_smooth(method = "auto") + 
  
  geom_vline(xintercept = 2006, size = 2, color = "black") + 
  ggtitle("clade7_N1_China") + 
  xlab("")


## clade232_HA_China ----------------

clade232_N1_China <- 
  
  ggplot( data = n1_table_CH[which(n1_table_CH$clade232 == 1), ], 
          aes(x = time, y = distance) ) + 
  
  theme_bw() + 
  
  # gray background: all GsGD
  geom_point( data = n1_table_CH[which(n1_table_CH$sGsGD == 1), ], 
              aes(x = time, y = distance ), 
              color = "gray", 
              alpha = 0.2) +
  
  geom_point(data = n1_table_CH[which(n1_table_CH$clade232 == 1), ], 
             aes(x = time, y = distance), color = "blue") +
  
  stat_smooth(method = "auto") + 
  
  geom_vline(xintercept = 2012, size = 2, color = "black") + 
  ggtitle("clade232_N1_China") + 
  xlab("")



## combine ----------------

multiplot(clade234_N1_China, clade7_N1_China, clade232_N1_China, ncol = 1)





