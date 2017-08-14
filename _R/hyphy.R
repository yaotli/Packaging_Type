
library(seqinr)
library(stringr)
library(ggtree)
library(ggplot2)
library(readr)

setwd("~/Desktop/Hyphy/")


### SLAC  --------------------------------

gene       <- c()  
clade      <- c()  
year       <- c()  
w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

for(g in 1: 5)
{
  g_cat   <- c("full", "ha1", "ha2", "nah", "nas")
  
  log_ls  <- grep(".log", list.files( paste0("./result/src_slac/", g_cat[g], "/done") ), 
                  value = T)
  log_dir <- paste0("./result/src_slac/", g_cat[g], "/done/", log_ls)

  txt_ls  <- list.files( paste0("./result/src_slac/", g_cat[g], "/txt") )
  txt_dir <- paste0("./result/src_slac/", g_cat[g], "/txt/", txt_ls)
             
  
  if(g == 1)
  {
    str_log <- str_match(string = log_ls, 
                          "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    
    str_log = cbind(str_log, str_log[,3] )
    
    
    str_txt <- str_match(string = txt_ls, 
                         "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    
    str_txt = cbind(str_txt, str_txt[,3] )
    
  }else
  {
    str_log <- str_match(string = log_ls, 
                          "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" )
    
    str_txt <- str_match(string = txt_ls, 
                          "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" )
    }

# parse log 
  
  for(l in 1: length(log_ls) )
  {
    file_log <- read_file(log_dir[l])
    mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
    U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
    L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
    
    mean <- as.numeric(mean)
    U    <- as.numeric(U)
    L    <- as.numeric(L)
    
    
    file_txt <- read.table(txt_dir[l], sep = "\t", header = T)
    N        <- mean( file_txt[,8] ) 
    S        <- mean( file_txt[,7] )
    dnds     <- mean( file_txt[, 9] ) 
    dnds_sd  <- sd( file_txt[,9] )
    po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
    ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
      
    if ( TRUE %in% (str_log[,1] != str_txt[,1]) )
    {
      print("ERROR")
      print(g)
      print(l)
      
    }else
    {
      gene  <- append(gene, str_log[,5][l] )  
      clade <- append(clade, str_log[,2][l] )    
      year  <- append(year, str_log[,4][l] )
      
      w     <- append(w, mean)    
      w_ciu <- append(w_ciu, U)  
      w_cil <- append(w_cil, L) 
      
      dn_mean    <- append(dn_mean, N)    
      ds_mean    <- append(ds_mean, S)      
      po_site    <- append(po_site, po)      
      ne_site    <- append(ne_site, ne) 
      dn_ds_mean <- append(dn_ds_mean, dnds) 
      dn_ds_sd   <- append(dn_ds_sd, dnds_sd) 
      }
  }  

}


slac_0810   <- data.frame( gene, clade, year, w, w_ciu, w_cil, 
                           dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd,
                           stringsAsFactors = FALSE )

slac_0810_p <- slac_0810[ - which( slac_0810$gene == "h5" | slac_0810$gene  == "n1" ),]

## plot ----------------

# dN - dS
ggplot(data = slac_0810_p, aes(x = year, y = dn_ds_mean, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("dN - dS") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )
  
# dN
ggplot(data = slac_0810_p, aes(x = year, y = dn_mean, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("dN") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )

# dS
ggplot(data = slac_0810_p, aes(x = year, y = ds_mean, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("dS") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )


# w 

ggplot(data = slac_0810, aes(x = year, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1.5) ) +
  facet_wrap(~gene, ncol = 2) + 
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  xlab("") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 12, face = "bold")) 
  
  
  
# ne_count 

ggplot(data = slac_0810_p, aes(x = year, y = ne_site, colour = clade, group = clade)) + 
geom_point( size = 2 ) +
geom_line( size = 1.5 ) +
facet_wrap(~gene, ncol = 2) + 
theme_bw() + 
theme(panel.grid.minor = element_blank(), 
      panel.grid.major.y = element_blank(), 
      panel.grid.major.x = element_blank() ) + 
xlab("") + ylab("site under ne. selection") +
scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
theme(
  panel.border = element_rect(color = "black", fill = NA, size = 1), 
  strip.background = element_rect(fill = "white", color = "white"),
  strip.text = element_text(size = 12, face = "bold") )


### FUBAR  --------------------------------

gene  <- c()
clade <- c()
year  <- c()
codon <- c()
po    <- c()
ne    <- c()
dS    <- c()
dN    <- c()
dN_dS <- c()
  
  
for(g in 1: 5)
{
  g_cat   <- c("full", "ha1", "ha2", "nah", "nas")
  
  csv_ls  <- list.files( paste0("./result/src_fubar/", g_cat[g], "/csv") )
  csv_dir <- paste0("./result/src_fubar/", g_cat[g], "/csv/", csv_ls)
  
  if(g == 1)
  {
    str_csv <- str_match(string = csv_ls, 
                         "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    
    str_csv = cbind(str_csv, str_csv[,3])
    
  }else
  {
    str_csv <- str_match(string = csv_ls, 
                         "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" )
  }
  
  # parse csv 
  
  for(l in 1: length(csv_ls) )
  {
    readin <- read.csv( csv_dir[l] )

    gene    <- append( gene, str_csv[,5][l] )  
    clade   <- append( clade, str_csv[,2][l] )    
    year    <- append( year, str_csv[,4][l] )
    
    po      <- append( po, sum( readin[, 5] >= 0.9 ) )
    ne      <- append( ne, sum( readin[, 6] >= 0.9 ) )
    dS      <- append( dS, mean( readin[, 2] ) )
    dN      <- append( dN, mean( readin[, 3] ) )
    dN_dS   <- append( dN_dS, mean( readin[, 4] ) )

  }
      
}    
  
  
  fubar_0810   <- data.frame(gene, clade, year, po, ne, 
                             dS, dN, dN_dS, stringsAsFactors = FALSE)
  
  fubar_0810_p <- fubar_0810[ - which( fubar_0810$gene == "h5" | fubar_0810$gene  == "n1" ), ]

  
## plot ----------------  

  # ne_count 

ggplot(data = fubar_0810_p, aes(x = year, y = ne, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("site under ne. selection") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )
  
  # po_count
  
ggplot(data = fubar_0810_p, aes(x = year, y = po, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("site under po. selection") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )


### MEME --------------------------------

gene       <- c()  
clade      <- c()  
year       <- c()  
w          <- c()  
w_ciu      <- c()  
w_cil      <- c()



log_ls  <- grep(".log", list.files( paste0("./result/meme_log/") ), value = T)
log_dir <- paste0("./result/meme_log/", log_ls)

for (i in 1: length(log_ls) )
{
  if( grepl("p-", log_ls[i]) == FALSE)
  {
    str_log <- str_match(log_ls[i], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    str_log <- c(str_log, str_log[, 3])
    
  }else
  {
    str_log <- str_match(log_ls[i], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" ) 
    }
  
  file_log <- read_file( log_dir[i] )
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean <- as.numeric(mean)
  U    <- as.numeric(U)
  L    <- as.numeric(L)
  
  gene       <- c(gene, str_log[5] )  
  clade      <- c(clade, str_log[2] )  
  year       <- c(year, str_log[4] )  
  w          <- c(w, mean)  
  w_ciu      <- c(w_ciu, U)  
  w_cil      <- c(w_cil, L)
  
}


meme_w <- data.frame(gene, clade, year, w, w_ciu, w_cil, stringsAsFactors = FALSE)


## plot ----------------

ggplot(data = meme_w, aes(x = year, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1.5) ) +
  facet_wrap(~gene, ncol = 2) + 
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  xlab("") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 12, face = "bold")) 



### test for sample size --------------------------------

gene       <- c()  
clade      <- c()  
year       <- c()  
n          <- c() 

seq_list <- list.files("~/Desktop/data_souce/seq_2yr_CNHK/full/")
seq_dir  <- paste0("~/Desktop/data_souce/seq_2yr_CNHK/full/", seq_list)

for(j in 1: length(seq_list))
{
  str_seq <- str_match( seq_list[j], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
  n_j     <- length( read.fasta(seq_dir[j]) )
  
  gene  = c(gene, str_seq[3])
  clade = c(clade, str_seq[2])
  year  = c(year, str_seq[4])
  
  n     = c(n, n_j)
}

testsamplesize <- data.frame(gene, clade, year, n, stringsAsFactors = FALSE)


ggplot(data = testsamplesize, aes(x = year, y = n, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("no. of sequence") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )












