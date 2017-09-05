# this file is used to parse hyphy outcome

library(seqinr)
library(stringr)
library(ggtree)
library(ggplot2)
library(readr)
library(dplyr)

setwd("~/Desktop/Hyphy/")
source("~/Packaging_Type/_R/Function.R")

### SLAC (src) --------------------------------

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

# write.csv(slac_0810, "./slac_0810.csv")

slac_0810_p <- slac_0810[ - which( slac_0810$gene == "h5" | slac_0810$gene  == "n1" ),]


## plot:SLAC (src) ----------------

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


### FUBAR (src) --------------------------------

gene  <- c()
clade <- c()
year  <- c()
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

# write.csv(fubar_0810, "./fubar_0810.csv")

fubar_0810_p <- fubar_0810[ - which( fubar_0810$gene == "h5" | fubar_0810$gene  == "n1" ), ]

  
## plot: FUBAR (src) ----------------  

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


### MEME (src) --------------------------------

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


## plot:MEME (src) ----------------

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



### test for sample size (src) --------------------------------

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

### Sampling 2344 2nd time point --------------------------------

setwd("~/Desktop/data_souce/")

# ha1 (n = 112)
rSeq(n = 28, s_times = 10, 
     filedir = "./seq_2yr_CNHK/ha1/p-clade234_h5-2006-2008-ha1.fasta")


# nah (n = 97)
rSeq(n = 24, s_times = 10, 
     filedir = "./seq_2yr_CNHK/nah/p-clade234_n1-2006-2008-nah.fasta")


# ha2 
rSeq(n = 28, s_times = 10, 
     filedir = "./seq_2yr_CNHK/ha2/p-clade234_h5-2006-2008-ha2.fasta")


# nas 
rSeq(n = 24, s_times = 10, 
     filedir = "./seq_2yr_CNHK/nas/p-clade234_n1-2006-2008-nas.fasta")


## SLAC result (s_2344_0608) ----------------

gene       <- c() 
rp         <- c() 
w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/s_2344_06_08/sampling_slac/log/")
log_dir <- paste0("./result/s_2344_06_08/sampling_slac/log/", log_ls) 
txt_ls  <- list.files("./result/s_2344_06_08/sampling_slac/txt/")
txt_dir <- paste0("./result/s_2344_06_08/sampling_slac/txt/", txt_ls)

for(i in 1: length(log_ls) )
{
  str_i   <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})_(r[0-9]{1,2})" )
  
  
  file_log <- read_file(log_dir[i])
  file_txt <- read.table(txt_dir[i], sep = "\t", header = T)
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  
  gene  <- append(gene, str_i[,5] )  
  rp    <- append(rp, str_i[,6] )  
  
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

sampling_2344_0608_slac <- data.frame(gene, rp, w, w_ciu, w_cil, 
                                       dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd)


## FUBAR result (s_2344_0608) ----------------


gene  <- c()
rp    <- c()
po    <- c()
ne    <- c()
dS    <- c()
dN    <- c()
dN_dS <- c()


csv_ls  <- list.files("./result/s_2344_06_08/sampling_fubar/csv/")
log_dir <- paste0("./result/s_2344_06_08/sampling_fubar/csv/", csv_ls) 


for(i in 1: length(csv_ls))
{
  str_i    <- str_match( csv_ls[i], 
                        "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})_(r[0-9]{1,2})" )
  
  file_csv <- read.csv( log_dir[i] )
  
  gene  <- append(gene, str_i[,5] )  
  rp    <- append(rp, str_i[,6] )  
  
  po      <- append( po, sum( file_csv[, 5] >= 0.9 ) )
  ne      <- append( ne, sum( file_csv[, 6] >= 0.9 ) )
  dS      <- append( dS, mean( file_csv[, 2] ) )
  dN      <- append( dN, mean( file_csv[, 3] ) )
  dN_dS   <- append( dN_dS, mean( file_csv[, 4] ) )
  
}

sampling_2344_0608_fubar <- data.frame(gene, rp, 
                                       po, ne, dS, dN, dN_dS, stringsAsFactors = FALSE)


## curate the potential biased results (s_2344_0608) ----------------

sum_sampling_2344_0608_slac = 
sampling_2344_0608_slac %>%
group_by(gene)          %>%
summarise(m.w          = median(w), 
          m.u          = median(w_ciu),
          m.l          = median(w_cil),
          m.dn_mean    = median(dn_mean),
          m.ds_mean    = median(ds_mean),
          m.po         = median(po_site),
          m.ne         = median(ne_site),
          m.dn_ds_mean = median(dn_ds_mean))


sum_sampling_2344_0608_fubar = 
sampling_2344_0608_fubar %>%
group_by(gene)           %>%
summarise( m.po         = median(po),
           m.ne         = median(ne),
           m.ds_mean    = median(dS),
           m.dn_mean    = median(dN),
           m.dn_ds_mean = median(dN_dS))  
  
## fixed version (s_2344_0608) ----------------

F_slac_0810_p = slac_0810_p

F_slac_0810_p[which(F_slac_0810_p$year == "2006-2008" & F_slac_0810_p$clade == "234"), ][, c(4:11)] = 
  sum_sampling_2344_0608_slac[c(2:9)]

F_fubar_0810_p = fubar_0810_p
  
F_fubar_0810_p[which(F_fubar_0810_p$year == "2006-2008" & F_fubar_0810_p$clade == "234"), ][, c(4:8)] = 
  sum_sampling_2344_0608_fubar[c(2, 6)]

## plot: SLAC (s_2344_0608)----------------

# dN - dS
ggplot(data = F_slac_0810_p, aes(x = year, y = dn_ds_mean, colour = clade, group = clade)) + 
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
ggplot(data = F_slac_0810_p, aes(x = year, y = dn_mean, colour = clade, group = clade)) + 
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
ggplot(data = F_slac_0810_p, aes(x = year, y = ds_mean, colour = clade, group = clade)) + 
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

ggplot(data = F_slac_0810_p, aes(x = year, y = w, colour = clade)) +
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

ggplot(data = F_slac_0810_p, aes(x = year, y = ne_site, colour = clade, group = clade)) + 
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



## plot: FUBAR (s_2344_0608) ----------------

# ne_count 

ggplot(data = F_fubar_0810_p, aes(x = year, y = ne, colour = clade, group = clade)) + 
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

ggplot(data = F_fubar_0810_p, aes(x = year, y = po, colour = clade, group = clade)) + 
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





### test for sample size (cut07bird) --------------------------------

gene       <- c()  
clade      <- c()  
year       <- c()  
n          <- c() 

seq_list <- list.files("~/Desktop/data_souce/str_time_CNHK_H5N1/cut2007_bird/")
seq_dir  <- paste0("~/Desktop/data_souce/str_time_CNHK_H5N1/cut2007_bird/", seq_list)

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


### SLAC result (cut07bird) ----------------------------

gene       <- c() 
year       <- c()
clade      <- c()
w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/cut2007_bird_slac/log/")
log_dir <- paste0("./result/cut2007_bird_slac/log/", log_ls) 
txt_ls  <- list.files("./result/cut2007_bird_slac/txt/")
txt_dir <- paste0("./result/cut2007_bird_slac/txt/", txt_ls)

for(i in 1: length(log_ls) )
{
  if( grepl( pattern = "p-", x = log_ls[i]) == FALSE )
  {
    str_i <- str_match(log_ls[i], 
                         "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    
    str_i <- cbind(str_i, str_i[3])
    
  }else
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" )
    }
  
  
  file_log <- read_file( log_dir[i] )
  file_txt <- read.table( txt_dir[i], sep = "\t", header = T )
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  clade <- append(clade, str_i[,2] )
  year  <- append(year, str_i[,4] ) 
  gene  <- append(gene, str_i[,5] )
  
  
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

cut07bird_slac <- data.frame(gene, year, clade, 
                             w, w_ciu, w_cil, 
                             dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd, 
                             stringsAsFactors = FALSE)



## plot: SLAC (cut07bird) ----------------

# dN - dS
ggplot(data = cut07bird_slac, aes(x = year, y = dn_ds_mean, colour = clade, group = clade)) + 
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
ggplot(data = cut07bird_slac, aes(x = year, y = dn_mean, colour = clade, group = clade)) + 
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
ggplot(data = cut07bird_slac, aes(x = year, y = ds_mean, colour = clade, group = clade)) + 
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

ggplot(data = cut07bird_slac, aes(x = year, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1) ) +
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

ggplot(data = cut07bird_slac, aes(x = year, y = ne_site, colour = clade, group = clade)) + 
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





### test for sample size (2yr060810) --------------------------------

gene       <- c()  
clade      <- c()  
year       <- c()  
n          <- c() 

seq_list <- list.files("~/Desktop/data_souce/seq_2yr_060810_CNHK/avian/")
seq_dir  <- paste0("~/Desktop/data_souce/seq_2yr_060810_CNHK/avian/", seq_list)

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




### Sampling 2344 2nd time point (2yr060810) --------------------------------

setwd("~/Desktop/data_souce/")
# 2006 - 2008

# ha1 (n = 97)
rSeq(n = 24, s_times = 20, 
     filedir = "./seq_2yr_060810_CNHK/avian_domain/p-clade234_h5-2006-2008-ha1.fasta")


# nah (n = 83)
rSeq(n = 21, s_times = 20, 
     filedir = "./seq_2yr_060810_CNHK/avian_domain/p-clade234_n1-2006-2008-nah.fasta")


# ha2 
rSeq(n = 24, s_times = 20, 
     filedir = "./seq_2yr_060810_CNHK/avian_domain/p-clade234_h5-2006-2008-ha2.fasta")


# nas 
rSeq(n = 21, s_times = 20, 
     filedir = "./seq_2yr_060810_CNHK/avian_domain/p-clade234_n1-2006-2008-nas.fasta")

### SLAC result (2yr060810) ----------------------------

gene       <- c() 
year       <- c()
clade      <- c()
w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/avian_domain_slac/log/")
log_dir <- paste0("./result/avian_domain_slac/log/", log_ls) 
txt_ls  <- list.files("./result/avian_domain_slac/txt/")
txt_dir <- paste0("./result/avian_domain_slac/txt/", txt_ls)

for(i in 1: length(log_ls) )
{
  if( grepl( pattern = "p-", x = log_ls[i]) == FALSE )
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    
    str_i <- cbind(str_i, str_i[3])
    
  }else
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" )
  }
  
  
  file_log <- read_file( log_dir[i] )
  file_txt <- read.table( txt_dir[i], sep = "\t", header = T )
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  clade <- append(clade, str_i[,2] )
  year  <- append(year, str_i[,4] ) 
  gene  <- append(gene, str_i[,5] )
  
  
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

p2yr060810_slac <- data.frame(gene, year, clade, 
                              w, w_ciu, w_cil, 
                              dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd, 
                              stringsAsFactors = FALSE)



## plot: SLAC (2yr060810) ----------------

# dN - dS
ggplot(data = p2yr060810_slac, aes(x = year, y = dn_ds_mean, colour = clade, group = clade)) + 
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
ggplot(data = p2yr060810_slac, aes(x = year, y = dn_mean, colour = clade, group = clade)) + 
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
ggplot(data = p2yr060810_slac, aes(x = year, y = ds_mean, colour = clade, group = clade)) + 
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

ggplot(data = p2yr060810_slac, aes(x = year, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1) ) +
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

ggplot(data = p2yr060810_slac, aes(x = year, y = ne_site, colour = clade, group = clade)) + 
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


## SLAC result (s_234_06-08_avian_slac) ----------------

gene       <- c() 
rp         <- c() 
w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/sampling_234_06-08_avian_slac/log/")
log_dir <- paste0("./result/sampling_234_06-08_avian_slac/log/", log_ls) 
txt_ls  <- list.files("./result/sampling_234_06-08_avian_slac/txt/")
txt_dir <- paste0("./result/sampling_234_06-08_avian_slac/txt/", txt_ls)

for(i in 1: length(log_ls) )
{
  str_i   <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})_(r[0-9]{1,2})" )
  
  
  file_log <- read_file(log_dir[i])
  file_txt <- read.table(txt_dir[i], sep = "\t", header = T)
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  
  gene  <- append(gene, str_i[,5] )  
  rp    <- append(rp, str_i[,6] )  
  
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

s_234_0608_avian_slac <- data.frame(gene, rp, w, w_ciu, w_cil, 
                                      dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd)



## curate the potential biased results (s_234_06-08_avian_slac) ----------------

sum_s_234_0608_avian_slac = 
  s_234_0608_avian_slac %>%
  group_by(gene)          %>%
  summarise(m.w          = median(w), 
            m.u          = median(w_ciu),
            m.l          = median(w_cil),
            m.dn_mean    = median(dn_mean),
            m.ds_mean    = median(ds_mean),
            m.po         = median(po_site),
            m.ne         = median(ne_site),
            m.dn_ds_mean = median(dn_ds_mean))

F_p2yr060810_slac = p2yr060810_slac

F_p2yr060810_slac[which(F_p2yr060810_slac$year == "2006-2008" & F_p2yr060810_slac$clade == "234"), ][, c(4:11)] = 
  sum_s_234_0608_avian_slac[c(2:9)]


## plot: SLAC (s_234_06-08_avian_slac)----------------

# dN - dS
ggplot(data = F_p2yr060810_slac, aes(x = year, y = dn_ds_mean, colour = clade, group = clade)) + 
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
ggplot(data = F_p2yr060810_slac, aes(x = year, y = dn_mean, colour = clade, group = clade)) + 
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
ggplot(data = F_p2yr060810_slac, aes(x = year, y = ds_mean, colour = clade, group = clade)) + 
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

ggplot(data = F_p2yr060810_slac, aes(x = year, y = w, colour = clade)) +
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

ggplot(data = F_p2yr060810_slac, aes(x = year, y = ne_site, colour = clade, group = clade)) + 
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












## SLAC result (str_0826 - avian_SW) ----------------

gene       <- c() 
clade      <- c()
region     <- c()

w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/str_0826/avian_SW/log/")
log_dir <- paste0("./result/str_0826/avian_SW/log/", log_ls) 
txt_ls  <- list.files("./result/str_0826/avian_SW/txt/")
txt_dir <- paste0("./result/str_0826/avian_SW/txt/", txt_ls)


for(i in 1: length(log_ls) )
{
  str_i    <- str_match(log_ls[i], "([0-9]{1,3})_([a-z0-9]+)_([a-zA-Z]{2,3})" )
  
  file_log <- read_file( log_dir[i] )
  file_txt <- read.table( txt_dir[i], sep = "\t", header = T )
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  clade  <- append(clade, str_i[,2] )
  gene   <- append(gene, str_i[,3] )
  region <- append(region, str_i[,4] )
  
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

str_0826_avian_SW_slac <- 
  data.frame( gene, clade, region, 
              w, w_ciu, w_cil, 
              dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd, stringsAsFactors = FALSE)



## plot: SLAC (str_0826 - avian_SW) ----------------


# w 

ggplot(data = str_0826_avian_SW_slac, aes(x = region, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1) ) +
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

# 234 only
str_0826_avian_SW_slac_234 <- str_0826_avian_SW_slac[which(str_0826_avian_SW_slac$gene == "h5" & str_0826_avian_SW_slac$clade == "234"), ]

str_0826_avian_SW_slac_234$region[1] = "Others"

ggplot(data = str_0826_avian_SW_slac_234, aes(x = region, y = w, color = region) ) +
  geom_point(size = 2) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), width = 0.1, size = 1.2) +
  scale_color_manual( values = c("black", "chartreuse2") ) + 
  xlab("") + ylab(expression(omega)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.ticks = element_blank(), 
        panel.border = element_blank(),
        axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 20), 
        legend.position = "none")


## SLAC result (str_0826 - avian_2yr) ----------------

gene       <- c() 
year       <- c()
clade      <- c()

w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/str_0826/avian_2yr/log/")
log_dir <- paste0("./result/str_0826/avian_2yr/log/", log_ls) 
txt_ls  <- list.files("./result/str_0826/avian_2yr/txt/")
txt_dir <- paste0("./result/str_0826/avian_2yr/txt/", txt_ls)

for(i in 1: length(log_ls) )
{
  if( grepl( pattern = "p-", x = log_ls[i]) == FALSE )
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    
    str_i <- cbind(str_i, str_i[3])
    
  }else
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" )
  }
  
  
  file_log <- read_file( log_dir[i] )
  file_txt <- read.table( txt_dir[i], sep = "\t", header = T )
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  clade <- append(clade, str_i[,2] )
  year  <- append(year, str_i[,4] ) 
  gene  <- append(gene, str_i[,5] )
  
  
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

str_0826_avian_2yr_slac <- 
  data.frame( gene, year, clade, 
              w, w_ciu, w_cil, 
              dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd, stringsAsFactors = FALSE)


## plot: SLAC (str_0826 - avian_2yr) ----------------

# w 

ggplot(data = str_0826_avian_2yr_slac, aes(x = year, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1) ) +
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


str_0826_avian_2yr_slac_23_full <- str_0826_avian_2yr_slac[which(str_0826_avian_2yr_slac$clade == "232" | str_0826_avian_2yr_slac$clade == "234"), ]
str_0826_avian_2yr_slac_23_full <- str_0826_avian_2yr_slac_23_full[which(str_0826_avian_2yr_slac_23_full$gene == "h5"| str_0826_avian_2yr_slac_23_full$gene == "n1"),]

# 232/234 only - h5
ggplot(data = str_0826_avian_2yr_slac_23_full, 
       aes(x = year, y = w, colour = clade)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1, size = 1.2) +
  facet_wrap(~gene, ncol = 2) + 
  scale_y_continuous(limits = c(0.1,0.4) ) +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) +
  scale_x_discrete(labels = c("2004-2005", "2006-2007", "2008-2009")) +
  xlab("") + ylab(expression(omega)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 22), 
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 20, face = "bold"))
  
## SLAC result (str_0826 - avian_SW_1yr_234) ----------------

year       <- c()
region     <- c()
gene       <- c()

w          <- c()  
w_ciu      <- c()  
w_cil      <- c()
dn_mean    <- c()  
ds_mean    <- c() 
dn_ds_mean <- c() 
dn_ds_sd   <- c() 

po_site    <- c()  
ne_site    <- c()  

log_ls  <- list.files("./result/str_0826/avian_SW_1yr_234/log/")
log_dir <- paste0("./result/str_0826/avian_SW_1yr_234/log/", log_ls) 
txt_ls  <- list.files("./result/str_0826/avian_SW_1yr_234/txt/")
txt_dir <- paste0("./result/str_0826/avian_SW_1yr_234/txt/", txt_ls)

for(i in 1: length(log_ls) )
{
  if( grepl( pattern = "p-", x = log_ls[i]) == FALSE )
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)_([a-zA-Z]{2,3})-([0-9]{4}-[0-9]{4})" )
    
    str_i <- cbind(str_i, str_i[3])
    
  }else
  {
    str_i <- str_match(log_ls[i], 
                       "([0-9]{1,3})_([a-z0-9]+)_([a-zA-Z]{2,3})-([0-9]{4}-[0-9]{4})-([a-z0-9]{3})" )
  }
  
  file_log <- read_file( log_dir[i] )
  file_txt <- read.table( txt_dir[i], sep = "\t", header = T )
  
  mean  <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,2]
  U     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,3]
  L     <- str_match(file_log, "Using dN/dS=([0-9.]+)\\(Estimated 95% CI = \\[([0-9.]+),([0-9.]+)" )[,4]
  
  mean  <- as.numeric(mean)
  U     <- as.numeric(U)
  L     <- as.numeric(L)
  
  N        <- mean( file_txt[,8] ) 
  S        <- mean( file_txt[,7] )
  dnds     <- mean( file_txt[, 9] ) 
  dnds_sd  <- sd( file_txt[,9] )
  po       <- sum( file_txt[,10][ which(file_txt[,10] != 0) ] < 0.1)  
  ne       <- sum( file_txt[,11][ which(file_txt[,11] != 0) ] < 0.1)  
  
  year   <- append(year, str_i[,5] ) 
  gene   <- append(gene, str_i[,6] )
  region <- append(region, str_i[,4] )
  
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

str_0826_avian_SW_1yr_234_slac <- 
  data.frame( year, gene, region,
              w, w_ciu, w_cil, 
              dn_mean, ds_mean, po_site, ne_site, dn_ds_mean, dn_ds_sd, stringsAsFactors = FALSE)

## plot: SLAC (str_0826 - avian_2yr) ----------------

# w 

ggplot(data = str_0826_avian_SW_1yr_234_slac, aes(x = year, y = w, colour = region)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1) +
  
  scale_y_continuous(limits = c(0,1) ) +
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



str_0826_avian_SW_1yr_234_slac_h5 = str_0826_avian_SW_1yr_234_slac[which(str_0826_avian_SW_1yr_234_slac$gene == "h5"), ]
str_0826_avian_SW_1yr_234_slac_h5$region[c(1,2,3,4)] = "Others"

# ha only 

ggplot(data = str_0826_avian_SW_1yr_234_slac_h5, aes(x = year, y = w, colour = region, group = region)) +
  geom_point( position = position_dodge(0.8) ) +
  geom_errorbar( aes(ymin = w_cil, ymax = w_ciu), 
                 position = position_dodge(0.8), width = 0.1, size = 1.2) +
  
  scale_x_discrete(labels = c("2004-2005", "2006", "2007", "2008-2009")) +
  # geom_line(position = position_dodge(0.8)) +
  scale_y_continuous(limits = c(0,0.5) ) +
  scale_color_manual(values = c("black", "chartreuse2") ) + 
  xlab("") + ylab(expression(omega)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 20), 
        panel.border = element_rect(color = "black", fill = NA, size = 1) )





