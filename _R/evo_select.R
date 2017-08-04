# load-in 

require(adegenet)
require(pegas)
library(seqinr)
library(stringr)
library(ggtree)
library(ggplot2)
library(PopGenome)

setwd("~/Desktop/data_souce/")
source("~/Packaging_Type/_R/Function.R")


### tajima's D --------------------------------

# THIS IS FUNCTION IS MODIFIED FROM PACKAGES: pegas & ape

tajimaGap <- 
  function (x){
  require(ape)
  require(seqinr)
  require(pegas)
  
  
  
  n <- if (is.list(x)) 
    length(x)
  
  else dim(x)[1]
  
  if (n < 4) {
    warning("Tajima test requires at least 4 sequences")
    return(list(D = NaN, Pval.normal = NaN, Pval.beta = NaN))
  }
  
  khat   <- mean(dist.dna(x, "N"))
  ss     <- seg.sites(x)  
  BF.col <- matrix(NA_real_, length(ss), 4)
  for (i in seq_along(ss)) BF.col[i, ] <- base.freq(x[, ss[i]])
  tmpi <- apply(BF.col, 1, function(x) sum(x > 0))
  S    <- sum( tmpi > 1)  
  
  if (!S) {
    warning("no segregating sites")
    return(list(D = NaN, Pval.normal = NaN, Pval.beta = NaN))
  }
  
  tmp <- 1:(n - 1)
  a1  <- sum(1/tmp)
  a2  <- sum(1/tmp^2)
  b1  <- (n + 1)/(3 * (n - 1))
  b2  <- 2 * (n^2 + n + 3)/(9 * n * (n - 1))
  c1  <- b1 - 1/a1
  c2  <- b2 - (n + 2)/(a1 * n) + a2/a1^2
  e1  <- c1/a1
  e2  <- c2/(a1^2 + a2)
  D   <- (khat - S/a1)/sqrt(e1 * S + e2 * S * (S - 1))
  
  Dmin <- (2/n - 1/a1)/sqrt(e2)
  Dmax <- ((n + 1)/(2 * n) - 1/a1)/sqrt(e2)
  tmp1 <- 1 + Dmin * Dmax
  tmp2 <- Dmax - Dmin
  
  a <- -tmp1 * Dmax/tmp2
  b <- tmp1 * Dmin/tmp2
  p <- pbeta((D - Dmin)/tmp2, b, a)
  p <- if (p < 0.5) 
    2 * p
  else 2 * (1 - p)
  list(D = D, Pval.normal = 2 * pnorm(-abs(D)), Pval.beta = p, pi = nuc.div(x))
}

## full HA and NA ----------------

# read-in

time_str_list = list.files("./str_time_CNHK_H5N1/pool_p2yr/")[ 
  grep(".fasta", list.files("./str_time_CNHK_H5N1/pool_p2yr/") ) ]

list_name = str_match(time_str_list, 
                      pattern = "(clade[0-9]{1,3})_([a-z][0-9])-([0-9]{4}-[0-9]{4})" )

clade = c()
gene  = c()
time  = c()
D     = c()
pN    = c()
pB    = c()
pi    = c()


for(k in 1: length(time_str_list) )
{
  bin <- fasta2DNAbin( paste0("./str_time_CNHK_H5N1/pool_p2yr//", time_str_list[k]) )
  out <- tajimaGap(bin)
  
  D[ length(D) + 1 ]         = out$D 
  pN[ length(pN) + 1 ]       = out$Pval.normal
  pB[ length(pB) + 1 ]       = out$Pval.beta
  pi[ length(pi) + 1]        = out$pi
  
  clade[ length(clade) + 1 ] = list_name[,2][k]
  gene[ length(gene) + 1 ]   = list_name[,3][k]
  time[ length(time) + 1 ]   = list_name[,4][k]
  
}

D_result1 <- data.frame(clade, gene, time, D, pN, pB, pi,stringsAsFactors = FALSE)

ggplot(data = D_result1, aes(x = time, y = D, colour = clade, group = clade) ) + 
  facet_wrap(~ gene, ncol = 1) + geom_line(size = 2) + theme_bw()


## HA1  ----------------

time_str_list = list.files("./str_time_CNHK_H5N1/pool_p2yr/HA1/")

list_name = str_match(time_str_list, 
                      pattern = "(clade[0-9]{1,3})_([a-z][0-9])-([0-9]{4}-[0-9]{4})" )

clade = c()
time  = c()
gene  = c()
D     = c()
pN    = c()
pB    = c()
pi    = c()


for(k in 1: length(time_str_list) )
{
  bin <- fasta2DNAbin( paste0("./str_time_CNHK_H5N1/pool_p2yr/HA1/", time_str_list[k]) )
  out <- tajimaGap(bin)
  
  D[ length(D) + 1 ]         = out$D 
  pN[ length(pN) + 1 ]       = out$Pval.normal
  pB[ length(pB) + 1 ]       = out$Pval.beta
  pi[ length(pi) + 1]        = out$pi
  
  clade[ length(clade) + 1 ] = list_name[,2][k]
  gene[ length(gene) + 1 ]   = list_name[,3][k]
  time[ length(time) + 1 ]   = list_name[,4][k]
  
}

D_result_ha1 <- data.frame(clade, time, D, pN, pB, pi, gene, stringsAsFactors = FALSE)


## NAh  ----------------

time_str_list = list.files("./str_time_CNHK_H5N1/pool_p2yr/NAh//")

list_name = str_match(time_str_list, 
                      pattern = "(clade[0-9]{1,3})_([a-z][0-9])-([0-9]{4}-[0-9]{4})" )

clade = c()
gene  = c()
time  = c()
D     = c()
pN    = c()
pB    = c()
pi    = c()


for(k in 1: length(time_str_list) )
{
  bin <- fasta2DNAbin( paste0("./str_time_CNHK_H5N1/pool_p2yr/NAh/", time_str_list[k]) )
  out <- tajimaGap(bin)
  
  D[ length(D) + 1 ]         = out$D 
  pN[ length(pN) + 1 ]       = out$Pval.normal
  pB[ length(pB) + 1 ]       = out$Pval.beta
  pi[ length(pi) + 1]        = out$pi
  
  clade[ length(clade) + 1 ] = list_name[,2][k]
  gene[ length(gene) + 1 ]   = list_name[,3][k]
  time[ length(time) + 1 ]   = list_name[,4][k]
  
}

D_result_nah <- data.frame(clade, time, D, pN, pB, pi, gene, stringsAsFactors = FALSE)

D_result_epi <- rbind(D_result_ha1, D_result_nah)

ggplot(data = D_result_epi, aes(x = time, y = pi, colour = clade, group = clade) ) + 
  facet_wrap(~ gene, ncol = 1) + geom_line(size = 2) + theme_bw()





