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


## tajima's D function --------------------------------

# THIS IS FUNCTION IS MODIFIED FROM PACKAGES: pegas & ape

tajimaGap <- 
  function (x)
{
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

### 2yr_CNHK --------------------------------

# read-in

fas_ls  <- list.files( "./pool_seq_2yr_CNHK/" )
log_dir <- paste0( "./pool_seq_2yr_CNHK/", fas_ls )

clade = c()
gene  = c()
year  = c()
D     = c()
pN    = c()
pB    = c()
pi    = c()


for(k in 1: length(log_dir) )
{
  
  if( grepl("p-", fas_ls[k]) == FALSE)
  {
    str_fas <- str_match(fas_ls[k], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    str_fas <- c(str_fas, str_fas[, 3])
    
  }else
  {
    str_fas <- str_match(fas_ls[k], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" ) 
  }
  
  bin <- fasta2DNAbin( log_dir[k] )
  out <- tajimaGap(bin)
  
  D[ length(D) + 1 ]         = out$D 
  pN[ length(pN) + 1 ]       = out$Pval.normal
  pB[ length(pB) + 1 ]       = out$Pval.beta
  pi[ length(pi) + 1]        = out$pi
  
  gene  <- c( gene, str_fas[5] )  
  clade <- c( clade, str_fas[2] )  
  year  <- c( year, str_fas[4] )  
  
}

D_full <- data.frame(clade, gene, year, D, pN, pB, pi,stringsAsFactors = FALSE)

D_full_f <- D_full[ which( D_full$gene == "h5" | D_full$gene  == "n1" ),]


## plot ----------------

ggplot(data = D_full_f, aes(x = year, y = D, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("Tajima's D") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )

### per year --------------------------------

# read-in

fas_ls  <- list.files( "./pool_seq_pyr_CNHK/" )
log_dir <- paste0( "./pool_seq_pyr_CNHK/", fas_ls )

clade = c()
gene  = c()
year  = c()
D     = c()
pN    = c()
pB    = c()
pi    = c()


for(k in 1: length(log_dir) )
{
  
  if( grepl("p-", fas_ls[k]) == FALSE)
  {
    str_fas <- str_match(fas_ls[k], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})" )
    str_fas <- c(str_fas, str_fas[, 3])
    
  }else
  {
    str_fas <- str_match(fas_ls[k], "([0-9]{1,3})_([a-z0-9]+)-([0-9]{4}-[0-9]{4})-([a-z0-9]{3,4})" ) 
  }
  
  testN <- length( attributes( read.fasta( log_dir[k] ) )$names )
  
  if( testN > 3)
  {
    bin <- fasta2DNAbin( log_dir[k] )
    out <- tajimaGap(bin)
    
    D[ length(D) + 1 ]         = out$D 
    pN[ length(pN) + 1 ]       = out$Pval.normal
    pB[ length(pB) + 1 ]       = out$Pval.beta
    pi[ length(pi) + 1]        = out$pi
    
    gene  <- c( gene, str_fas[5] )  
    clade <- c( clade, str_fas[2] )  
    year  <- c( year, str_fas[4] )  
    
  }
  
}

D_full_p1y   <- data.frame(clade, gene, year, D, pN, pB, pi,stringsAsFactors = FALSE)

D_full_f_p1y <- D_full_p1y[ which( D_full_p1y$gene == "h5" | D_full_p1y$gene  == "n1" ),]

D_full_f_p1y <- D_full_f_p1y[ which( D_full_f_p1y$year == "2002-2003" | D_full_f_p1y$year == "2003-2004" |
                                     D_full_f_p1y$year == "2004-2005" | D_full_f_p1y$year == "2005-2006" | 
                                     D_full_f_p1y$year == "2006-2007" | D_full_f_p1y$year == "2007-2008" |
                                     D_full_f_p1y$year == "2008-2009" | D_full_f_p1y$year == "2009-2010"), ]

## plot ----------------

ggplot(data = D_full_f_p1y, aes(x = year, y = D, colour = clade, group = clade)) + 
  geom_point( size = 2 ) +
  geom_line( size = 1.5 ) +
  facet_wrap(~gene, ncol = 2) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank() ) + 
  xlab("") + ylab("Tajima's D") +
  scale_color_manual(values = c("#228B22", "#EE2C2C", "#000000") ) + 
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 12, face = "bold")
  )

