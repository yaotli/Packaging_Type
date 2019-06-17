source( "./function.coevo.R" )

require( ape )

# extract gsgd-like and rebuild a tree 

ls.gsgd <- tagExtra( "./raw_data/processed/tree/fasttree_pH5_8334_e0921.tre" )
ls.gsgd <- ls.gsgd$id[ !is.na( ls.gsgd$tag ) ]
leafEx( "./raw_data/processed/pH5_8334_trim2.2.fasta", ls.gsgd )

# gsgd-strict sequences 

ls.gsgd_s <- read.beast( "./gsgd/processed/tree/fasttree_pH5_gsgd_7252_e0921.tre" )
ls.gsgd_s <- fortify( ls.gsgd_s )
colnames( ls.gsgd_s )[ ncol( ls.gsgd_s  ) ] = "br_col"
 
id.gsgd_s <- ls.gsgd_s$label[ getDes( tre.d = ls.gsgd_s, node = which( !is.na( ls.gsgd_s$br_col ) ) ) ]
id.gsgd_s <- id.gsgd_s[ which( !is.na(id.gsgd_s) ) ]
leafEx( "./raw_data/processed/pH5_8334_trim2.2.fasta",  id.gsgd_s )

rmDup( "./gsgd/processed/pH5_8334_trim2.2_7245.fasta", rmdup = TRUE )
rmdup_plus( "./gsgd/processed/pH5_8334_trim2.2_7245_rmd.fasta" )



# extract c2344 and remove duplicates 

ls.c2344 <- tagExtra( "./gsgd/processed/tree/fasttree_pH5_gsgd_7252_e0921.tre" )
ls.c2344 <- ls.c2344$id[ !is.na( ls.c2344$tag ) ]
leafEx( "./raw_data/processed/pH5_8334_trim2.2.fasta", ls.c2344 )

rmDup( "./gsgd/processed/pH5_8334_trim2.2_2659.fasta", rmdup = TRUE )
rmdup_plus( "./gsgd/processed/pH5_8334_trim2.2_2659_rmd.fasta" )

# edit IQ-TREE

raw_iqtree           <- read.tree( "./gsgd/processed/tree/iqtree_c2344_rm2/pH5_8334_trim2.2_2659_rmd_rmd2.tre" )
raw_iqtree$tip.label <- sub( sub( raw_iqtree$tip.label, pattern = "__", replacement = "_\\|" ), pattern = "__", replacement = "\\|_" ) 
write.tree( raw_iqtree, file = "./gsgd/processed/tree/iqtree_c2344_rm2/iqtree_pH5_c2344_1600.tre" )

# find 3 outliers in TempEst 
# JX878683_duck_HuBei_03_2010_|China|_H5N5_2010.236
# EPI237995_turkey_Minnesota_9845_4_2015_|USA|_H5N2_2005.000
# LC316683_tundra_swan_Niigata_1511C003_2016_|Japan|_H5N6_2016.904
# remove the 1 and 3 of the above an corrected the 2 time to 2015.000

n1600.seqname = fastaEx( "./gsgd/processed/pH5_8334_trim2.2_2659_rmd_rmd2.fasta" )$id
n1600.seq     = fastaEx( "./gsgd/processed/pH5_8334_trim2.2_2659_rmd_rmd2.fasta" )$seq

n1600.seqname[ grep( "EPI237995", n1600.seqname ) ] = 
  sub( "_2005.000", "_2015.000", n1600.seqname[ grep( "EPI237995", n1600.seqname ) ] )

rm <- grep( "JX878683|LC316683", n1600.seqname )
write.fasta( n1600.seq[-rm], n1600.seqname[-rm], file.out = "./gsgd/processed/pH5_c2344_1598.fasta" )


# edit IQ-TREE ( correct '__' back )

raw_iqtree           <- read.nexus( "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1014.tre" )
raw_iqtree$tip.label <- sub( sub( raw_iqtree$tip.label, pattern = "__", replacement = "_\\|" ), pattern = "__", replacement = "\\|_" ) 
write.nexus( raw_iqtree, file = "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1014_.tre" )


# searching for the reassortment events 

# tipls_n1958 <- taxaInfo( useTree = TRUE, file = "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1014_.tre", makecsv = FALSE )
# 
# cladeSampling( trefile   = "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1014_.tre",
#                fasfile   = "./raw_data/processed/pH5_8334_trim2.2.fasta",
#                suppList  = TRUE,
#                listinput = tipls_n1958,
#                grid      = 100,
#                list.x    = c( 6, 4, 3), saveFasta = FALSE )



# extract early strains for ancestral reconstruction 

tipls_n1958 <- taxaInfo( useTree = TRUE, file = "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1110_.tre", makecsv = FALSE )

anc_recon   <- c( tipls_n1958[[6]][ which( tipls_n1958[[4]] < 2014 ) ],
                  tipls_n1958[[6]][ which( tipls_n1958[[7]] == "ff6fcf" &  tipls_n1958[[4]] < 2014.1 &  
                                         ( tipls_n1958[[2]] == "South_Korea" | tipls_n1958[[2]] == "China" ) ) ]  ) #n139

leafEx( filedir = "./raw_data/processed/pH5_8334_trim2.2.fasta", consenseq = FALSE, leaflist = anc_recon )
timeDice( "./gsgd/processed/pH5_c2344anc_139.fasta", ecotable = FALSE,
          oldfas.dir = c( "./raw_data/sources/pH5_G_2690_20180912.fasta", 
                          "./raw_data/sources/pH5_N_6452_20180912.fasta") )

# exclude old H5N1

old_h5n1 <- tipls_n1958[[6]][ which( tipls_n1958[[7]] == 808080 ) ] #n124
anc_2344 <- setdiff( anc_recon, old_h5n1 )
leafEx( filedir = "./raw_data/processed/pH5_8334_trim2.2.fasta", consenseq = FALSE, leaflist = anc_2344 )

b2014      <- tipls_n1958[[6]][ which( tipls_n1958[[4]] < 2014 ) ]
b2014_2344 <- setdiff( b2014, old_h5n1 )
leafEx( filedir = "./raw_data/processed/pH5_8334_trim2.2.fasta", consenseq = FALSE, leaflist = b2014_2344 )



# examine each clade 

table(tipls_n1958[[7]])
cladeSampling( trefile   = "./gsgd/processed/tree/iqtree_c2344_1598/pH5_c2344_1598_e1014_.tre",
               fasfile   = "./raw_data/processed/pH5_8334_trim2.2.fasta",
               suppList  = TRUE,
               listinput = tipls_n1958,
               grid      = 0.5,
               list.x    = c( 6, 4, 7), saveFasta = TRUE )

for( l in 1:4 )
{
  col = unique( tipls_n1958[[7]] )[ c(1,2,5,6) ]
  group_i <- intersect( tipls_n1958[[6]][ which( tipls_n1958[[7]] == col[l]) ], fastaEx("./gsgd/processed/pH5_c2344_1598_cs2.fasta")$id )
  
  leafEx( filedir = "./raw_data/processed/pH5_8334_trim2.2.fasta", consenseq = FALSE, leaflist = group_i )
}



# prepare h5 big tree ( with redundant srquences )
# curate as previous nonredundant alignment  

n2659.seqname = fastaEx( "./gsgd/processed/pH5_8334_trim2.2_2659.fasta" )$id
n2659.seq     = fastaEx( "./gsgd/processed/pH5_8334_trim2.2_2659.fasta" )$seq

n2659.seqname[ grep( "EPI237995", n2659.seqname ) ] = 
  sub( "_2005.000", "_2015.000", n2659.seqname[ grep( "EPI237995", n2659.seqname ) ] )

rm <- grep( "JX878683|LC316683", n2659.seqname )
write.fasta( n2659.seq[-rm], n2659.seqname[-rm], file.out = "./gsgd/processed/pH5_c2344_2657.fasta" )


# prepare big tree for ancestral aa reconstruction

n2657.seqname = fastaEx( "./gsgd/processed/pH5_c2344_2657.fasta" )$id
n2657.seq     = fastaEx( "./gsgd/processed/pH5_c2344_2657.fasta" )$seq
mix <- grep( "mix", n2657.seqname, ignore.case = TRUE )
write.fasta( n2657.seq[-mix], n2657.seqname[-mix], file.out = "./gsgd/processed/pH5_c2344_2629.fasta" )

rmDup( "./gsgd/processed/pH5_c2344_2629.fasta", rmdup = TRUE )
rmdup_plus( "./gsgd/processed/pH5_c2344_2629_rmd.fasta" )


# subsampling 
# criteria: one sample with the same country and the same serotype in time period - 

n_1547 = "./gsgd/processed/pH5_c2344_2629_rmd_rmd2.fasta"

in.seq = fastaEx( n_1547 )$seq
in.id  = fastaEx( n_1547 )$id

set.seed( 1000 )
con.i  = str_match( in.id, "\\|([A-Za-z_]+)\\|" )[,2]
sero.i = str_match( in.id, "\\|_(H5N[0-9])_" )[,2]
time.i = str_match( in.id, "[0-9.]+$" )

# set grid
time.t = round( as.numeric(time.i), digit = 2 )

samp= c()
for( y in 1: length( unique(as.character(time.t)) ) )
{
  ii = which( time.t == unique( as.character(time.t) )[y] )  # time 
  if( length(ii) == 1 ){ samp = c( samp, ii ) }else
  {
    con.ii   = con.i[ ii ]
    sero.ii  = sero.i[ ii ]
    con_sero = paste0( con.ii, sero.ii )
    
    for( z in 1: length( unique( con_sero ) ) )
    {
      jj = which( con_sero == unique(con_sero)[z] ) # sero + country
      if( length(jj) == 1 ){ samp = c( samp, ii[jj] ) }else
      {
        samp = c( samp, sample( ii[jj], size = 1) )
      } 
    }
  }
}

write.fasta( in.seq[samp], in.id[samp], file.out = "./gsgd/processed/pH5_c2344_1547ts.fasta" )






















