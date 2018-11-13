source( "./function.coevo.R" )
require( seqinr )
require( tidyverse )
require( stringr )
require( bio3d )

# sequence from prepared dataset 
rawseqH5     <- fastaEx( "./raw_data/processed/pH5_8334_trim2.2.fasta" )$seq
rawidH5      <- fastaEx( "./raw_data/processed/pH5_8334_trim2.2.fasta" )$id
vn1203.data  <- translate( rawseqH5[[ grep( "viet_nam_1203_2004", rawidH5, ignore.case = TRUE )[1] ]] )

# sequence from the table
numtable     <- read.csv( "./aa_numbering//numbering.csv", stringsAsFactors = FALSE )
vn1203.table <-  as.vector( na.omit( str_match( numtable$A.Vietnam.1203.2004.H5N1, "([0-9]{1,3}) ([A-Z]{1})" )[,3], 2 ) )

# 2fk0 & 4jul structure 
pdb.2fk0   <- read.pdb( "./aa_numbering/2fk0.pdb" )
pdb.4jul   <- read.pdb( "./aa_numbering/4jul.pdb" )

# extract real residues in A,B chains 
st.files = c( "pdb.2fk0", "pdb.4jul" )

st.seq = list()
st.no  = list()
for( i in 1:2 )
{
  st = get( st.files[i] )
  
  tem.table = 
  st$atom %>%
    filter( type == "ATOM" ) %>%
    filter( chain == "A" | chain == "B" ) %>%
    select( resid, resno, insert, chain ) 
  
  for( j in 1:2 )
  {
    st.1 = 
      tem.table %>% 
      filter( chain == c("A", "B")[j] )
    
    st.no[[ length(st.no) + 1 ]] = gsub( "NA$", "", unique( paste0( st.1$resno, st.1$insert ) ) )
    
    aa.i = st.1$resid[ which(  !duplicated( paste0( st.1$resno, st.1$insert ) ) ) ]
    st.seq[[ length(st.seq) +1 ]] = 
      sapply( as.list( aa.i ), 
              function(x)
              {
                y    = tolower( s2c(x) )
                y[1] = toupper( y[1] )  
                y    = paste0( y, collapse = "" )
                return( a(y) )
              })
  }
}

# write a ref. fasta and csv
seqinr::write.fasta( sequences = list( vn1203.table, st.seq[[1]], st.seq[[2]], st.seq[[3]], st.seq[[4]], vn1203.data ), 
             names     = c( "vn1203.table", "vn1203.c1", "vn1203.c2", "anhui01.c1", "anhui01.c2", "vn1203.data" ), 
             file.out  = "./aa_numbering/aa_ref.fasta" )

anno_numtable <- 
data.frame( numtable, 
            p2fk0_1 = c( st.no[[1]], rep("-", dim(numtable)[1] - length(st.no[[1]])) ), 
            p2fk0_2 = c( st.no[[2]], rep("-", dim(numtable)[1] - length(st.no[[2]])) ),
            p4jul_1 = c( st.no[[3]], rep("-", dim(numtable)[1] - length(st.no[[3]])) ),
            p4jul_2 = c( st.no[[4]], rep("-", dim(numtable)[1] - length(st.no[[4]])) ) )
 
write.csv( anno_numtable, file = "./aa_numbering/ref_numbering.csv", row.names = FALSE)                  

 

#ha_num( ref_csv  = "./aa_numbering/ref_numbering.csv", ref_fas = "./aa_numbering/align_aa_ref.fasta", 
#        data_pos = c( 224, 226, 170) )



