source( "./function.coevo.R" )
require( stringr )

raw_na_seq = "./raw_data/processed/pNA_8334.fasta"

raw_na_seq_nt   <- fastaEx( raw_na_seq )$seq
raw_na_seq_id   <- fastaEx( raw_na_seq )$id
raw_na_seq_sero <- str_match( raw_na_seq_id, "\\|_(H5N[0-9]+)_" )[,2]

table( raw_na_seq_sero  )

for( s in 1: length( unique( raw_na_seq_sero ) ) )
{
  type.i = unique( raw_na_seq_sero )[ s ]
  ii     = which( raw_na_seq_sero == type.i )
  
  write.fasta( sequences = raw_na_seq_nt[ii],
               names     = raw_na_seq_id[ii],
               file.out  = paste0( "./aa_numbering/na_type/raw_",  type.i, "_", length(ii),  ".fasta" ) )
}


# manually trim to coding sequence

folder = "./aa_numbering/na_type/"

lsConSeq = list()   
for( d in 1: 9 )
{
  trim_file <- paste0( folder, grep( "trim_", list.files( folder ), value = TRUE )[d] )
  seq0      <- fastaEx( trim_file )$seq
  
  seq_matrix <- do.call( rbind, seq0 )
  
  lsConSeq[[ d ]] =
    apply( seq_matrix, 2, 
           function(x){ 
             con_nt <- which.max( as.data.frame( table( x[ which( x !="~" & x !="-") ] ) , stringsAsFactors = FALSE )[,2] )
             y      <- as.data.frame( table( x[which( x !="~" & x !="-")] ) , stringsAsFactors = FALSE)[,1][ con_nt ]
             
             return(y)
           }) 
  print( d )
}
  
write.fasta( lsConSeq, 
             names    = paste0( "con_N", seq(1,9) ), 
             file.out = "./aa_numbering/na_type/con_Nx.fasta")   

  
# translatorX
# perl translatorx_vLocal.pl -i con_Nx.fasta -p F -t F
# 
# but the web-base performs better than local script 
# probabily becasue 'muscle' is optimized for translatorX
# and muscle is unable to run locally


# link to the working alignments 

n1_align <- "./raw_data/processed/pN1_4696_trim2.fasta"
n2_align <- "./raw_data.n2/processed/pN2_7158.fasta"
n5_align <- "./raw_data.n5/processed/pN5_565.fasta"
n6_align <- "./raw_data.n6/processed/pN6_4286_trim2.fasta"
n8_align <- "./raw_data.n8/processed/pN8_3467.fasta"

nx_align = c( n1_align, n2_align, n5_align, n6_align, n8_align )
for( a in 1:5 )
{
  leafEx( filedir   = nx_align[ a ], 
          leaflist  = fastaEx( nx_align[ a ] )$id, 
          consenseq = TRUE )
}


# pair the working aa sequence to the original aa sequence 

working_fd   = "./aa_numbering/na_type/working_align/"
working_list = list.files( working_fd )
lsConSeq     = "./aa_numbering/na_type/con_Nx.fasta"

for( j in 1:length(working_list) )
{
  working_list.i = paste0( working_fd, working_list[j] )
  sero           = as.numeric( str_match( string = working_list.i, pattern = "pN([0-9])_" )[,2] )
  
  write.fasta( sequences = list(  translate(  fastaEx( working_list.i )$seq[[1]] ), 
                                  translate(  fastaEx( lsConSeq )$seq[[ sero ]] ) ),
               names     = c( fastaEx( working_list.i )$id, fastaEx( lsConSeq )$id[ sero ] ),
               file.out  = paste0( "./aa_numbering/na_type/working_align/paring_", fastaEx( lsConSeq )$id[ sero ], ".fasta" )
               )
  
}










  