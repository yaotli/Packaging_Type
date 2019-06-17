# check if 2 codons represent the same amino acid
.check_rep = function( codon1, codon2 )
{
  a1 = translate( codon1 )
  a2 = translate( codon2 )
  
  if( a1 == "X" | a2 == "X" ){ out =  FALSE  }else
  {
    if( a1 != a2 ){ out = TRUE }else
    {
      out = FALSE
    }
  }
  return( out )
}



# proportional site counting method (Bhatt et al., 2010)

prop_site_count_F = c( 0, 1, 0, 0.5, 0, 1/3, 0 )
prop_site_count_P = c( 0, 0, 1, 0.5, 1, 2/3, 1 )
prop_site_count_x = c( "01", "10", "15", "20", "25", "30", "35" )

.rep.p = function( mx, out_conseq = out_conseq )
{
  R_P = 
    apply( mx, 2, 
           function(x)
           {
             
             char.i = grep( pattern = "a|t|c|g", x, ignore.case = TRUE, value = TRUE )
             char   = length( char.i )
             m.i    = char.i[ which( char.i != out_conseq[ as.numeric(x[1]) ] ) ]
             
             if( char < 1)
             {
               P.i = 0
               R.i = 0
             }else
             {
               # P.i - polymorphism score
               
               out.p  = ( char - length( m.i ) )/char
               out.p  = ifelse( ( out.p<1 ) & (out.p>0), 5, out.p  )
               
               type.i = length( unique( m.i ) )
               x.i    = paste0( type.i, out.p )
               
               P.i = prop_site_count_P[ match( x.i, prop_site_count_x ) ]   
               
               # R.i - replacement score 
               
               codon.n = ceiling( (as.numeric(x[1]))/3 ) - 1
               codon   = out_conseq[ c( codon.n*3+1, codon.n*3+2, codon.n*3+3) ]
               codon.p = as.numeric( x[1] ) %%3
               codon.p = ifelse( codon.p == 0, 3, codon.p )
               
               R.0 = 
                 length(
                   which(
                     sapply( as.list( seq(1, length(char.i) ) ), 
                             function(x)
                             {
                               codon.2 = codon
                               codon.2[ codon.p ] = char.i[x]
                               
                               return( .check_rep( codon, codon.2 ) )
                             }) == TRUE) )
               
               R.i = R.0/char
             }
             return( R.i* P.i)
           }
    )
  return( sum(R_P) )
}

.s.p   = function( mx, out_conseq = out_conseq )
{
  S_P = 
    apply( mx, 2, 
           function(x)
           {
             
             char.i = grep( pattern = "a|t|c|g", x, ignore.case = TRUE, value = TRUE )
             char   = length( char.i )
             m.i    = char.i[ which( char.i != out_conseq[ as.numeric(x[1]) ] ) ]
             
             if( char < 1)
             {
               P.i = 0
               S.i = 0
             }else
             {
               # P.i - polymorphism score
               
               out.p  = ( char - length( m.i ) )/char
               out.p  = ifelse( ( out.p<1 ) & (out.p>0), 5, out.p  )
               
               type.i = length( unique( m.i ) )
               x.i    = paste0( type.i, out.p )
               
               P.i = prop_site_count_P[ match( x.i, prop_site_count_x ) ]   
               
               
               # S.i - silent score 
               
               codon.n = ceiling( (as.numeric(x[1]))/3 ) - 1
               codon   = out_conseq[ c( codon.n*3+1, codon.n*3+2, codon.n*3+3) ]
               codon.p = as.numeric( x[1] ) %%3
               codon.p = ifelse( codon.p == 0, 3, codon.p )
               
               R.0 = 
                 length(
                   which(
                     sapply( as.list( seq(1, length(char.i) ) ), 
                             function(x)
                             {
                               codon.2 = codon
                               codon.2[ codon.p ] = char.i[x]
                               
                               return( .check_rep( codon, codon.2 ) )
                             }) == TRUE) )
               
               R.i = R.0/char
               S.i = (1 - R.i)
             }
             return( S.i* P.i)
           }
    )
  return( sum(S_P) )
}

.rep.f = function( mx, out_conseq = out_conseq )
{
  R_F = 
    apply( mx, 2, 
           function(x)
           {
             
             char.i = grep( pattern = "a|t|c|g", x, ignore.case = TRUE, value = TRUE )
             char   = length( char.i )
             m.i    = char.i[ which( char.i != out_conseq[ as.numeric(x[1]) ] ) ]
             
             if( char < 1)
             {
               F.i = 0
               R.i = 0
             }else
             {
               # F.i - fixation score
               
               out.p  = ( char - length( m.i ) )/char
               out.p  = ifelse( ( out.p<1 ) & (out.p>0), 5, out.p  )
               
               type.i = length( unique( m.i ) )
               x.i    = paste0( type.i, out.p )
               
               F.i = prop_site_count_F[ match( x.i, prop_site_count_x ) ]   
               
               # R.i - replacement score 
               
               codon.n = ceiling( (as.numeric(x[1]))/3 ) - 1
               codon   = out_conseq[ c( codon.n*3+1, codon.n*3+2, codon.n*3+3) ]
               codon.p = as.numeric( x[1] ) %%3
               codon.p = ifelse( codon.p == 0, 3, codon.p )
               
               R.0 = 
                 length(
                   which(
                     sapply( as.list( seq(1, length(char.i) ) ), 
                             function(x)
                             {
                               codon.2 = codon
                               codon.2[ codon.p ] = char.i[x]
                               
                               return( .check_rep( codon, codon.2 ) )
                             }) == TRUE) )
               
               R.i = R.0/char
             }
             return( R.i* F.i)
           }
    )
  return( sum(R_F) )
}

.s.f   = function( mx, out_conseq = out_conseq )
{
  S_F = 
    apply( mx, 2, 
           function(x)
           {
             
             char.i = grep( pattern = "a|t|c|g", x, ignore.case = TRUE, value = TRUE )
             char   = length( char.i )
             m.i    = char.i[ which( char.i != out_conseq[ as.numeric(x[1]) ] ) ]
             
             if( char < 1)
             {
               F.i = 0
               S.i = 0
             }else
             {
               # F.i - fixation score
               
               out.p  = ( char - length( m.i ) )/char
               out.p  = ifelse( ( out.p<1 ) & (out.p>0), 5, out.p  )
               
               type.i = length( unique( m.i ) )
               x.i    = paste0( type.i, out.p )
               
               F.i = prop_site_count_F[ match( x.i, prop_site_count_x ) ]   
               
               # S.i - silent score 
               
               codon.n = ceiling( (as.numeric(x[1]))/3 ) - 1
               codon   = out_conseq[ c( codon.n*3+1, codon.n*3+2, codon.n*3+3) ]
               codon.p = as.numeric( x[1] ) %%3
               codon.p = ifelse( codon.p == 0, 3, codon.p )
               
               R.0 = 
                 length(
                   which(
                     sapply( as.list( seq(1, length(char.i) ) ), 
                             function(x)
                             {
                               codon.2 = codon
                               codon.2[ codon.p ] = char.i[x]
                               
                               return( .check_rep( codon, codon.2 ) )
                             }) == TRUE) )
               
               R.i = R.0/char
               S.i = ( 1 - R.i )
             }
             return( S.i* F.i )
           }
    )
  return( sum(S_F) )
}





# main function to calculate parameters for adaptive rate

freq_class = function( outgroup_ls,
                       ingroup_ls,
                       def_low_freq = 0.15,
                       def_hi_freq  = 0.75 )
{
  ### the method is adopted from Bhatt et al. (2011), Mol. Biol. Evol., 28(9) ###
  
  require( seqinr )
  
  ingroup  = do.call( rbind, ingroup_ls )
  
  if( is.list(outgroup_ls) )
  {
    outgroup   = do.call( rbind, outgroup_ls )
    out_conseq = 
      apply( outgroup, 2, 
             function(x)
             {
               con = max( grep( pattern = "a|t|c|g", x, ignore.case = TRUE, value = TRUE ) )
               con = ifelse( grepl("a|t|c|g", con, ignore.case = TRUE), con, "-" )
               
               return( con )
             }
             )
    
  }else{ out_conseq = outgroup_ls }
  
  
  if( length(out_conseq) != dim( ingroup )[2] ){ stop( "not consistent seq. length" ) }
  
  ingroup_n = rbind( seq( 1:length(out_conseq) ), ingroup )
  
  # classification 
  freq_mutation = 
    apply( ingroup_n ,
           2,
           function(x)
           {
             char.i = grep( pattern = "a|t|c|g", x, ignore.case = TRUE, value = TRUE )
             n.i    = length( char.i )
             m.i    = length( which( char.i != out_conseq[ as.numeric(x[1]) ] ) )
             
             if( n.i < 1){ p = 0 }else{ p = m.i/n.i }
             return( p )
           })
  
  pos.low = which( freq_mutation < def_low_freq  )
  pos.mid = which( freq_mutation >= def_low_freq & freq_mutation < def_hi_freq  )
  pos.hi  = which( freq_mutation >= def_hi_freq & freq_mutation < 1  )
  pos.fix = which( freq_mutation == 1 )
  
  ingroup_l = ingroup_n[ ,pos.low ]
  ingroup_m = ingroup_n[ ,pos.mid ]
  ingroup_h = ingroup_n[ ,pos.hi ]
  ingroup_f = ingroup_n[ ,pos.fix ]
  
  
  S.f = .s.f( ingroup_f, out_conseq = out_conseq ) 
  R.f = .rep.f( ingroup_f, out_conseq = out_conseq ) 
  
  S = 
  sapply( list( ingroup_h, ingroup_m, ingroup_l ),
          function(x){ .s.p( x, out_conseq = out_conseq )  }
           )
  R = 
  sapply( list( ingroup_h, ingroup_m, ingroup_l ),
          function(x){ .rep.p( x, out_conseq = out_conseq )  }
  )
  
  freq = c( S.f, R.f, S[1], R[1], S[2], R[2], S[3], R[3] )
  names( freq ) = c( "Sf", "Rf", "Sh", "Rh", "Sm", "Rm", "Sl", "Rl"  )
  
  return(freq)
}




# test 

h5_sample = "./gsgd/processed/pH5_c2344_1547ts.fasta"
all       = getSequence( read.fasta( h5_sample ) )
all_id    = attributes( read.fasta( h5_sample ) )$name

all_time  = as.numeric( str_extract( all_id, "([0-9.]+)$" ) )
#range(all_time)

h5_t0 = all[ which( floor( all_time ) < 2008 )  ]
h5_t1 = all[ which( floor( all_time ) < 2011 & floor( all_time ) >= 2008  ) ]
h5_t2 = all[ which( floor( all_time ) < 2014 & floor( all_time ) >= 2011  ) ]  
h5_t3 = all[ which( floor( all_time ) == 2014 )  ]  
h5_t4 = all[ which( floor( all_time ) == 2015 )  ]  
h5_t5 = all[ which( floor( all_time ) == 2016 )  ]  
h5_t6 = all[ which( floor( all_time ) == 2017 )  ]
h5_t7 = all[ which( floor( all_time ) == 2018 )  ]


t1 = freq_class( h5_t0, h5_t1 )
t2 = freq_class( h5_t0, h5_t2 )
t3 = freq_class( h5_t0, h5_t3 )
t4 = freq_class( h5_t0, h5_t4 )
t5 = freq_class( h5_t0, h5_t5 )
t6 = freq_class( h5_t0, h5_t6 )
t7 = freq_class( h5_t0, h5_t7 )




Mi = c( t1[6]/t1[5], t2[6]/t2[5], t3[6]/t3[5], t4[6]/t4[5], t5[6]/t5[5], t6[6]/t6[5], t7[6]/t7[5]  )
Ri = 
sapply( list( t1, t2, t3, t4, t5, t6, t7 ), 
        function(x)
        {
          af = x[2] - x[1]*( x[6]/x[5] )
          ah = x[4] - x[3]*( x[6]/x[5] )
          return(af+ah)
        } 
        )


test_df = data.frame( x = c( 2009, 2 012, 2014, 2015, 2016, 2017, 2018 ),
                      r = Ri )

ggplot( data = test_df, aes( x =  x, y = r) ) + geom_point() + geom_smooth(method = lm)
 



