source("functions.R")


# input --------

rawfas.234 <- "./bigmcc/raw/pH5_c234_2429.fasta"
rawfas.232 <- "./bigmcc/raw/pH5_c232_1296.fasta"
rawfas.7   <- "./bigmcc/raw/pH5_c7_80.fasta"


# light clean --------

rmDup( rawfas.234, year = c( 1996, 3000), rmdup = TRUE )
rmDup( rawfas.232, year = c( 1996, 3000), rmdup = TRUE )
rmDup( rawfas.7, year = c( 1996, 3000), rmdup = TRUE )

rmdup_plus( "./bigmcc/raw/pH5_c234_2429_rmd.fasta" )
rmdup_plus( "./bigmcc/raw/pH5_c232_1296_rmd.fasta" )
rmdup_plus( "./bigmcc/raw/pH5_c7_80_rmd.fasta" )


# geo dist. --------

faslist_234 <- taxaInfo( "./bigmcc/raw/pH5_c234_2429_rmdP.fasta" )
faslist_232 <- taxaInfo( "./bigmcc/raw/pH5_c232_1296_rmdP.fasta" )
faslist_7   <- taxaInfo( "./bigmcc/raw/pH5_c7_80_rmdP.fasta" )

faslist_234[[8]] <- geoID( faslist_234[[2]] )
faslist_232[[8]] <- geoID( faslist_232[[2]] )
faslist_7[[8]]   <- geoID( faslist_7[[2]] )

# remove 1 from 232
# subject to phyml with spr 


# subsampling (1) with geo. as distinct property --------

# clade 7
cladeSampling( trefile   = "./bigmcc/subsampling/ph5_c7_80_rmdp.tre", 
               fasfile   = "./bigmcc/raw/pH5_c7_80_rmdP.fasta", 
               suppList  = TRUE, 
               listinput = faslist_7,
               grid      = 0.5, 
               list.x    = c(6,4,8), 
               saveFasta = TRUE )
# clade 234
cladeSampling( trefile   = "./bigmcc/subsampling/ph5_c234_2429_rmdp.tre", 
               fasfile   = "./bigmcc/raw/pH5_c234_2429_rmdP.fasta", 
               suppList  = TRUE, 
               listinput = faslist_234,
               grid      = 0.5, 
               list.x    = c(6,4,8), 
               saveFasta = TRUE )
# clade 232
cladeSampling( trefile   = "./bigmcc/subsampling/ph5_c232_1296_rmdpe.tre", 
               fasfile   = "./bigmcc/raw/pH5_c232_1296_rmdPe.fasta", 
               suppList  = TRUE, 
               listinput = faslist_232,
               grid      = 0.5, 
               list.x    = c(6,4,8), 
               saveFasta = TRUE )


# prepare BEAST (1) --------

# clade 7
# a trait table
eco.big.7 <- data.frame( id     = fastaEx( "./bigmcc/beast/pH5_c7_80_cs.fasta" )$id, 
                         states = faslist_7[[8]][ match( fastaEx( "./bigmcc/beast/pH5_c7_80_cs.fasta" )$id, faslist_7[[6]] ) ], 
                         stringsAsFactors = FALSE )

eco.big.7$states <- gsub( "^g", "", eco.big.7$states )
write.table( x = eco.big.7, file = "./bigmcc/beast/eco.big.7", sep = "\t", quote = FALSE, row.names = FALSE )

# b state mx 
jumpMx( sort( unique( eco.big.7$states ) ) )

# c time sampling
timeDice( "./bigmcc/beast/pH5_c7_80_cs.fasta", "./bigmcc/beast/eco.big.7", "./time/time_ha_7326" )


# clade 234
# a trait table
eco.big.234 <- data.frame( id     = fastaEx( "./bigmcc/beast/pH5_c234_2429_cs.fasta" )$id, 
                           states = faslist_234[[8]][ match( fastaEx( "./bigmcc/beast/pH5_c234_2429_cs.fasta" )$id, faslist_234[[6]] ) ], 
                           stringsAsFactors = FALSE )

write.table( x = eco.big.234, file = "./bigmcc/beast/eco.big.234", sep = "\t", quote = FALSE, row.names = FALSE )

# b state mx 
jumpMx( sort( unique( eco.big.234$states ) ) )

# c time sampling
timeDice( "./bigmcc/beast/pH5_c234_2429_cs.fasta", "./bigmcc/beast/eco.big.234", "./time/time_ha_7326" )


# clade 232
# a trait table
eco.big.232 <- data.frame( id     = fastaEx( "./bigmcc/beast/pH5_c232_1296_cs.fasta" )$id, 
                           states = faslist_232[[8]][ match( fastaEx( "./bigmcc/beast/pH5_c232_1296_cs.fasta" )$id, faslist_232[[6]] ) ], 
                           stringsAsFactors = FALSE )

write.table( x = eco.big.232, file = "./bigmcc/beast/eco.big.232", sep = "\t", quote = FALSE, row.names = FALSE )

# b state mx 
jumpMx( sort( unique( eco.big.232$states ) ) )

# c time sampling
timeDice( "./bigmcc/beast/pH5_c232_1296_cs.fasta", "./bigmcc/beast/eco.big.232", "./time/time_ha_7326" )



# subsampling (2) --------

faslist_232[[9]] <- rep( "dum", length( faslist_232[[1]]) )
cladeSampling( trefile   = "./bigmcc/subsampling/ph5_c232_1296_rmdpe.tre", 
               fasfile   = "./bigmcc/raw/pH5_c232_1296_rmdPe.fasta", 
               suppList  = TRUE, 
               listinput = faslist_232,
               grid      = 0.5, 
               list.x    = c(6,4,9), 
               saveFasta = TRUE )

faslist_234.tem <- taxaInfo( "./bigmcc/raw/pH5_c232_1296_rmdPe_s.fasta" )
faslist_234.tem[[8]] <- geoID( faslist_234.tem[[2]] )
table( floor(faslist_234.tem[[4]]), faslist_234.tem[[8]] )

eco.big.232t <- data.frame( id     = fastaEx( "./bigmcc/beast/c232_testing/pH5_c232_1296_dum.fasta" )$id, 
                           states = faslist_232[[8]][ match( fastaEx( "./bigmcc/beast/c232_testing/pH5_c232_1296_dum.fasta" )$id, faslist_232[[6]] ) ], 
                           stringsAsFactors = FALSE )

write.table( x = eco.big.232t, file = "./bigmcc/beast/c232_testing/eco.big.232t", sep = "\t", quote = FALSE, row.names = FALSE )

