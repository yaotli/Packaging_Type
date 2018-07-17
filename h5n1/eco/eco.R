source("./functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(ggpubr)


raw_pH5_c234   <- "./dynamics/pH5_c234_159.fasta"
raw_pN1_c234   <- "./dynamics/pN1_c234_152.fasta"
raw_pH5_c232   <- "./dynamics/pH5_c232_204.fasta"
raw_pN1_c232   <- "./dynamics/pN1_c232_189.fasta"

raw_pH5Nx_c234 <- "./dynamics/pH5Nx_c234_240.fasta"

csv_eco_c234   <- "./state/eco_234.csv"
csv_eco_c232   <- "./state/eco_232.csv"

time_ha        <- "./time/time_ha_7326"
time_na        <- "./time/time_na_7138"


# c234 --------
sub_pH5_c234_159        <- taxaInfo( file = raw_pH5_c234, useTree = FALSE )
sub_pH5_c234_159[[ 8 ]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5_c234_159[[6]], read.csv( csv_eco_c234, stringsAsFactors = FALSE )$name ) ]

sub_pH5_c234_159[[ 8 ]] <- ifelse( startsWith( sub_pH5_c234_159[[ 8 ]], "D" ), "D", "W" )
eco.234 <- data.frame( id = sub_pH5_c234_159[[ 6 ]], states = sub_pH5_c234_159[[ 8 ]], stringsAsFactors = FALSE )
write.table( x = eco.234, file = "./eco/eco.234", sep = "\t", quote = FALSE, row.names = FALSE )

#c234_n1
sub_pN1_c234_152      <- taxaInfo( file = raw_pN1_c234, useTree = FALSE )
sub_pN1_c234_152[[8]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE)$states[ 
  match( gsub( pattern = "^[A-Z0-9]+", "", sub_pN1_c234_152[[6]] ), 
         gsub( pattern = "^[A-Z0-9]+", "", read.csv( csv_eco_c234, stringsAsFactors = FALSE)$name ) ) ]

sub_pN1_c234_152[[8]] <- ifelse( startsWith( sub_pN1_c234_152[[8]], "D"), "D", "W")
eco.234.n1            <- data.frame( id = sub_pN1_c234_152[[6]], states = sub_pN1_c234_152[[8]], stringsAsFactors = FALSE )
write.table( x = eco.234.n1, file = "./eco/eco.234.n1", sep = "\t", quote = FALSE, row.names = FALSE )

#c234_h5nx
sub_pH5Nx_c234_240        <- taxaInfo( file = raw_pH5Nx_c234, useTree = FALSE )
sub_pH5Nx_c234_240[[ 8 ]] <- read.csv( csv_eco_c234, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5Nx_c234_240[[6]], read.csv( csv_eco_c234, stringsAsFactors = FALSE )$name ) ]

sub_pH5Nx_c234_240[[ 8 ]] <- ifelse( startsWith( sub_pH5Nx_c234_240[[ 8 ]], "D" ), "D", "W" )
eco.234nx <- data.frame( id = sub_pH5Nx_c234_240[[ 6 ]], states = sub_pH5Nx_c234_240[[ 8 ]], stringsAsFactors = FALSE )
write.table( x = eco.234nx, file = "./eco/eco.234nx", sep = "\t", quote = FALSE, row.names = FALSE )

# eco-stratified 
ac_c234h5_D <- acSearch( faslist = sub_pH5_c234_159, keyword.dir = 8, keyword = "D")
ac_c234h5_W <- acSearch( faslist = sub_pH5_c234_159, keyword.dir = 8, keyword = "W")
subfastaSeq( AC = TRUE, AC_list = ac_c234h5_D, filedir = raw_pH5_c234 )
subfastaSeq( AC = TRUE, AC_list = ac_c234h5_W, filedir = raw_pH5_c234 )





#c232 --------
sub_pH5_c232_204        <- taxaInfo( file = raw_pH5_c232, useTree = FALSE )
sub_pH5_c232_204[[ 8 ]] <- read.csv( csv_eco_c232, stringsAsFactors = FALSE )$states[ 
  match( sub_pH5_c232_204[[6]], read.csv( csv_eco_c232, stringsAsFactors = FALSE )$name ) ]

sub_pH5_c232_204[[ 8 ]] <- ifelse( startsWith( sub_pH5_c232_204[[ 8 ]], "D" ), "D", "W" )

eco.232 <- data.frame( id = sub_pH5_c232_204[[ 6 ]], states = sub_pH5_c232_204[[ 8 ]], stringsAsFactors = FALSE )
write.table( x = eco.232, file = "./eco/eco.232", sep = "\t", quote = FALSE, row.names = FALSE )

eco.232_2014 <- eco.232[ which( floor( as.numeric( str_match( eco.232$id, "[0-9.]+$" ) ) ) < 2014 ), ]
write.table( x = eco.232_2014, file = "./eco/eco.232_2014", sep = "\t", quote = FALSE, row.names = FALSE )
#eco.232_2012 <- eco.232[ which( floor( as.numeric( str_match( eco.232$id, "[0-9.]+$" ) ) ) < 2012 ), ]
#write.table( x = eco.232_2012, file = "./eco/eco.232_2012", sep = "\t", quote = FALSE, row.names = FALSE )


#c232_n1
sub_pN1_c232_189      <- taxaInfo( file = raw_pN1_c232, useTree = FALSE )
sub_pN1_c232_189[[8]] <- read.csv( csv_eco_c232, stringsAsFactors = FALSE)$states[ 
  match( gsub( pattern = "^[A-Z0-9]+", "", sub_pN1_c232_189[[6]] ), 
         gsub( pattern = "^[A-Z0-9]+", "", read.csv( csv_eco_c232, stringsAsFactors = FALSE)$name ) ) ]

sub_pN1_c232_189[[8]] <- ifelse( startsWith( sub_pN1_c232_189[[8]], "D"), "D", "W")
eco.232.n1            <- data.frame( id = sub_pN1_c232_189[[6]], states = sub_pN1_c232_189[[8]], stringsAsFactors = FALSE )
write.table( x = eco.232.n1, file = "./eco/eco.232.n1", sep = "\t", quote = FALSE, row.names = FALSE )


# eco-stratified c232
rmDup( fasfile = raw_pH5_c232, year = c( 2000, 2014 ), rmdup = FALSE )
rmDup( fasfile = raw_pN1_c232, year = c( 2000, 2014 ), rmdup = FALSE )

ac_c232h5_D <- acSearch( faslist = sub_pH5_c232_204, keyword.dir = 8, keyword = "D", range = c( 2000, 2011 ), range.dir = 4 )
ac_c232h5_W <- acSearch( faslist = sub_pH5_c232_204, keyword.dir = 8, keyword = "W", range = c( 2000, 2011 ), range.dir = 4 )
subfastaSeq( AC = TRUE, AC_list = ac_c232h5_D, filedir = "./dynamics/ride_2012/pH5_c232_204_2012.fasta" )
subfastaSeq( AC = TRUE, AC_list = ac_c232h5_W, filedir = "./dynamics/ride_2012/pH5_c232_204_2012.fasta" )





# time-sampling

timeDice( "./dynamics/pH5Nx_c234_240.fasta", "./eco/eco.234nx", "time/time_ha_7326" )

timeDice( fas.dir = "./Nx/pH5_c234_81.fasta", ecotab.dir = "./eco/eco.234nx", timetab.dir = "time/time_ha_7326",ecotable = TRUE)
timeDice( fas.dir = "./clade7/Nx/pH5nx_c7_65.fasta", timetab.dir = "time/time_ha_7326",ecotable = FALSE)


timeDice( "./dynamics/pH5_c234_159.fasta", "./eco/eco.234", "./time/time_ha_7326" )
timeDice( "./dynamics/ride_2012/pH5_c232_204_2012.fasta", "./eco/eco.232", "./time/time_ha_7326" )
timeDice( "./dynamics/pN1_c234_152.fasta", "./eco/eco.234.n1", "./time/time_na_7138" )
timeDice( "./dynamics/ride_2012/pN1_c232_189_2012.fasta", "./eco/eco.232.n1", "./time/time_na_7138" )

timeDice( "./dynamics/pH5Nx_c234_240.fasta", "./eco/eco.234nx", "./time/time_ha_7326" )
timeDice( "./dynamics/ride_2014/pH5_c232_135_2014.fasta", "./eco/eco.232_2014", "./time/time_ha_7326" )

timeDice( "./eco/pH5_c232_204_2012D.fasta", "./eco/eco.232", "./time/time_ha_7326" )
timeDice( "./eco/pH5_c232_204_2012W.fasta", "./eco/eco.232", "./time/time_ha_7326" )
timeDice( "./eco/pH5_c234_159D.fasta", "./eco/eco.234", "./time/time_ha_7326" )
timeDice( "./eco/pH5_c234_159W.fasta", "./eco/eco.234", "./time/time_ha_7326" )




