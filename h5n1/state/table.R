source("functions.R")

require(stringr)
require(car)

# input 

rawtable_G     <- "./data/raw/G_2135_20171124.csv"
eco_table_c232 <- "./state/eco_232.csv"
eco_table_c234 <- "./state/eco_234.csv"

eco_sampling_c232 <- "./dynamics/pH5_c232_204.fasta"
eco_sampling_c234 <- "./dynamics/pH5Nx_c234_240.fasta"

# readin and compile 

# c232

add    <- setdiff( fastaEx( eco_sampling_c232 )$id, read.csv( eco_table_c232, stringsAsFactors = FALSE )$name )
name   <- read.csv( eco_table_c232, stringsAsFactors = FALSE )$name
name   <- c( name, add )

name_c <- str_match( name, pattern = "_([A-Za-z0-9_]+)_" )[,2]

ac     <- str_match( name, pattern = "^[A-Z0-9]+" )
source <- c()

for( i in 1:length(ac) )
{
  if( grepl( pattern = "EPI", ac[i] ) )
  {
    source[ length(source) + 1 ] = "GISAID" 
    temp = gsub( "[A-Z]", "", ac[i] )
    idx  = grep( temp, gsub( "[A-Z]", "", read.csv(rawtable_G, stringsAsFactors = FALSE)$Isolate_Id ) )
    if( length(idx) >1 ) stop() else{ ac[i] = str_match( read.csv(rawtable_G, stringsAsFactors = FALSE)$HA.Segment_Id, "[EPI0-9]+")[idx] }
    
  }else{ source[ length(source) + 1 ] = "GenBank"  }
  
}

eco_species   <- c( read.csv( eco_table_c232, stringsAsFactors = FALSE )$states, "W_a" )
eco_species_c <- recode( eco_species, 
                         "'D_m'='domestic_mammal'; 'D_a'='domestic_avian'; 'D_e'='domestic_environment';
                          'W_m'='wild_mammal'; 'W_a'='wild_avian'; 'W_e'='wild_environment'" )

df232 <- data.frame( accession_no = ac, database = source, name = name_c, eco_species = eco_species_c, clade = "2.3.2" )


# c234

add    <- setdiff( fastaEx( eco_sampling_c234 )$id, read.csv( eco_table_c234, stringsAsFactors = FALSE )$name )
name   <- read.csv( eco_table_c234, stringsAsFactors = FALSE )$name
name   <- c( name, add )

name_c <- str_match( name, pattern = "_([A-Za-z0-9_]+)_" )[,2]

ac     <- str_match( name, pattern = "^[A-Z0-9]+" )
source <- c()

for( i in 1:length(ac) )
{
  if( grepl( pattern = "EPI", ac[i] ) )
  {
    source[ length(source) + 1 ] = "GISAID" 
    temp = gsub( "[A-Z]", "", ac[i] )
    idx  = grep( temp, gsub( "[A-Z]", "", read.csv(rawtable_G, stringsAsFactors = FALSE)$Isolate_Id ) )
    if( length(idx) >1 ) stop() else{ ac[i] = str_match( read.csv(rawtable_G, stringsAsFactors = FALSE)$HA.Segment_Id, "[EPI0-9]+")[idx] }
  }else{ source[ length(source) + 1 ] = "GenBank" }
}

eco_species   <- c( read.csv( eco_table_c234, stringsAsFactors = FALSE )$states, "W_a" )
eco_species_c <- recode( eco_species, 
                         "'D_m'='domestic_mammal'; 'D_a'='domestic_avian'; 'D_e'='domestic_environment';
                          'W_m'='wild_mammal'; 'W_a'='wild_avian'; 'W_e'='wild_environment';
                          'U_e'='Unknown_environment'" )

df234 <- data.frame( accession_no = ac, database = source, name = name_c, eco_species = eco_species_c, clade = "2.3.4" )

write.csv( rbind( df232, df234 ), row.names = FALSE, "Supp_table.csv" )

