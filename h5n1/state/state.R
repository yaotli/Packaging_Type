source("functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)



## c234 Eco ----------------

pH5_pri_cn_234      <- taxaInfo( "./curation/raw/pH5_c234_2429.fasta" )
pH5_pri_cn_234[[8]] <- data.frame( name   = pH5_pri_cn_234[[6]][ which( pH5_pri_cn_234[[2]] == "China" | 
                                                                        pH5_pri_cn_234[[2]] == "Hong_Kong")], 
                                   states = NA, ref = NA, stringsAsFactors = FALSE)

pH5_pri_cn_234[[8]]$states <- geoID( strings = pH5_pri_cn_234[[6]][ which( pH5_pri_cn_234[[2]] == "China" | 
                                                                           pH5_pri_cn_234[[2]] == "Hong_Kong")], 
                                     host    = TRUE)

pH5_pri_cn_234[[8]] <- pH5_pri_cn_234[[8]][ sort( str_match( pH5_pri_cn_234[[8]]$name, "^[A-Z0-9]+"), index.return = TRUE )$ix, ]


# human + mammal 
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ][ c( seq(1,25), seq(28, 29), seq(32, 35), seq(37, 46), seq(49, 61), seq(64, 68) ) ] <- 
  "D_m"

pH5_pri_cn_234[[8]]$ref[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ]    <- 
  c( rep("Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", 4 ),
     "Novel Reassortant Avian Influenza A(H5N6) Viruses in Humans, Guangdong, China, 2015*",
     NA, NA, #7
     NA, 
     "GenBank",
     rep( "Highly pathogenic H5N6 influenza A viruses recovered from wild birds in Guangdong, southern China, 2014–2015", 3),
     NA, NA, #14
     "GenBank", #15
     "Genetic Characterization of Continually Evolving Highly Pathogenic H5N6 Influenza Viruses in China, 2012–2016",
     "GenBank",
     NA, NA )
     
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ] <- 
  c( rep("D_m", 2), rep("W_a", 2), 
     "D_e", 
     rep("U_e", 2), #7
     "W_a", 
     "D_m", 
     rep("W_a", 3), #12
     rep("D_m", 2), #14
     "D_m", #15
     "D_a",
     "W_a", #17
     rep("D_m", 2) ) 
  
# environment
pH5_pri_cn_234[[8]]$states[ grep( "environm", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c(
    rep( "D_e", 5 ),
    "D_e",
    rep( "D_e", 5 ), #11
    rep( "D_e", 2 ), #13
    "D_e", #14
    rep("D_e", 68 ), #82
    rep("D_e", 10 ), #92
    rep("D_e", 3 ), #95
    rep("D_e", 57 ), #152
    "D_e", #153
    "D_e", #154
    rep("W_e", 2), #156
    "D_e", #157
    "D_e", #158
    rep( "D_e", 3 ), #161
    rep( "D_e", 7 ), #168
    "W_e", #169
    rep( "D_e", 3),  #172
    rep( "D_e", 3), #175
    "D_e", #176
    "D_e"
  )

pH5_pri_cn_234[[8]]$ref[ grep( "environment", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c(
    rep( "GenBank", 5 ), #5
    "Sustained live poultry market surveillance contributes to early warnings for human infection with avian influenza viruses*",
    rep( "A fatal case of infection with a further reassortant, highly pathogenic avian influenza (HPAI) H5N6 virus in Yunnan, China", 5),
    rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 2), #13
    "Novel Highly Pathogenic Avian H5 Influenza A Viruses in Live Poultry Markets, Wuxi City, China, 2013−2014",
    rep( "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", 68), #82
    rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 10), #92
    rep( "Sustained live poultry market surveillance contributes to early warnings for human infection with avian influenza viruses", 3),
    rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 57), #152
    "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", #153
    "Molecular Evolution and Emergence of H5N6 Avian Influenza Virus in Central China", #154
    rep( "GenBank", 2), #156
    "GenBank", #157
    "Two novel reassortants of avian influenza A (H5N6) virus in China", #158
    rep( "Emergence and evolution of H10 subtype influenza viruses in poultry in China*", 3), #161
    rep( "Continuing Reassortant of H5N6 Subtype Highly Pathogenic Avian Influenza Virus in Guangdong", 7),
    "Dispersal and Transmission of Avian Paramyxovirus Serotype 4 among Wild Birds and Domestic Poultry*", 
    rep( "Emergence of triple-subtype reassortants of fatal human H5N6 avian influenza virus in Yunnan, China", 3),
    rep( "Continuing Reassortant of H5N6 Subtype Highly Pathogenic Avian Influenza Virus in Guangdong", 3),
    "Aerosolized avian influenza A (H5N6) virus isolated from a live poultry market, China", #176
    "GenBank"
  )

pH5_pri_cn_234[[8]]$states[ grep( "water", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c( rep( "D_e", 5 ), rep( "W_a", 2) )

pH5_pri_cn_234[[8]]$ref[ grep( "water", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]    <-
  c( rep( "GenBank", 5), rep( NA, 2 ) )

# confusing species
pH5_pri_cn_234[[8]]$states[ grep( "pigeon|dove", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep( "D_a", 11 ), 
     "D_a",
     rep("D_a", 2 ),
     "W_a", 
     rep("D_a", 6 ) )

pH5_pri_cn_234[[8]]$ref[ grep( "pigeon|dove", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]         <-
  c( rep( "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", 11 ), 
     "GenBank*", 
     rep( "Two novel reassortants of avian influenza A (H5N6) virus in China", 2 ),
     NA, 
     rep("Diversity and evolution of avian influenza viruses in live poultry markets, free-range poultry and wild wetland birds in China", 6 ) )  
    
# poultry 
pH5_pri_cn_234[[8]]$states[ grep( "[A-Za-z]_duck|[A-Za-z]_chicken|[A-Za-z]_goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c( "W_a",
     rep( "W_a", 35 ), #36
     rep( "D_a", 37), #73
     rep( "W_a", 4), #77
     rep( "D_a", 2), #79
     rep( "W_a", 3), #82
     rep( "D_a", 2)
     ) 
    
pH5_pri_cn_234[[8]]$ref[ grep( "[A-Za-z]_duck|[A-Za-z]_chicken|[A-Za-z]_goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c( "GISAID",
     rep( NA, 35 ),
     rep( NA, 37 ),
     rep( NA, 4 ), #77
     rep( "GenBank", 2), #79
     rep( "GenBank", 3), #82
     rep( "GenBank", 2) 
     )
pH5_pri_cn_234[[8]]$states[ grep( "quail|pheasant|peacock", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ intersect( grep( "chicken|duck|goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ), 
                            which( pH5_pri_cn_234[[8]]$states == "nonML" ) )  ]                           <- "D_a"


# wild 
pH5_pri_cn_234[[8]]$states[ grep( "mallard", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- c( "W_a", rep("D_a", 2), "W_a" )
pH5_pri_cn_234[[8]]$ref[ grep( "mallard", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]    <- 
  c( "GISAID",  
     rep("Characterization of duck H5N1 influenza viruses with differing pathogenicity in mallard (Anas platyrhynchos) ducks", 2), 
     "GenBank" )

pH5_pri_cn_234[[8]]$states[ grep( "domestic", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ grep( "wild", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]     <- "W_a"

pH5_pri_cn_234[[8]]$ref[ which( pH5_pri_cn_234[[8]]$states == "nonML") ][ seq(44, 48) ]    <- "GenBank"
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "nonML") ][ seq(44, 48) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "nonML") ]                <- "W_a"


pH5_pri_cn_234[[8]]      <- data.frame( pH5_pri_cn_234[[8]], time = NA)     
pH5_pri_cn_234[[8]]$time <- floor( as.numeric( str_match( pH5_pri_cn_234[[8]]$name, "[0-9.]+$" ) ) )

write.csv( pH5_pri_cn_234[[8]], "eco_234.csv")



## c232 Eco ----------------

pH5_pri_cn_232      <- taxaInfo( "./curation/raw/pH5_c232_1296.fasta" )
pH5_pri_cn_232[[8]] <- data.frame( name   = pH5_pri_cn_232[[6]][ which( pH5_pri_cn_232[[2]] == "China" | 
                                                                        pH5_pri_cn_232[[2]] == "Hong_Kong")], 
                                   states = NA,
                                   ref    = NA, stringsAsFactors = FALSE)

pH5_pri_cn_232[[8]]$states <- geoID( strings = pH5_pri_cn_232[[6]][ which( pH5_pri_cn_232[[2]] == "China" | 
                                                                           pH5_pri_cn_232[[2]] == "Hong_Kong")], 
                                     host    = TRUE)

pH5_pri_cn_232[[8]] <- pH5_pri_cn_232[[8]][ sort( str_match( pH5_pri_cn_232[[8]]$name, "^[A-Z0-9]+"), index.return = TRUE )$ix, ]

# human + mammal
pH5_pri_cn_232[[8]]$states[ which( pH5_pri_cn_232[[8]]$states == "Unknown" ) ] <- 
  c( rep( "D_a", 4 ), 
     rep( "D_m", 3 ), #7
     "W_e", #8
     "D_m", #9
     rep( "W_m", 4 ), #13
     rep( "D_m", 6 ), #19
     rep( "D_a", 2 )
  )

pH5_pri_cn_232[[8]]$ref[ which( pH5_pri_cn_232[[8]]$states == "Unknown" ) ]   <- 
  c( rep( NA, 4 ), 
     rep( NA, 3 ),
     "Highly Pathogenic Avian Influenza A(H5N1) Virus Struck Migratory Birds in China in 2015",
     NA, #9
     rep( "GenBank", 4 ),  #13
     rep( NA, 6 ), #19
     rep( "Characterization of H5N1 highly pathogenic mink influenza viruses in eastern China", 2 ) 
  )

# environment
pH5_pri_cn_232[[8]]$states[ grep( "environment", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "D_e", 
     rep( "D_e", 38 ), # 2-39
     rep( "D_e", 5 ), # 40-44
     rep( "W_e", 4 ), # 45-48
     "W_e", #49
     rep( "D_e", 7) )

pH5_pri_cn_232[[8]]$ref[ grep( "environment", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "Coexistence of influenza H7N9 and H9N2 in poultry linked to human H7N9 infection and their genome characteristics",
     rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 38 ),
     rep( "Novel Highly Pathogenic Avian H5 Influenza A Viruses in Live Poultry Markets, Wuxi City, China, 2013−2014*", 5 ),
     rep( "GenBank", 4 ), #48
     "Full Genome Sequence of an Avian Influenza H5N1 Virus Isolated from the Environment in Hunan Province, China", #49
     rep( "GenBank", 7 ) )

pH5_pri_cn_232[[8]]$ref[ grep( "water", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]    <- "GenBank"
pH5_pri_cn_232[[8]]$states[ grep( "water", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "D_e"



# confusing species
pH5_pri_cn_232[[8]]$states[ grep( "pigeon|dove", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "W_a", 
     "W_a", 
     rep( "D_a", 3) )

pH5_pri_cn_232[[8]]$ref[ grep( "pigeon", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]     <- 
  c( NA,
     "The Survey of H5N1 Flu Virus in Wild Birds in 14 Provinces of China from 2004 to 2007",
     rep( "Genetic and molecular characterization of H9N2 and H5 avian influenza viruses from live poultry markets in Zhejiang Province, eastern China", 3) 
  )


# poultry 
pH5_pri_cn_232[[8]]$states[ grep( "[A-Za-z]_duck|[A-Za-z]_chicken|[A-Za-z]_goose", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "D_a", "W_a", rep("W_a", 5) )
  
pH5_pri_cn_232[[8]]$states[ grep( "quail|pheasant|peacock", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_232[[8]]$states[ intersect( grep( "chicken|duck|goose", pH5_pri_cn_232[[8]]$name, ignore.case = T ), 
                                       which( pH5_pri_cn_232[[8]]$states == "nonML" ) )  ]                <- "D_a"

# wild
pH5_pri_cn_232[[8]]$states[ grep( "wild", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "W_a"
pH5_pri_cn_232[[8]]$states[ which( pH5_pri_cn_232[[8]]$states == "nonML") ]             <- "W_a"


pH5_pri_cn_232[[8]]      <- data.frame( pH5_pri_cn_232[[8]], time = NA)     
pH5_pri_cn_232[[8]]$time <- floor( as.numeric( str_match( pH5_pri_cn_232[[8]]$name, "[0-9.]+$" ) ) )

write.csv( pH5_pri_cn_232[[8]], "eco_232.csv")



## c7 Eco ----------------

pH5_pri_cn_7      <- taxaInfo( "./clade7/pH5_c7_80.fasta" )
pH5_pri_cn_7[[8]] <- data.frame( name = pH5_pri_cn_7[[6]][ which( pH5_pri_cn_7[[2]] == "China" |
                                                                  pH5_pri_cn_7[[2]] == "Hong_Kong" ) ],
                                 states = NA, ref = NA, stringsAsFactors = FALSE )

pH5_pri_cn_7[[8]]$states <- geoID( strings = pH5_pri_cn_7[[6]][ which( pH5_pri_cn_7[[2]] == "China" |
                                                                       pH5_pri_cn_7[[2]] == "Hong_Kong" ) ],
                                   host = TRUE )

pH5_pri_cn_7[[8]] <- pH5_pri_cn_7[[8]][ sort( str_match( pH5_pri_cn_7[[8]]$name, "^[A-Z0-9]+" ), 
                                                         index.return = TRUE)$ix,  ]

# human + mammal 
pH5_pri_cn_7[[8]]$states[ which( pH5_pri_cn_7[[8]]$states == "Unknown") ] <- "D_m"

# environment 
pH5_pri_cn_7[[8]]$ref[ grep("envir", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ] <- 
  c( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 
     "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses*",
     "Novel Highly Pathogenic Avian H5 Influenza A Viruses in Live Poultry Markets, Wuxi City, China, 2013−2014",
     "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 
     rep( "GenBank*", 4) )

pH5_pri_cn_7[[8]]$states[ grep("envir", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ] <- 
  c( rep( "D_e", 4 ), rep( "W_e", 4) )


# poultry
pH5_pri_cn_7[[8]]$states[ grep("chicken|duck|goose|pheasant|peacock", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ] <- "D_a"


# exceptions
pH5_pri_cn_7[[8]]$ref[ grep("mallard", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ] <- "Characterization of duck H5N1 influenza viruses with differing pathogenicity in mallard (Anas platyrhynchos) ducks"
pH5_pri_cn_7[[8]]$states[ grep("mallard", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ] <- "D_a"

# wild
#pH5_pri_cn_7[[8]]$name[ which( pH5_pri_cn_7[[8]]$states == "nonML") ]

pH5_pri_cn_7[[8]]$states[ which( pH5_pri_cn_7[[8]]$states == "nonML") ]                <- "W_a"
pH5_pri_cn_7[[8]]$states[ grep("wild", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ]     <- "W_a"
pH5_pri_cn_7[[8]]$states[ grep("domestic", pH5_pri_cn_7[[8]]$name, ignore.case = T ) ] <- "D_a"

pH5_pri_cn_7[[8]]      <- pH5_pri_cn_7[[8]][ sort( as.numeric( rownames( pH5_pri_cn_7[[8]] ) ), index.return = TRUE )$ix, ]
pH5_pri_cn_7[[8]]      <- data.frame( pH5_pri_cn_7[[8]], time = NA)
pH5_pri_cn_7[[8]]$time <- floor( as.numeric( str_match( pH5_pri_cn_7[[8]]$name, "[0-9.]+$") ) )

write.csv( pH5_pri_cn_7[[8]], "eco_7.csv")