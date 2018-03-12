source("functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)



## c234 Eco ----------------

pH5_pri_cn_234      <- taxaInfo( "./bigtree/raw/pH5_c234_1699_rmdP.fasta" )
pH5_pri_cn_234[[8]] <- data.frame( name   = pH5_pri_cn_234[[6]][ which( pH5_pri_cn_234[[2]] == "China" | 
                                                                        pH5_pri_cn_234[[2]] == "Hong_Kong")], 
                                   states = NA, ref = NA, stringsAsFactors = FALSE)

pH5_pri_cn_234[[8]]$states <- geoID( strings = pH5_pri_cn_234[[6]][ which( pH5_pri_cn_234[[2]] == "China" | 
                                                                           pH5_pri_cn_234[[2]] == "Hong_Kong")], 
                                     host    = TRUE)

pH5_pri_cn_234[[8]] <- pH5_pri_cn_234[[8]][ sort( str_match( pH5_pri_cn_234[[8]]$name, "^[A-Z0-9]+"), index.return = TRUE )$ix, ]

pH5_pri_cn_234[[8]]$ref[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ][c(28, 31, 56, 57, 58, 60)] <-     
  c( "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China",
     "Novel Reassortant Avian Influenza A(H5N6) Viruses in Humans, Guangdong, China, 2015",
     rep("Highly pathogenic H5N6 influenza A viruses recovered from wild birds in Guangdong, southern China, 2014–2015", 3),
     "Genetic Characterization of Continually Evolving Highly Pathogenic H5N6 Influenza Viruses in China, 2012–2016" )
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ] <- 
  c( rep("D_m", 27), "W_a", rep("D_m", 2), "D_e", rep("D_m", 7), rep("U_e", 2), rep("D_m", 15),
     rep("W_a", 3), "D_m", "D_a", rep("D_m", 2) )


pH5_pri_cn_234[[8]]$states[ grep( "environment", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c( 
    rep( "D_e", 3), #1
    "D_e", #2
    rep( "D_e", 4), #3 8
    rep( "D_e", 2), #4 10
    "D_e", #5 11
    rep( "D_e", 35),#6 46
    rep( "D_e", 10), #7 56
    rep( "D_e", 2), #8 58
    rep( "D_e", 51),#9 109
    "D_e", #10
    "D_e", #11
    rep( "W_e", 2), #12 113
    "D_e", #13
    "D_e", #14
    rep( "D_e", 6),#15 121
    "W_e", #16 122
    rep( "D_e", 3), #17
    rep( "D_e", 3), #18
    "D_e", #19
    "D_e" )

pH5_pri_cn_234[[8]]$ref[ grep( "environment", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c( 
    rep( "GenBank", 3 ), #1
    "Sustained live poultry market surveillance contributes to early warnings for human infection with avian influenza viruses", 
    rep( "A fatal case of infection with a further reassortant, highly pathogenic avian influenza (HPAI) H5N6 virus in Yunnan, China", 4),
    rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 2), #4
    "Novel Highly Pathogenic Avian H5 Influenza A Viruses in Live Poultry Markets, Wuxi City, China, 2013−2014",
    rep( "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", 35), #6
    rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 10), #7
    rep( "Sustained live poultry market surveillance contributes to early warnings for human infection with avian influenza viruses", 2),
    rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 51), #9
    "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", #10
    "Molecular Evolution and Emergence of H5N6 Avian Influenza Virus in Central China", #11
    rep( "GenBank", 2), #12
    "GenBank", #13
    "Two novel reassortants of avian influenza A (H5N6) virus in China", #14
    rep( "Continuing Reassortant of H5N6 Subtype Highly Pathogenic Avian Influenza Virus in Guangdong", 6), #15
    "Dispersal and Transmission of Avian Paramyxovirus Serotype 4 among Wild Birds and Domestic Poultry", #16
    rep( "Emergence of triple-subtype reassortants of fatal human H5N6 avian influenza virus in Yunnan, China", 3),
    rep( "Continuing Reassortant of H5N6 Subtype Highly Pathogenic Avian Influenza Virus in Guangdong", 3), #18
    "Aerosolized avian influenza A (H5N6) virus isolated from a live poultry market, China", #19
    "GenBank" )

pH5_pri_cn_234[[8]]$states[ grep( "water", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep( "D_e", 3 ), rep( "W_a", 2 ) )
pH5_pri_cn_234[[8]]$ref[ grep( "water", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]    <-
  c( rep( "GenBank", 3), rep( NA, 2 ) )

pH5_pri_cn_234[[8]]$states[ grep( "pigeon", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep("D_a", 5), "D_a", "D_a", "W_a" )

pH5_pri_cn_234[[8]]$ref[ grep( "pigeon", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep("Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", 5), 
     "Two novel reassortants of avian influenza A (H5N6) virus in China",
     "Complete Genome Sequences of an H5N1 Highly Pathogenic Avian Influenza Virus Isolated from Pigeon in China in 2012",
     "GenBank")

pH5_pri_cn_234[[8]]$states[ grep( "dove", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ][c(2,3)] <- rep("D_a", 2)
pH5_pri_cn_234[[8]]$ref[ grep( "dove", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  rep("Diversity and evolution of avian influenza viruses in live poultry markets, free-range poultry and wild wetland birds in China", 2) 

pH5_pri_cn_234[[8]]$states[ grep( "[A-Za-z]_duck", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- c( rep("D_a", 27), rep("W_a",6) )
pH5_pri_cn_234[[8]]$states[ grep( "[0-9]_duck", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]    <- "D_a"

pH5_pri_cn_234[[8]]$states[ grep( "chicken", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ grep( "quail", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]   <- "D_a"

pH5_pri_cn_234[[8]]$states[ grep( "[A-Za-z]_goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "W_a"
pH5_pri_cn_234[[8]]$states[ grep( "[0-9]_goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]    <- "D_a"

pH5_pri_cn_234[[8]]$ref[ which( pH5_pri_cn_234[[8]]$states == "nonML") ][c(26,27,28)]    <- "GenBank"
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "nonML") ][c(26,27,28)] <- "D_a"

pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "nonML") ]                 <- "W_a"
pH5_pri_cn_234[[8]]$states[ grep( "domestic", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ grep( "wild", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]     <- "W_a"
pH5_pri_cn_234[[8]]$states[ grep( "bar_head", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "W_a"

pH5_pri_cn_234[[8]]      <- pH5_pri_cn_234[[8]][ sort( as.numeric( rownames(pH5_pri_cn_234[[8]]) ), index.return = TRUE )$ix, ]
pH5_pri_cn_234[[8]]      <- data.frame( pH5_pri_cn_234[[8]], time = NA)     
pH5_pri_cn_234[[8]]$time <- floor( as.numeric( str_match( pH5_pri_cn_234[[8]]$name, "[0-9.]+$" ) ) )

write.csv( pH5_pri_cn_234[[8]], "eco_234.csv")



## c232 Eco ----------------

pH5_pri_cn_232      <- taxaInfo( "./bigtree/raw/pH5_c232_1007_rmdP.fasta" )
pH5_pri_cn_232[[8]] <- data.frame( name   = pH5_pri_cn_232[[6]][ which( pH5_pri_cn_232[[2]] == "China" | 
                                                                        pH5_pri_cn_232[[2]] == "Hong_Kong")], 
                                   states = NA,
                                   ref    = NA, stringsAsFactors = FALSE)

pH5_pri_cn_232[[8]]$states <- geoID( strings = pH5_pri_cn_232[[6]][ which( pH5_pri_cn_232[[2]] == "China" | 
                                                                           pH5_pri_cn_232[[2]] == "Hong_Kong")], 
                                     host    = TRUE)

pH5_pri_cn_232[[8]] <- pH5_pri_cn_232[[8]][ sort( str_match( pH5_pri_cn_232[[8]]$name, "^[A-Z0-9]+"), index.return = TRUE )$ix, ]

pH5_pri_cn_232[[8]]$ref[ which( pH5_pri_cn_232[[8]]$states == "Unknown" ) ][c(18,19)] =
  "Characterization of H5N1 highly pathogenic mink influenza viruses in eastern China"

pH5_pri_cn_232[[8]]$states[ which( pH5_pri_cn_232[[8]]$states == "Unknown" ) ] <- 
  c( rep( "D_a", 4), rep( "D_m", 4), rep("W_m", 3), rep("D_m", 8) )


pH5_pri_cn_232[[8]]$states[ grep( "environment", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "D_e", 
     rep( "D_e", 35),  # 2-36
     rep( "D_e", 5), # 37-41
     rep( "W_e", 4), # 42-45
     "W_e", 
     rep( "D_e", 7) )

pH5_pri_cn_232[[8]]$ref[ grep( "environment", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "Coexistence of influenza H7N9 and H9N2 in poultry linked to human H7N9 infection and their genome characteristics",
     rep( "Genesis and dissemination of highly pathogenic H5N6 avian influenza viruses", 35),
     rep( "Novel Highly Pathogenic Avian H5 Influenza A Viruses in Live Poultry Markets, Wuxi City, China, 2013−2014", 5),
     rep( "GenBank", 4 ),
     "Full Genome Sequence of an Avian Influenza H5N1 Virus Isolated from the Environment in Hunan Province, China",
     rep( "GenBank", 7 ) )

pH5_pri_cn_232[[8]]$ref[ grep( "water", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]    <- "GenBank"
pH5_pri_cn_232[[8]]$states[ grep( "water", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "D_e"


pH5_pri_cn_232[[8]]$ref[ grep( "pigeon", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]  <- 
  c( "GenBank", 
     "The Survey of H5N1 Flu Virus in Wild Birds in 14 Provinces of China from 2004 to 2007",
     rep("Genetic and molecular characterization of H9N2 and H5 avian influenza viruses from live poultry markets in Zhejiang Province,eastern China", 2) )
pH5_pri_cn_232[[8]]$states[ grep( "pigeon", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- 
  c( "W_a", "W_a", "D_a", "D_a" )

pH5_pri_cn_232[[8]]$states[ grep( "chicken", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_232[[8]]$states[ grep( "duck", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]    <- "D_a"
pH5_pri_cn_232[[8]]$states[ grep( "quail", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]   <- "D_a"
pH5_pri_cn_232[[8]]$states[ grep( "goose", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]   <- "D_a"
pH5_pri_cn_232[[8]]$states[ grep( "pheasant", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ]<- "D_a"
pH5_pri_cn_232[[8]]$states[ grep( "peacock", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "D_a"

pH5_pri_cn_232[[8]]$states[ which( pH5_pri_cn_232[[8]]$states == "nonML") ]             <- "W_a"
pH5_pri_cn_232[[8]]$states[ grep( "wild", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "W_a"
pH5_pri_cn_232[[8]]$states[ grep( "bar_head", pH5_pri_cn_232[[8]]$name, ignore.case = T ) ] <- "W_a"

pH5_pri_cn_232[[8]] <- pH5_pri_cn_232[[8]][ sort( as.numeric( rownames(pH5_pri_cn_232[[8]]) ), index.return = TRUE )$ix, ]

pH5_pri_cn_232[[8]]      <- data.frame( pH5_pri_cn_232[[8]], time = NA)     
pH5_pri_cn_232[[8]]$time <- floor( as.numeric( str_match( pH5_pri_cn_232[[8]]$name, "[0-9.]+$" ) ) )

write.csv( pH5_pri_cn_232[[8]], "eco_232.csv")

