library(seqinr)
library(dplyr)
library(stringr)
library(ggtree)
library(ape)

setwd("~/Desktop/data_souce/Nov2017_hana/")
source("~/Packaging_Type/_R/Function.R")


### 232 primary data (CN_HK) Species --------------------------------

rmDup( fasfile = "./c232/rmd/pH5_c232_1007_rmdP.fasta", geo = c( "China", "Hong_Kong" ) )




pH5_pri_cn_232      <- taxaInfo( "./c232/rmd/cn/pH5_c232_232.fasta" )
pH5_pri_cn_232[[8]] <- data.frame( name   = pH5_pri_cn_232[[6]], 
                                   states = geoID( strings = pH5_pri_cn_232[[6]], host = TRUE ),
                                   ref    = NA, stringsAsFactors = FALSE )

pH5_pri_cn_232[[8]] <- pH5_pri_cn_232[[8]][ sort( str_match( pH5_pri_cn_232[[8]]$name, "^[A-Z0-9]+"), index.return = TRUE )$ix, ]

pH5_pri_cn_232[[8]]$ref[ which( pH5_pri_cn_232[[8]]$states == "Unknown" ) ][c(18,19)] =
  "Characterization of H5N1 highly pathogenic mink influenza viruses in eastern China"

pH5_pri_cn_232[[8]]$states[ which( pH5_pri_cn_232[[8]]$states == "Unknown" ) ] <- 
  c( rep( "D_a", 4), rep( "D_m", 4), rep("W_m", 3), rep("D_m", 8) )



# environment
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


# other avian 
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

# write.csv( pH5_pri_cn_232[[8]], "eco_232.csv")

### 234 primary data (CN_HK) Species  --------------------------------

rmDup( fasfile = "./c234/rmd/pH5_c234_1699_rmdP.fasta", geo = c( "China", "Hong_Kong" ) )




pH5_pri_cn_234      <- taxaInfo( "./c234/rmd/cn/pH5_c234_858.fasta" )
pH5_pri_cn_234[[8]] <- data.frame( name   = pH5_pri_cn_234[[6]], 
                                   states = geoID( strings = pH5_pri_cn_234[[6]], host = TRUE ),
                                   ref    = NA, stringsAsFactors = FALSE)
pH5_pri_cn_234[[8]] <- pH5_pri_cn_234[[8]][ sort( str_match( pH5_pri_cn_234[[8]]$name, "^[A-Z0-9]+"), index.return = TRUE )$ix, ]

pH5_pri_cn_234[[8]]$ref[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ][c(28, 31, 56, 57, 58, 60)] <-     
  c( "Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China",
     "Novel Reassortant Avian Influenza A(H5N6) Viruses in Humans, Guangdong, China, 2015",
     rep("Highly pathogenic H5N6 in uenza A viruses recovered from wild birds in Guangdong, southern China, 2014–2015", 3),
     "Genetic Characterization of Continually Evolving Highly Pathogenic H5N6 Influenza Viruses in China, 2012–2016" )
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "Unknown" ) ] <- 
  c( rep("D_m", 27), "W_a", rep("D_m", 2), "D_e", rep("D_m", 7), rep("U_e", 2), rep("D_m", 15),
     rep("W_a", 3), "D_m", "D_a", rep("D_m", 2) )
     
  
# environment 
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
     
     
# water      
pH5_pri_cn_234[[8]]$states[ grep( "water", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep( "D_e", 3), rep( "W_a", 2) )
pH5_pri_cn_234[[8]]$ref[ grep( "water", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep( "GenBank", 3), rep(NA, 2) )
     
# other avian
pH5_pri_cn_234[[8]]$states[ grep( "pigeon", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep("D_a", 5),
     "D_a",
     "D_a",
     "W_a" )
     
pH5_pri_cn_234[[8]]$ref[ grep( "pigeon", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <-
  c( rep("Genesis, Evolution and Prevalence of H5N6 Avian Influenza Viruses in China", 5), 
     "Two novel reassortants of avian influenza A (H5N6) virus in China",
     "Complete Genome Sequences of an H5N1 Highly Pathogenic Avian Influenza Virus Isolated from Pigeon in China in 2012",
     "GenBank")

pH5_pri_cn_234[[8]]$states[ grep( "[A-Za-z]_duck", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- 
  c( rep("D_a", 27), rep("W_a",6) )
pH5_pri_cn_234[[8]]$states[ grep( "[0-9]_duck", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"

pH5_pri_cn_234[[8]]$states[ grep( "chicken", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ grep( "quail", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]   <- "D_a"

pH5_pri_cn_234[[8]]$states[ grep( "[A-Za-z]_goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "W_a"
pH5_pri_cn_234[[8]]$states[ grep( "[0-9]_goose", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]    <- "D_a"

# 'avian'
pH5_pri_cn_234[[8]]$ref[ which( pH5_pri_cn_234[[8]]$states == "nonML") ][c(31, 32, 33)]    <- "GenBank"
pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "nonML") ][c(31, 32, 33)] <- "D_a"

pH5_pri_cn_234[[8]]$states[ which( pH5_pri_cn_234[[8]]$states == "nonML") ]                 <- "W_a"
pH5_pri_cn_234[[8]]$states[ grep( "domestic", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "D_a"
pH5_pri_cn_234[[8]]$states[ grep( "wild", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ]     <- "W_a"
pH5_pri_cn_234[[8]]$states[ grep( "bar_head", pH5_pri_cn_234[[8]]$name, ignore.case = T ) ] <- "W_a"


pH5_pri_cn_234[[8]]      <- pH5_pri_cn_234[[8]][ sort( as.numeric( rownames(pH5_pri_cn_234[[8]]) ), index.return = TRUE )$ix, ]

pH5_pri_cn_234[[8]]      <- data.frame( pH5_pri_cn_234[[8]], time = NA)     
pH5_pri_cn_234[[8]]$time <- floor( as.numeric( str_match( pH5_pri_cn_234[[8]]$name, "[0-9.]+$" ) ) )

# write.csv( pH5_pri_cn_234[[8]], "eco_234.csv")

    
### figure species --------------------------------


# 232
table_232_species        <- as.data.frame( 
  prop.table( table( pH5_pri_cn_232[[8]][,2], pH5_pri_cn_232[[8]][,4]), margin = 2) )


table_232_species$Var1 <- factor( table_232_species$Var1, levels = c("D_a", "D_m", "D_e", "W_e", "W_m", "W_a") )
table_232_species$Var2 <- as.numeric( as.character(table_232_species$Var2) )

s1 <- 
ggplot( data =  table_232_species, 
        aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2002.5, 2016.5), breaks = seq(2003, 2016), labels = seq(2003, 2016),
                      expand = c(0,0.1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Sequences") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e") ) + 
  theme( legend.title = element_blank(), 
         axis.title.y = element_text(size = 16) )

  
     
# 234
table_234_species        <- as.data.frame( 
  prop.table( table( pH5_pri_cn_234[[8]][,2], pH5_pri_cn_234[[8]][,4]), margin = 2) )

table_234_species$Var1 <- factor( table_234_species$Var1, levels = c("D_a", "D_m", "D_e", "U_e","W_e", "W_a") )
table_234_species$Var2 <- as.numeric( as.character(table_234_species$Var2) )

s2 <- 
ggplot( data =  table_234_species, 
        aes( x = Var2, y = Freq, fill = Var1) ) +
  geom_bar(stat = "identity") + 
  scale_x_continuous( limit = c(2002.5, 2016.5), breaks = seq(2003, 2016), labels = seq(2003, 2016),
                      expand = c(0,0.1)) +
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Sequences") +
  scale_fill_manual( values = c( "#8c510a", "#d8b365", "#f6e8c3", "white", "#c7eae5", "#01665e") ) +
  theme( legend.title = element_blank(), 
         axis.title.y = element_text(size = 16) )

multiplot( s2, s1, ncol = 1)


     
### figure counts (jump) --------------------------------

jumpcount <- data.frame( time  = rep( c("0405", "0405", "0607", "0607", "0809", "0809", "1011", "1011"), 2), 
                         clade = c( rep( "232", 8), rep( "234", 8) ),
                         count = c( 1, 0, 2, 3, 4, 8, 14, 13, 6, 4, 9, 2, 7, 8, NA, NA),
                         type  = rep( c("DtoW", "WtoD"), 8) )     
jumpcount <- jumpcount[-c(15, 16), ]
jumpcount$clade <- factor( jumpcount$clade, levels = c("234", "232") )

j1 <- 
ggplot( data = jumpcount[c(9:16), ] ) + 
  geom_line( aes( x = time, y = count, color = type, group = type ), size = 2) + 
  geom_point( aes( x = time, y = count, color = type, group = type ), size = 4) +
  theme_bw() + 
  xlab("") + ylab("Count") +
  scale_x_discrete(labels = c("2005", "2007", "2009", "2011") ) +
  scale_y_continuous( limits = c(0,15), breaks = seq(0,14, by = 2), labels = seq(0,14, by = 2) ) +
  scale_color_manual(values = c( "#8c564b", "#17becf" ) ) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.title       = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 12),
        legend.position = "none" )

j2 <- 
  ggplot( data = jumpcount[c(1:8), ] ) + 
  geom_line( aes( x = time, y = count, color = type, group = type ), size = 2) + 
  geom_point( aes( x = time, y = count, color = type, group = type ), size = 4) +
  theme_bw() + 
  xlab("") + ylab("Count") +
  scale_x_discrete(labels = c("2005", "2007", "2009", "2011") ) +
  scale_y_continuous( limits = c(0,15), breaks = seq(0,14, by = 2), labels = seq(0,14, by = 2) ) +
  scale_color_manual(values = c( "#8c564b", "#17becf" ) ) +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.major.x = element_blank(),
        axis.title       = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none" )


multiplot( j1, j2, ncol = 1)



### figure geo --------------------------------

# 232 
pH5_pri_cn_232[[9]] <- geoID( pH5_pri_cn_232[[6]] )
pH5_pri_cn_232[[9]][ which( pH5_pri_cn_232[[9]] == "Unknown") ] <- 
  c( rep(NA, 3), rep(NA, 2), "cnS", rep("cnC", 3), "cnS")

# inner_mongolia
pH5_pri_cn_232[[9]][ which(pH5_pri_cn_232[[9]] == "nA") ]   <- "cnN"
pH5_pri_cn_232[[9]][ which(pH5_pri_cn_232[[9]] == "cnBH") ] <- "cnE"
pH5_pri_cn_232[[9]][ which(pH5_pri_cn_232[[9]] == "cnYZ") ] <- "cnE"
pH5_pri_cn_232[[9]][ which(pH5_pri_cn_232[[9]] == "cnNE" | 
                           pH5_pri_cn_232[[9]] == "cnNW") ] <- "cnN"

pH5_pri_cn_232[[9]]      <- data.frame( geo = pH5_pri_cn_232[[9]], time = NA)     
pH5_pri_cn_232[[9]]$time <- floor( as.numeric( pH5_pri_cn_232[[4]] ) )
pH5_pri_cn_232[[9]]$geo  <- factor( pH5_pri_cn_232[[9]]$geo, levels = c("cnN", "cnE", "cnC", "cnSW", "cnS") )

table_232_geo        <- as.data.frame( prop.table( table( pH5_pri_cn_232[[9]] ), margin = 2) )
table_232_geo$time   <- as.numeric( as.character(table_232_geo$time) )


# 234 
pH5_pri_cn_234[[9]] <- geoID( pH5_pri_cn_234[[6]] )
pH5_pri_cn_234[[9]][ which( pH5_pri_cn_234[[9]] == "Unknown") ] <- 
  c( rep("cnE", 8), "cnN", rep("cnE", 25), "cnS")

pH5_pri_cn_234[[9]][ which(pH5_pri_cn_234[[9]] == "cnBH") ] <- "cnE"
pH5_pri_cn_234[[9]][ which(pH5_pri_cn_234[[9]] == "cnYZ") ] <- "cnE"
pH5_pri_cn_234[[9]][ which(pH5_pri_cn_234[[9]] == "cnNE" | 
                           pH5_pri_cn_234[[9]] == "cnNW") ] <- "cnN"

pH5_pri_cn_234[[9]]      <- data.frame( geo = pH5_pri_cn_234[[9]], time = NA)     
pH5_pri_cn_234[[9]]$time <- floor( as.numeric( pH5_pri_cn_234[[4]] ) )
pH5_pri_cn_234[[9]]$geo  <- factor( pH5_pri_cn_234[[9]]$geo, levels = c("cnN", "cnE", "cnC", "cnSW", "cnS") )

table_234_geo        <- as.data.frame( prop.table( table( pH5_pri_cn_234[[9]] ), margin = 2) )
table_234_geo$time   <- as.numeric( as.character(table_234_geo$time) )


# figures 

g1 <- 
ggplot( data =  table_232_geo, 
        aes( x = time, y = Freq, fill = geo) ) +
  scale_x_continuous( limit = c(2002.5, 2016.5), breaks = seq(2003, 2016), labels = seq(2003, 2016),
                      expand = c(0, 0.1)) +
  geom_bar(stat = "identity") + 
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Sequences") +
  scale_fill_manual( values = pyCol( c( "gray", "green", "brown", "blue", "red") ),
                     breaks = c("cnN", "cnE", "cnC", "cnSW", "cnS"), 
                     labels = c("N", "E", "C", "SW", "S") ) + 
  theme( legend.title = element_blank(), 
         axis.title.y = element_text(size = 16) )

g2 <- 
ggplot( data =  table_234_geo, 
        aes( x = time, y = Freq, fill = geo) ) +
  scale_x_continuous( limit = c(2002.5, 2016.5), breaks = seq(2003, 2016), labels = seq(2003, 2016),
                      expand = c(0,0.1)) +
  geom_bar(stat = "identity") + 
  theme( panel.background = element_blank() ) +
  xlab("") + ylab("% Sequences") +
  scale_fill_manual( values = pyCol( c( "gray", "green", "brown", "blue", "red") ),
                     breaks = c("cnN", "cnE", "cnC", "cnSW", "cnS"), 
                     labels = c("N", "E", "C", "SW", "S") ) + 
  theme( legend.title = element_blank(), 
         axis.title.y = element_text(size = 16) )

multiplot(g2, g1, ncol = 1) #4*6in

