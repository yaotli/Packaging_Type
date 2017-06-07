source("~/Packaging_Type/_R/Function.R")
setwd("~/Desktop/Geo/")

library(seqinr)
library(stringr)
library(ggtree)
library(ggplot2)
library(dplyr)
library(plyr)


### classidied countries into different area --------------------------------

# sub_2255 is the old version sub_2553 which remains duplicated sequences

sub2255_file <- read.tree("./tree/sub_2255")

sub2255_geo  <- 
  gsub("\\|", "", 
       str_match(
         gsub("'", "", 
              fortify(sub2255_file)$label[fortify(sub2255_file)$isTip]
              ), "\\|([A-Za-z_])+\\|")[,1])


## manually curate ----------------

countrylist <- sort(unique(sub2255_geo))

# write.csv(countrylist, "geo.csv")

geo_China <- c("China", "Hong_Kong")
geo_NA    <- c("Korea", "South_Korea", "Japan", "Mongolia")
geo_SSEA  <- c("Bangladesh", "Cambodia", "India", "Indonesia", "Laos", "Malaysia", "Myanmar", 
              "Singapore", "Taiwan", "Thailand", "Viet_Nam", "Vietnam")

geo_EU    <- c("Austria", "Azerbaijan", "Belgium", "Bosnia_and_Herzegovina", 
            "Bulgaria", "Croatia", "Czech_Republic", "Denmark", "France", 
            "Germany", "Hungary", "Iceland", "Ireland","Italy", "Netherlands", "Poland", 
            "Portugal", "Romania", "Russia", "Russian_Federation", "Slovakia",
            "Slovenia", "Sweden", "Switzerland", "Ukraine", "United_Kingdom")

geo_ME    <- c("Afghanistan", "Egypt", "Iran", "Iraq", "Israel", "Kuwait", "Lebanon", 
            "Pakistan", "Saudi_Arabia", "Turkey")

geo_NAm   <- c("United_States", "USA")


# check:
# countrylist[which( countrylist %in% c(geo_China, geo_NA, geo_SSEA, geo_EU, geo_ME, geo_NAm) == FALSE)] 



### map Geo info on ML tree --------------------------------

# read-in tree

sub571_file <- read.tree("./tree/sub_571")
sub571_T0   <- ggtree(sub571_file, size = 0.5)

Geo_5       <- c(geo_China, geo_NA, geo_SSEA, geo_EU, geo_ME, geo_NAm)

geo5_color  <- 
  c( rep("gray", length(geo_China)), 
     rep("chartreuse3", length(geo_NA) ), 
     rep("darkorange", length(geo_SSEA) ),
     rep("deepskyblue2", length(geo_EU) ),
     rep("goldenrod", length(geo_ME) ),  
     rep("#9400D3", length(geo_NAm)) )

sub571_color_geo5 <- 
  findtaxa(type = 1, sub571_file, targetid = Geo_5, target = geo5_color)

# basic geo tree 

sub571_geo5 <- 
  sub571_T0 %<+%     
  sub571_color_geo5 + 
  aes(color = I(colorr)) 



### residue and Nx correlation stratified by area --------------------------------

# whole available sequences

file_5p124     <- read.fasta("./align_cut/5p124nt_H5_merged_8602.fasta")
file_5p124_id  <- attributes(file_5p124)$names
file_5p124_seq <- getSequence(file_5p124)

# read-in from tree file

gsgd1883_file    <- read.tree("./tree/sub_1883_GsGD")

taxa_gsgd1883    <- 
  gsub( "'", "", 
        fortify(gsgd1883_file)$label[
          which( fortify(gsgd1883_file)$isTip == TRUE)] )

# prepare data.frame

gsgd1883_5p_seq  <- file_5p124_seq[ match(taxa_gsgd1883, file_5p124_id) ]
gsgd1883_seqMx   <- do.call(rbind, gsgd1883_5p_seq)


gsgd1883_subtype <- gsub("H5", "", 
                         str_match(pattern = "(H[0-9]+N[0-9]+)", taxa_gsgd1883)[,1])

gsgd1883_year    <- 
  
  round(as.numeric(
    gsub("_", "", 
         str_match(taxa_gsgd1883, "_([0-9]{4})\\.([0-9]+)")[,1])), 0)


gsgd1883_country <- 
  
  gsub("\\|", "", 
       str_match(
         gsub("'", "", 
              fortify(gsgd1883_file)$label[fortify(gsgd1883_file)$isTip]
         ), "\\|([A-Za-z_])+\\|")[,1])


# dataframe

df_gsgd1883 <- data.frame(id = taxa_gsgd1883, 
                          subtype = gsgd1883_subtype, 
                          year = gsgd1883_year, 
                          country = gsgd1883_country,
                          p24 = gsgd1883_seqMx[,24],
                          p35 = gsgd1883_seqMx[,35],
                          p72 = gsgd1883_seqMx[,72],
                          p108 = gsgd1883_seqMx[,108], stringsAsFactors = FALSE)
  
# stratified by area

df_gsgd1883_China                  <-
  df_gsgd1883                      %>%
  filter( country %in% geo_China ) %>%
  filter( year >= 2008 )            
  
df_gsgd1883_NA                  <-
  df_gsgd1883                   %>%
  filter( country %in% geo_NA ) %>%
  filter( year >= 2008 )            

df_gsgd1883_SSEA                  <-
  df_gsgd1883                     %>%
  filter( country %in% geo_SSEA ) %>%
  filter( year >= 2008 )            

df_gsgd1883_EU                  <-
  df_gsgd1883                   %>%
  filter( country %in% geo_EU ) %>%
  filter( year >= 2008 )            


### figure strat --------------------------------

# n = 1 - 4 (above 4)


df_area_name = c("df_gsgd1883_China", "df_gsgd1883_NA", "df_gsgd1883_SSEA", "df_gsgd1883_EU")


df_area_data = count(get0(df_area_name[n]), c("year", "subtype") )

if( length(unique(df_area_data$year)) < 10 )
{
   woYr        <- which( c(2008: 2017) %in% unique(df_area_data$year) == FALSE)
   
   tobebinddf  <- data.frame(year = c(2008: 2017)[woYr], 
                            subtype = rep(NA, length(woYr)),
                            freq = rep(0, length(woYr)) )
   
   df_area_data <- rbind(df_area_data, tobebinddf)
}else
{
    
  woYr = NA
  }

  
# p24 

p24_table <- 
  data.frame( table(
    get0( df_area_name[n] ) %>% 
    select(year, p24)  
    ))                      %>% 
    group_by(year)

p24_c <-   
  p24_table          %>%
  filter(p24 == "u") %>%
  select(Freq)

p24_sum <-   
  p24_table          %>%
  filter(p24 != "-") %>%
  dplyr::summarise(sum = sum(Freq))


p24_p    <- p24_c$Freq / p24_sum$sum *80


p24_line      <- data.frame(year = p24_c$year, p24_p, stringsAsFactors = FALSE)
p24_line$year <- factor(p24_line$year, levels = c(2008:2017))

tobebinddf2   <- data.frame(year = c(2008: 2017)[woYr], 
                            p24_p = rep(0, length(woYr) ) )

p24_line      <-  rbind(p24_line, tobebinddf2)

p24_ggplot <- ggplot() + 
  
  geom_bar(data = df_area_data, 
       aes(x = factor(year), y = freq, fill = subtype), stat = "identity")  + 
  theme_bw() + 
  theme(legend.position="none") +
  xlab("")  + ylab("No. Seq") +
  scale_y_continuous(limits = c(0, 80), 
                     sec.axis = sec_axis(~./80, name = "U %") ) + 
  
  scale_fill_manual(values = c("gray" , rep("red", 8 ) ) ) + 
  geom_line( aes(x = as.numeric(year) , y = p24_p), data = p24_line, size = 2) + 
  ggtitle("24")


# p35 


p35_table <- 
  data.frame( table(
    get0( df_area_name[n] ) %>% 
      select(year, p35)  
  ))                        %>% 
  group_by(year)

p35_c <-   
  p35_table          %>%
  filter(p35 == "g") %>%
  select(Freq)

p35_sum <-   
  p35_table          %>%
  filter(p35 != "-") %>%
  dplyr::summarise(sum = sum(Freq))


p35_p    <- p35_c$Freq / p35_sum$sum *80

p35_line <- data.frame(year = p35_c$year, p35_p, stringsAsFactors = FALSE)
p35_line$year <- factor(p35_line$year, levels = c(2008:2017))


tobebinddf2   <- data.frame(year = c(2008: 2017)[woYr], 
                            p35_p = rep(0, length(woYr) ) )

p35_line      <-  rbind(p35_line, tobebinddf2)


p35_ggplot <- ggplot() + 
  
  geom_bar(data = df_area_data, 
           aes(x = as.factor(year), y = freq, fill = subtype), stat = "identity")  + 
  theme_bw() + 
  theme(legend.position="none") +
  xlab("")  + ylab("No. Seq") +
  scale_y_continuous(limits = c(0, 80), 
                     sec.axis = sec_axis(~./80, name = "G %") ) + 
  
  scale_fill_manual(values = c("gray" , rep("red", 8 ) ) ) + 
  geom_line(aes(x = as.numeric(year), y = p35_p), data = p35_line, size = 2) + 
  ggtitle("35")


# p72 


p72_table <- 
  data.frame( table(
    get0( df_area_name[n] ) %>% 
      select(year, p72)  
  ))                        %>% 
  group_by(year)

p72_c <-   
  p72_table          %>%
  filter(p72 == "c") %>%
  select(Freq)

p72_sum <-   
  p72_table          %>%
  filter(p72 != "-") %>%
  dplyr::summarise(sum = sum(Freq))


p72_p    <- p72_c$Freq / p72_sum$sum *80

p72_line <- data.frame(year = p72_c$year, p72_p, stringsAsFactors = FALSE)
p72_line$year <- factor(p72_line$year, levels = c(2008:2017))


tobebinddf2   <- data.frame(year = c(2008: 2017)[woYr], 
                            p72_p = rep(0, length(woYr) ) )

p72_line      <-  rbind(p72_line, tobebinddf2)


p72_ggplot <- ggplot() + 
  
  geom_bar(data = df_area_data, 
           aes(x = as.factor(year), y = freq, fill = subtype), stat = "identity")  + 
  theme_bw() + 
  theme(legend.position="none") +
  xlab("")  + ylab("No. Seq") +
  scale_y_continuous(limits = c(0, 80), 
                     sec.axis = sec_axis(~./80, name = "C %") ) + 
  
  scale_fill_manual(values = c("gray" , rep("red", 8 ) ) ) + 
  geom_line(aes(x = as.numeric(year), y = p72_p), data = p72_line, size = 2) + 
  ggtitle("72")


# p108

p108_table <- 
  data.frame( table(
    get0( df_area_name[n] ) %>% 
      select(year, p108)  
  ))                        %>% 
  group_by(year)

p108_c <-   
  p108_table          %>%
  filter(p108 == "u" | p108 == "a" ) %>%
  dplyr::summarise(sum = sum(Freq))

p108_sum <-   
  p108_table          %>%
  filter(p108 != "-") %>%
  dplyr::summarise(sum = sum(Freq))


p108_p    <- p108_c$sum / p108_sum$sum *80

p108_line <- data.frame(year = p108_c$year, p108_p, stringsAsFactors = FALSE)
p108_line$year <- factor(p108_line$year, levels = c(2008:2017))

tobebinddf2   <- data.frame(year = c(2008: 2017)[woYr], 
                            p108_p = rep(0, length(woYr) ) )

p108_line      <-  rbind(p108_line, tobebinddf2)


p108_ggplot <- ggplot() + 
  
  geom_bar(data = df_area_data, 
           aes(x = as.factor(year), y = freq, fill = subtype), stat = "identity")  + 
  theme_bw() + 
  theme(legend.position="none") +
  xlab("")  + ylab("No. Seq") +
  scale_y_continuous(limits = c(0, 80), 
                     sec.axis = sec_axis(~./80, name = "U or A %") ) + 
  
  scale_fill_manual(values = c("gray" , rep("red", 8 ) ) ) + 
  geom_line(aes(x = as.numeric(year), y = p108_p), data = p108_line, size = 2) + 
  ggtitle("108")


## combine ----------------

multiplot(p24_ggplot, p35_ggplot, p72_ggplot, p108_ggplot, ncol = 1) 



