source("functions.R")

require(seqinr)
require(stringr)
require(ape)
require(ggtree)
require(dplyr)

raw.pH5_c234 <- "./curation/raw/pH5_c234_2429.fasta"
raw.pN1_c234 <- "./curation/raw/pN1_c234_607.fasta"
raw.pH5_c232 <- "./curation/raw/pH5_c232_1296.fasta"
raw.pN1_c232 <- "./curation/raw/pN1_c232_1283.fasta"
    
# keep only China (+HK) seq.

# c234
geo_234       <- read.table("./curation/geo/234rm.txt", header = FALSE, stringsAsFactors = FALSE)
geo_234.id    <- geo_234$V2[ which( geo_234$V1 == "-" ) ]
geo_234.id.na <- fastaEx( raw.pN1_c234 )$id[ match( gsub( "^[A-Za-z0-9]+", "", geo_234.id ), 
                                                    gsub( "^[A-Za-z0-9]+", "", fastaEx( raw.pN1_c234 )$id ) ) ]
       
leafEx( filedir = raw.pH5_c234, setdiff( fastaEx( raw.pH5_c234 )$id, geo_234.id ) )
leafEx( filedir = raw.pN1_c234, setdiff( fastaEx( raw.pN1_c234 )$id, geo_234.id.na ) )
# may need to move file to `curation` folder
rmDup( fasfile = "./curation/pH5_c234_2429_2421.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )
rmDup( fasfile = "./curation/pN1_c234_607_599.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )


# c232
geo_232    <- read.table("./curation/geo/232rm.txt", header = FALSE, stringsAsFactors = FALSE)
geo_232.id <- geo_232$V2[ which( geo_232$V1 == "-" ) ]
geo_232.id.na <- fastaEx( raw.pN1_c232 )$id[ match( gsub( "^[A-Za-z0-9]+", "", geo_232.id ), 
                                                    gsub( "^[A-Za-z0-9]+", "", fastaEx( raw.pN1_c232 )$id ) ) ]

leafEx( filedir = raw.pH5_c232, setdiff( fastaEx( raw.pH5_c232 )$id, geo_232.id ) )
leafEx( filedir = raw.pN1_c232, setdiff( fastaEx( raw.pN1_c232 )$id, geo_232.id.na ) )

rmDup( fasfile = "./curation/pH5_c232_1296_1273.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )
rmDup( fasfile = "./curation/pN1_c232_1283_1260.fasta", sero = "H5N1", geo = c( "China", "Hong_Kong" ), rmdup = FALSE )


# remove reassortants 

# tangleplot 
tre.pH5_7326_c23x  <- read.nexus( "./curation/reassort/pH5_c23x.tre" )
tre.pN1_4468_EA    <- read.nexus( "./curation/reassort/pN1_4468_EA.tre" )

cleantre.pH5_7326_c23x           <-  tre.pH5_7326_c23x
cleantre.pH5_7326_c23x$tip.label <-  sub( pattern = "^'[A-Z0-9]+_", replacement = "'", x = cleantre.pH5_7326_c23x$tip.label )

cleantre.pN1_4468_EA           <-  tre.pN1_4468_EA
cleantre.pN1_4468_EA$tip.label <-  sub( "^'[A-Z0-9]+_", "'", cleantre.pN1_4468_EA$tip.label )

cleantxt_pH5       <- fortify( cleantre.pH5_7326_c23x )
cleantxt_pH5       <- cleantxt_pH5[ which( cleantxt_pH5$isTip ), ]
cleantxt_pH5$label <- gsub( "'", "", cleantxt_pH5$label )
cleantxt_pH5       <- data.frame( cleantxt_pH5,  gene = "H5", stringsAsFactors = FALSE )
cleantxt_pH5$gene[ which( cleantxt_pH5$y >= 1319 ) ] = "h5_234"
cleantxt_pH5$gene[ which( cleantxt_pH5$y < 1319 & cleantxt_pH5$y >= 5 ) ] = "h5_232"

cleantxt_pN1       <- fortify( cleantre.pN1_4468_EA )
cleantxt_pN1       <- cleantxt_pN1[ which( cleantxt_pN1$isTip ), ]
cleantxt_pN1$label <- gsub( "'", "", cleantxt_pN1$label )
cleantxt_pN1       <- data.frame( cleantxt_pN1,  gene = "N1", stringsAsFactors = FALSE )

cleantxt <- rbind( cleantxt_pH5, cleantxt_pN1 )
cleantxt <- data.frame( cleantxt, geo  = str_match( cleantxt$label, "\\|([A-Za-z_]+)\\|" )[, 2],
                        stringsAsFactors = FALSE )

cleantxt %>%
  filter( gene == "h5_232" | gene == "N1" ) %>% 
  #filter( geo == "China" | geo == "Hong_Kong" ) %>% 
  ggplot( aes(x = gene, y = y, group = label) ) + 
  geom_point( size = 0.01 ) +
  geom_line( alpha = 0.5 ) + 
  theme_bw() + ylab("") + xlab("") + 
  geom_hline( yintercept = 3150, color = "blue" ) +
  geom_hline( yintercept = 2658, color = "blue" ) +
  geom_hline( yintercept = 2100, color = "blue" )

# c234
cleantxt_pHA_234 <- cleantxt_pH5[ which( cleantxt_pH5$gene == "h5_234" ), ]
tem.match        <- match( cleantxt_pHA_234$label, cleantxt_pN1$label )
cleantxt_pNA_234 <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
cleantxt_pHA_234 <- cleantxt_pHA_234[ which( cleantxt_pHA_234$label %in% cleantxt_pNA_234$label ), ]
tem.list         <- 
  intersect(  cleantxt_pHA_234[ which( cleantxt_pHA_234$y < 2000 ), ]$label, 
              cleantxt_pNA_234[ which( cleantxt_pNA_234$y < 3250 & cleantxt_pNA_234$y > 2500 ), ]$label )

ls.234   <- grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )
ls.234rm <- setdiff( grep( "\\|China\\||\\|Hong_Kong\\|", cleantxt_pHA_234$label, value = TRUE ), ls.234 )
write( ls.234rm, "./curation/reassort/234rm.txt")

leafEx( raw.pH5_c234, grep( paste0( gsub( "\\|", "\\\\|", ls.234 ), collapse = "|" ), fastaEx( raw.pH5_c234 )$id, value = TRUE ) )
leafEx( raw.pN1_c234, grep( paste0( gsub( "\\|", "\\\\|", ls.234 ), collapse = "|" ), fastaEx( raw.pN1_c234 )$id, value = TRUE ) )


# c232
cleantxt_pHA_232 <- cleantxt_pH5[ which( cleantxt_pH5$gene == "h5_232" ), ]
tem.match        <- match( cleantxt_pHA_232$label, cleantxt_pN1$label )
cleantxt_pNA_232 <- cleantxt_pN1[ tem.match[ !is.na( tem.match ) ], ]
cleantxt_pHA_232 <- cleantxt_pHA_232[ which( cleantxt_pHA_232$label %in% cleantxt_pNA_232$label ), ]
tem.list         <- #2700
  c( cleantxt_pNA_232[ which( cleantxt_pNA_232$y > 3150 ), ]$label, 
     cleantxt_pNA_232[ which( cleantxt_pNA_232$y < 2658 & cleantxt_pNA_232$y > 2100 ), ]$label )

ls.232   <- grep( "\\|China\\||\\|Hong_Kong\\|", tem.list, value = TRUE, ignore.case = TRUE )
ls.232rm <- setdiff( grep( "\\|China\\||\\|Hong_Kong\\|", cleantxt_pHA_232$label, value = TRUE ), ls.232 )
write( ls.232rm, "./curation/reassort/232rm.txt")

leafEx( raw.pH5_c232, grep( paste0( gsub( "\\|", "\\\\|", ls.232 ), collapse = "|" ), fastaEx( raw.pH5_c232 )$id, value = TRUE ) )
leafEx( raw.pN1_c232, grep( paste0( gsub( "\\|", "\\\\|", ls.232 ), collapse = "|" ), fastaEx( raw.pN1_c232 )$id, value = TRUE ) )


# combine and remove duplicated 

faslist  <- paste0( "./curation/", grep( ".fasta$", list.files("./curation/"), value = TRUE ) )
faslist0 <- c( "./curation/raw/pH5_c232_1296.fasta", "./curation/raw/pH5_c234_2429.fasta", 
               "./curation/raw/pN1_c232_1283.fasta", "./curation/raw/pN1_c234_607.fasta" )

for( i in 1:4 )
{
  tem.id <- Reduce( intersect, list( fastaEx(  faslist[ 3*( i-1 ) + 1 ] )$id, fastaEx(  faslist[ 3*( i-1 ) + 2 ] )$id, fastaEx(  faslist[ 3*( i-1 ) + 3 ] )$id) )
  leafEx( filedir = faslist0[ i ] , leaflist = tem.id )
}

faslist2 <- c( "./curation/pH5_c232_1296_231.fasta", "./curation/pN1_c232_1283_231.fasta",
               "./curation/pH5_c234_2429_205.fasta", "./curation/pN1_c234_607_205.fasta" )

for( j in 1: length(faslist2) )
{
  rmDup( faslist2[j]  )
  rmdup_plus( gsub( ".fasta", "_cr.fasta", faslist2[j] ) )
}

# examine and remove redundant sequences with ML trees
# for clade 234, remove the latest isolate due to its time gap to the previous one

leafEx( raw.pH5_c234, tagExtra( "./curation/trees/pH5_c234_2429_205_rmdP_phyml.tre" )$id[ is.na( tagExtra( "./curation/trees/pH5_c234_2429_205_rmdP_phyml.tre" )$tag ) ] )
leafEx( raw.pN1_c234, tagExtra( "./curation/trees/pN1_c234_607_205_rmdP_phyml.tre" )$id[ is.na( tagExtra( "./curation/trees/pN1_c234_607_205_rmdP_phyml.tre" )$tag ) ] )
leafEx( raw.pH5_c232, tagExtra( "./curation/trees/pH5_c232_1296_231_rmdP_phyml.tre" )$id[ is.na( tagExtra( "./curation/trees/pH5_c232_1296_231_rmdP_phyml.tre" )$tag ) ] )










