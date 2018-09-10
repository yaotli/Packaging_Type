source("functions.R")

# 7326 h5n1 big fasttree 
ls.gsgd        <- tagExtra( "./clade/raw/pH5_7326_e0823.tre" )
ls.gsgd.wo2344 <- ls.gsgd$id[ which( ls.gsgd$tag == "00ff00" ) ]

leafEx( "./clade/raw/trim_pH5_7326_3.fasta", ls.gsgd.wo2344 )
rmDup( "./clade/raw/trim_pH5_7326_3_4543.fasta", rmdup = TRUE)
rmdup_plus( "./clade/raw/trim_pH5_7326_3_4543_cr.fasta")


# internal validation 

oldc234 <- treeWalk( trefile = "./curation/ref/pH5_c234.tre", showTree = TRUE)
tb.234  <- read.table( "./curation/geo/234rm.txt", header = FALSE, stringsAsFactors = FALSE)
x.234   <- c( setdiff( oldc234[71:81], tb.234$V2), setdiff( tb.234$V2, oldc234[71:81] ))

oldc232 <- treeWalk( trefile = "./curation/ref/pH5_c232.tre" )
tb.232  <- read.table( "./curation/geo/232rm.txt", header = FALSE, stringsAsFactors = FALSE)
x.232   <- c( setdiff( oldc232, tb.232$V2), setdiff( tb.232$V2, oldc232 ))
# most are not to do with geo (same host)


# remove confusing China isolates
reIntro.id  <- treeWalk( trefile = "./clade/raw/fasttree_pH5_7326_3_4543_cr_rmdP.tre", pattern = "00", showTree = TRUE)

rmDup( "./clade/raw/trim_pH5_7326_3_4543_cr_rmdP.fasta", 
       geo  = c( "China", "Hong_Kong"), 
       sero = "H5N1" )

leafEx( "./clade/raw/trim_pH5_7326_3_4543_cr_rmdP_cr.fasta", 
        setdiff( fastaEx( "./clade/raw/trim_pH5_7326_3_4543_cr_rmdP_cr.fasta" )$id, reIntro.id ) )



# remove outlier with TempEst

outs <- c( "FJ602825_bar_headed_goose_Qinghai_10_2008_|China|_H5N1_2008.530",
           "FJ455820_environment_Qinghai_1_2008_|China|_H5N1_2008.445",
           "KX364460_swine_Shandong_SD1_2014_|China|_H5N1_2014.784",
           "FJ602820_bar_headed_goose_Qinghai_5_2008_|China|_H5N1_2008.456",
           "FJ390061_plateau_pika_Qinghai_04_2007_|China|_H5N1_2007.296",
           "KX364468_swine_Shandong_SD2_2014_|China|_H5N1_2014.838",
           "FJ602818_bar_headed_goose_Qinghai_3_2008_|China|_H5N1_2008.445" )

leafEx( "./clade/raw/pH5_737.fasta", 
        setdiff( fastaEx( "./clade/raw/pH5_737.fasta" )$id, outs ) ) # n = 730



cladeSampling( trefile = "./clade/raw/raxml_pH5_730e.tre", 
               fasfile = "./clade/raw/pH5_730.fasta", 
               saveFasta = TRUE, suppList = FALSE, grid = 0.5)



timeDice(fas.dir = "all_h5n1/pH5_730_cs.fasta", ecotable = FALSE, timetab.dir = "./time/time_ha_7326")




