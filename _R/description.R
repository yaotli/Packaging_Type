# aim: the description statistics of position of interests
#
# target position: (5') 35 and 108
# 

source("./_R/Function.R")
library(stringr)
library(seqinr)

    gsgd_fasta <- read.fasta("~/Desktop/RNA_GsGD.fas")
non_gsdd_fasta <- read.fasta("~/Desktop/RNA_nonGsGD.fas")

  seq_name0 = attributes(gsgd_fasta)$names
       seq0 = getSequence(gsgd_fasta)

  seq_name0_ng = attributes(non_gsdd_fasta)$names
       seq0_ng = getSequence(non_gsdd_fasta)
       
       
  # eliminate unsubtyped virus     

  # gsgd
  nonsubtype <- which(is.na(str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name0)[,1]) == TRUE)
  nonsubtype <- unique( sort( c(nonsubtype, grep("H5N0", x = seq_name0)) ) )

  seq_name <- seq_name0[-nonsubtype]
       seq <- seq0[-nonsubtype]
       
  # nongsgd     
  nonsubtype_ng <- which(is.na(str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name0_ng)[,1]) == TRUE)
  nonsubtype_ng <- unique( sort( c(nonsubtype_ng, grep("H5N0", x = seq_name0_ng)) ) )
       
  seq_name_ng <- seq_name0_ng[-nonsubtype_ng]
  seq_ng <- seq0_ng[-nonsubtype_ng]
      
  # input data        

     seq_matrix <- do.call(rbind, seq)
  seq_matrix_ng <- do.call(rbind, seq_ng)
  
  # extract subtype data  

     n1_id <- which( str_match(pattern = "_(H[0-9]+N[0-9]+)_", seq_name)[,2] == "H5N1" )
  nonN1_id <- seq(1:length(seq_name))[-n1_id]
  
  
# summary ####  
  
  X = 1735
  
  sub1651 <- seq_matrix[, X]
   id1651 <- seq_name[-which(sub1651 == "-")]   
  sub1651 <- sub1651[-which(sub1651 == "-")]
  
   NA1651 <- gsub(pattern = "H5", 
                  replacement = "", 
                  str_match(pattern = "_(H[0-9]+N[0-9]+)_", 
                            id1651)[,2])
   
   NA1651_N1 <- rep("nonN1", length(sub1651))
   NA1651_N1[which(NA1651 == "N1")]  = "N1"
  
  year1651 <- as.numeric( 
    gsub("_", "", str_match(pattern = "_([0-9]{4})\\.([0-9]+)", id1651)[,1]) 
    )
  
  # gsgd
  time1651 <- rep("2002", length(year1651))
  time1651[which(year1651 < 2018)] = "2017"
  time1651[which(year1651 < 2013)] = "2012"
  time1651[which(year1651 < 2008)] = "2007"
  time1651[which(year1651 < 2003)] = "2002"
  
    
  # nongsgd
  
  sub1651_ng <- seq_matrix_ng[, X]
   id1651_ng <- seq_name_ng[-which(sub1651_ng == "-")]
  sub1651_ng <- sub1651_ng[-which(sub1651_ng == "-")]
   
  
  # df   
  df_1651 = data.frame(id = id1651, 
                       nt = sub1651, 
                       Nx = NA1651, 
                       N1 = NA1651_N1,  
                       year = year1651)
  
  df_1651_1 = as.data.frame(table(df_1651$nt, df_1651$N1))
  
  df_1651_2 = as.data.frame(table(time1651, sub1651))
  
  df_1651_3 = as.data.frame(table(sub1651_ng))
  # df_1651_3 = df_1651_3[-4,] 
  
  df_1651_4 = as.data.frame(table(time1651, NA1651))
  
  
  # figure
p1 =  
ggplot(data = df_1651_1, aes(x = Var2, y = Freq, fill = Var1))  + 
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  guides(fill = guide_legend(title=NULL)) +
  xlab("") + 
  labs(title = "GsGD") + theme(legend.position="top")

p2 =   
ggplot(data = df_1651_2, aes(x = time1651, y = Freq, fill = sub1651))  + 
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  guides(fill = guide_legend(title=NULL)) +
  xlab("") + 
  labs(title = "GsGD") + theme(legend.position="top")


p3 =   
  ggplot(data = df_1651_4, aes(x = time1651, y = Freq, fill = NA1651))  + 
  geom_bar(stat = "identity", color = "black") + 
  theme_classic() + 
  guides(fill = guide_legend(title=NULL)) +
  xlab("")  + 
  scale_fill_brewer(palette="Blues") + 
  labs(title = "GsGD")


p4 = 
ggplot(data = df_1651_3, aes(x =sub1651_ng, y = Freq)) + 
  geom_bar(stat = "identity", color = "black", fill = "steelblue") + 
  theme_classic() + 
  labs(title = "non-GsGD") + xlab("")

multiplot(p1, p2, p4, ncol = 1)