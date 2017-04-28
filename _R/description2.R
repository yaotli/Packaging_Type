library(seqinr)
library(stringr)
source("./_R/Function.R")

# souce: 1. subSeq_GsGD_2085
#        2. subSeq_nonGsGD_644
#        3. subSeq_2344_498


  for (x in 1:3){
    
    # input
    
    filedir = c("~/Desktop/subSeq_GsGD_2085.fasta", 
                "~/Desktop/subSeq_nonGsGD_644.fasta", 
                "~/Desktop/subSeq_2344_498.fasta")
    
    file = read.fasta(filedir[x])
    seq_name0 = attributes(file)$names
         seq0 = getSequence(file)
 
         
    # year info (Year_x)
         
    if ( length( which(endsWith(seq_name0, "NA") == TRUE) ) > 0 ){
      
      seq_name = seq_name0[-which(endsWith(seq_name0, "NA") == TRUE)]
      seq  = seq0[-which(endsWith(seq_name0, "NA") == TRUE)]
      
    }else{
      
      seq_name = seq_name0
      seq = seq0
      
    }
    
    year0 <- str_match(seq_name, "_([0-9]{4})\\.([0-9]+)")[,2]
    
    Year = c()     
    
    for (k in 1: length(seq)){
      
      if (year0[k] < "1996"){
        
        Year[length(Year) + 1] = "1995"
        
      }else{
        
        Year[length(Year) + 1] = year0[k]
      }
    }
    
    
    assign(paste0("Year_", x), Year)
    
    
    # NA (Subtype_x)
    
    Subtype <- gsub(pattern = "H5", 
                    replacement = "", 
                    str_match(pattern = "_(H[0-9]+N[0-9]+)_", 
                              seq_name)[,2])
    
    assign(paste0("Subtype_", x), Subtype)
    
    
    # nucleotide info
    
    seq_matrix = do.call(rbind, seq )
    
    p24 = seq_matrix[,24]
    p35 = seq_matrix[,35]
    
    
    # make into dataframe

    groupnames <- c("GsGD", "nonGsGD", "Clade2344")
    
# df_subtype_y
    
    df = as.data.frame(
      table( Subtype, Year ), 
      stringsAsFactors = FALSE
    )
    df[,4] = groupnames[x]
    
    df2 = as.data.frame(
      prop.table ( table( Subtype, Year ), margin = 2), 
      stringsAsFactors = FALSE
    )
    df2[,4] = groupnames[x]
    
    colnames(df)[4] = "group"
    colnames(df2)[4] = "group"
    
    assign( paste0("df_subtype_y_", x), df)
    assign( paste0("df_Psubtype_y_", x), df2)
    
    
# df_p24_y    
    
    df = as.data.frame(
      prop.table( table( p24, Year, exclude = "-"), margin = 2),
      stringsAsFactors = FALSE
    )
    
             df[,4] = groupnames[x]
    colnames(df)[4] = "group"
    
    assign( paste0("df_p24_y_", x), df)
    
    
# df_p35_y
    
    df = as.data.frame(
      prop.table( table( p35, Year, exclude = "-"), margin = 2),
      stringsAsFactors = FALSE
    )
    
             df[,4] = groupnames[x]
    colnames(df)[4] = "group"
    
    assign( paste0("df_p35_y", x), df )
    
    
# df_p24_subtype
    
    df = as.data.frame(
      prop.table( table( p24, Subtype, exclude = "-"), margin = 2),
      stringsAsFactors = FALSE
    )
    
    df[,4] = groupnames[x]
    colnames(df)[4] = "group"
    
    assign( paste0("df_p24_subtype", x), df )
    

# df_p35_subtype
    
    df = as.data.frame(
      prop.table( table( p35, Subtype, exclude = "-"), margin = 2),
      stringsAsFactors = FALSE
    )
    
    df[,4] = groupnames[x]
    colnames(df)[4] = "group"
    
    assign( paste0("df_p35_subtype", x), df )
    
      }    

       
# combine df


   df_subtype_y <- rbind(df_subtype_y_1, df_subtype_y_2, df_subtype_y_3)
  df_Psubtype_y <- rbind(df_Psubtype_y_1, df_Psubtype_y_2, df_Psubtype_y_3)
         
      df_p24_y <- rbind(df_p24_y_1, df_p24_y_2, df_p24_y_3)     
  
      df_p35_y <- rbind(df_p35_y1, df_p35_y2, df_p35_y3)     
       
      df_p24_subtype <- rbind(df_p24_subtype1, df_p24_subtype2, df_p24_subtype3)
      df_p35_subtype <- rbind(df_p35_subtype1, df_p35_subtype2, df_p35_subtype3)
  
# ggplot
  
  
# year_subtype
  
library(ggplot2)  
  
  p = ggplot(data = df_subtype_y, aes(x = Year, y = Freq, fill = Subtype))  + 
    geom_bar(stat = "identity") +
    theme_classic() + 
    guides(fill = guide_legend(title=NULL)) +
    xlab("")  + 
    scale_fill_manual(values = c("gray" ,gg_color_hue(8)) ) +
    facet_wrap(~group, ncol = 1) + 
    theme(legend.position = "top")
  
  
 
  # prop. of subtype
  
  f1 = ggplot(data = df_Psubtype_y, aes(x = Year, y = Freq, group = Subtype, colour = Subtype))  + 
    geom_line(size = 1.5) + 
    guides(colour = guide_legend(title = "")) +
    facet_wrap(~ group, ncol = 1) +
    theme_bw() + 
    theme(strip.background = 
            element_rect(colour = "black", fill = "white", size = 1.5)) +
    xlab("") + 
    scale_color_manual(values = c("gray" ,gg_color_hue(8)) ) + 
    theme(legend.position = "top")
  
  
  
# nt_year
  
  # p24
  f2 = ggplot(data = df_p24_y, aes(x = Year, y = Freq, group = p24, colour = p24))  + 
    geom_line(size = 1.5) + 
    facet_wrap(~ group, ncol = 1) +
    theme_bw() + 
    theme(strip.background = 
            element_rect(colour = "black", fill = "white", size = 1.5)) +
    xlab("") + 
    theme(legend.position = "top")
  
  f3 = ggplot(data = df_p35_y, aes(x = Year, y = Freq, group = p35, colour = p35))  + 
    geom_line(size = 1.5) + 
    facet_wrap(~ group, ncol = 1) +
    theme_bw() + 
    theme(strip.background = 
            element_rect(colour = "black", fill = "white", size = 1.5)) +
    xlab("") + 
    theme(legend.position = "top")
    
  
library(ggtree)
  
  multiplot(p, f2, f3, ncol = 1)  
  
  

# nt_subtype
  
  # p24
  
  f5 = ggplot(data = df_p24_subtype, aes(x = Subtype, y = Freq, fill = p24))  + 
    geom_bar(stat = "identity") +
    theme_classic() + 
    guides(fill = guide_legend(title="p24")) +
    xlab("")  + 
    facet_wrap(~group, ncol = 1) + 
    theme(legend.position = "top")
  
  f6 = ggplot(data = df_p35_subtype, aes(x = Subtype, y = Freq, fill = p35))  + 
    geom_bar(stat = "identity") +
    theme_classic() + 
    guides(fill = guide_legend(title="p35")) +
    xlab("")  + 
    facet_wrap(~group, ncol = 1) + 
    theme(legend.position = "top")
  
  multiplot(f5, f6, ncol = 2)  
  
  
  
            
       
       
       
  
  