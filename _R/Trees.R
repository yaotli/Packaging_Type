# Aim of the project: 
# to understand the subtype dist. in Gs/GD and 
# inspect the essential cis-determinants
# for this section:
#
# 1. Tree editing
# 
# Data sources:
# 1. gisaid: 6278/ H5
# 2. ncbi: 6677/ H5
# 3. reference strains with clade label from Dr. Gavin Smith (238)

# Tree manipulation ####

  # p <- p %>% rotate(identify(p))
  # fig_s + geom_text(aes(label=node), size = 0.5)

library(ape)
library(ggtree)

  # read fron .nwk filw

  sub5496 <- read.tree("~/Desktop/TAC1/Sub5496_note")
  T_sub5496 <-  ggtree(sub5496)

  
  # time-isolate distribution 
  
  Year_id <- c("2014", "2015", "2016", "2017")
  yeardis <- findtaxa(type = 0, tree = sub5496, Year_id, rep("orange", 4))
  
  T_sub5496_yr <- T_sub5496 %<+% yeardis + 
    geom_tippoint(aes(color = I(shapee)), size = 1, alpha = 0.8)
  
  # clade 
  
  cladenode <- c("5523", "6930", "8014", "9056", "9521", "9955", "10606")
  cladenon <- c("2.3.4.4", "2.3.2", "2.2.1", "2.2.2", "2.1", "1", "0, 7")
  
  for (i in 1: 7){
    
    T_sub5496_yr <- T_sub5496_yr + 
      geom_cladelabel(node = cladenode[i], label = cladenon[i], align = FALSE, fontsize = 5)
    
  }
  
  
  
  # use findtaxa tip label

  serotype <- c("H5N8", "H5N6", "H5N2", "H5N3", "H5N5")
  serotype_col <- c("red", "red", "orange", "orange", "orange")

  sero_color <- findtaxa(type = 1, tree = sub5496, targetid = serotype, target = serotype_col)

  # label clade    

  fig_s <- ggtree(sub5496) %<+% sero_color + aes(color = I(colorr))

  fig_s1 <- fig_s + geom_cladelabel(identify(fig_s), label = "2.3.4.4", align = FALSE, fontsize = 5)

  # export 

pdf(file="large.pdf", width = 13, height = 18.3)
print(fig_s8)
dev.off()

# modify using ggplot2 color 

  gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n = 3
  cols = gg_color_hue(n)
  
  # plot(1:n, pch = 16, cex = 2, col = cols)

  # tip label (ggplot sytle color)
  
  serotype <- c("H5N8", "H5N6", "H5N2", "H5N3", "H5N5")
  serotype_col2 <- c("#F8766D", "#F8766D", "#00BA38", "#00BA38", "#00BA38")

  sero_color2 <- findtaxa(type = 1, tree = sub5496, targetid = serotype, target = serotype_col2)

  fig_p <- ggtree(sub5496) %<+% 
    sero_color2 + aes(color = I(colorr)) + 
    geom_tippoint() + geom_treescale()
  
  fig_p1 = fig_p

  # label clade

  cladenode <- c("5523", "6930", "8014", "9056", "9521", "9955", "10606")
  cladenon <- c("2.3.4.4", "2.3.2", "2.2.1", "2.2.2", "2.1", "1", "0, 7")

  for (i in 1: 7){
  
  fig_p1 <- fig_p1 + 
    geom_cladelabel(node = cladenode[i], label = cladenon[i], align = FALSE, fontsize = 5)
  
  }

# Small tree

  sub5496_2344 <- groupClade(sub5496, node = 5523)
  
  #  sub5496_2344_p <- ggtree(sub5496_2344, size = 0.2, aes(color=group)) + 
  #  scale_color_manual(values =  c("black", "#F8766D")) + geom_treescale()
  
  n = 5
  cols_subtype = c("black", gg_color_hue(n) )
  serotype = c("H5N1", "H5N2", "H5N3", "H5N5", "H5N6", "H5N8")
  
  sero_colo3 <- findtaxa(type = 1, sub5496_2344, targetid = serotype, target = cols_subtype)
  
  
  sub5496_2344_p <- ggtree(sub5496_2344, size = 0.3) %<+% 
    sero_colo3 + aes(color = I(colorr)) 
  
  sub5496_2344_p0 <- sub5496_2344_p + 
    geom_hilight(node = 5523, fill = "#fff7bc", alpha = 0.3) 
  
  sub5496_2344_p1 = sub5496_2344_p0
  
  for (i in 1: 7){
    
    sub5496_2344_p1 <- sub5496_2344_p1 + 
      geom_cladelabel(node = cladenode[i], label = cladenon[i], align = FALSE, fontsize = 12)
    
  }


# 2.3.4.4 clade tree
  
  sub888 <- read.tree(file.choose())
  T_sub888 <-  ggtree(sub888) + geom_treescale()
  
  sero_color_2344 <- findtaxa(type = 1, sub888, targetid = serotype, target = cols_subtype)
  
  # try to generate legend of subtype
  
  sero_color_2344[,12] = "N1"
  
  for (k in 1: length(sero_color_2344[,12])){
    
    sero_color_2344[,12][k] = 
    sub(pattern = "H5", replacement = "", serotype)[
      match(sero_color_2344[,11][k], cols_subtype)]
      
  }
  colnames(sero_color_2344)[12] = "Ntype"
  
  T_sub888_fornote = T_sub888 %<+% sero_color_2344 + 
    geom_tippoint(aes(color = Ntype), size = 0.8) + 
    scale_color_manual(values=cols_subtype) + 
    theme(legend.position = "left", legend.text = element_text(size=30)) +
    guides(colour=guide_legend("Ntype", override.aes = list(size =20)))
  
  
  T_sub888_p0 = T_sub888 %<+% sero_color_2344 + aes(color = I(colorr)) + 
    geom_tippoint(size = 0.8)

  
  
  

  
  
  
