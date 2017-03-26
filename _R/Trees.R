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

  sub5496 <- read.tree(file.choose())
  T_sub5496 <-  ggtree(sub5496)

  # use findtaxa tip label

  serotype <- c("H5N8", "H5N6", "H5N2", "H5N3", "H5N5")
  serotype_col <- c("red", "red", "orange", "orange", "orange")

  sero_color <- findtaxa(type = 1, tree = sub5496, targetid = serotype, target = serotype_col)

  # label clade    

  fig_s <- ggtree(sub5496) %<+% sero_color + aes(color = I(colorr))

  fig_s1 <- fig_s + geom_cladelabel(identify(fig_s), label = "2.3.4.4", align = FALSE, fontsize = 5)

  fig_s8 = fig_s7 + geom_treescale()

  # export 

pdf(file="large.pdf", width = 13, height = 18.3)
print(fig_s8)
dev.off()

  # modify 

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

  fig_p <- ggtree(sub5496) %<+% sero_color2 + aes(color = I(colorr)) + geom_tippoint() + geom_treescale()
  fig_p1 = fig_p

  # label clade

  cladenode <- c("5523", "6930", "8014", "9056", "9521", "9955", "10606")
  cladenon <- c("2.3.4.4", "2.3.2", "2.2.1", "2.2.2", "2.1", "1", "0, 7")

  for (i in 1: 7){
  
  fig_p1 <- fig_p1 + geom_cladelabel(node = cladenode[i], label = cladenon[i], align = FALSE, fontsize = 5)
  
  }







