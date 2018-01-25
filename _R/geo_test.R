require(ggmap)
require(maptools)
require(rgdal)
require(plyr)

setwd("/Volumes/Macintosh_HD/Users/yaotli/Dropbox/labexchange/CHN_adm_shp/")

chn_adm1 <- readShapePoly("CHN_adm1.shp")
plot(chn_adm1)


plot(chn_adm1, col = c("red", rep("blue", 30)) )

cnN  <- c( 5, 11, 17, 18, 19, 20, 21, 22, 28 )
cnS  <- c( 4, 6, 9 ) 
cnSW <- c( 3, 7, 8, 26, 29, 30 ) 
cnC  <- c( 1, 12, 13, 14, 16, 25 ) 
cnE  <- c( 2, 10, 15, 23, 24, 27, 31) 

col_area_5       <- rep( "white", 31 )
col_area_5[cnN]  <- "#7f7f7f"
col_area_5[cnE]  <- "#2ca02c"
col_area_5[cnC]  <- "#8c564b"
col_area_5[cnSW] <- "#1f77b4"
col_area_5[cnS]  <- "#d62728"

plot(x = chn_adm1, col = col_area_5, border = "white") #jpg 600^2


colfunc <- colorRampPalette(c( "#1f77b4", "#d62728"))

p234         <- c(40, 145, 327, 146, 200)
p234         <- floor( p234/sum(p234)*100 )
col_area_all_234 <- as.raster(matrix(colfunc(100), ncol=1))[p234]

col_area_5_234 <- rep( "white", 31 )
col_area_5_234[cnN]  <- col_area_all_234[1]
col_area_5_234[cnE]  <- col_area_all_234[2]
col_area_5_234[cnC]  <- col_area_all_234[3]
col_area_5_234[cnSW] <- col_area_all_234[4]
col_area_5_234[cnS]  <- col_area_all_234[5]

plot(x = chn_adm1, col = col_area_5_234, border = "white") 



p232         <- c(17, 44, 55, 71, 40)
p232         <- floor(p232/sum(p232)*100)
col_area_all_232 <- as.raster(matrix(colfunc(100), ncol=1))[p232]

col_area_5_232 <- rep( "white", 31 )
col_area_5_232[cnN]  <- col_area_all_232[1]
col_area_5_232[cnE]  <- col_area_all_232[2]
col_area_5_232[cnC]  <- col_area_all_232[3]
col_area_5_232[cnSW] <- col_area_all_232[4]
col_area_5_232[cnS]  <- col_area_all_232[5]

plot(x = chn_adm1, col = col_area_5_232, border = "white") 


# plot(1:100, 1:100, pch = 19, cex=2, col = colfunc(100))
layout( matrix(1:3, ncol=3), width = c(2,2,1), height = c(1,1,1) )


legend_image <- as.raster(matrix(colfunc(100), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, 0, 1,1)






# https://cosx.org/2009/07/drawing-china-map-using-r/
# https://stackoverflow.com/questions/13355176/gradient-legend-in-base
# https://rstudio-pubs-static.s3.amazonaws.com/79029_b56eaffe36ef44f29b8efc0a07d67208.html


