# Functions applied for this project

# phylo_date ####
# convert YYYY-MM-DD to YY.date

phylo_date <- function(x){
  
  library(stringr)
  
  d = "([0-9]{4})-([0-9]{2})-([0-9]{2})"
  
  yr <- as.numeric(str_match(x, d)[,2])
  mo <- as.numeric(str_match(x, d)[,3])
  dy <- as.numeric(str_match(x, d)[,4])
  
  yr.0 <- paste0(yr, "-01-01")
  daydifference <- as.numeric(difftime(strptime(x, format = "%Y-%m-%d"),
                                       strptime(yr.0, format = "%Y-%m-%d"), 
                                       units = "days"))/365
  
  yr.daydifference = yr + daydifference
  yr.daydifference <- format( round(yr.daydifference, 2), nsmall = 2)
  
return(yr.daydifference)
  
}



