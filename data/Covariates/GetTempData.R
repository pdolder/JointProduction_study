#################################################################################
#########  Code to bulk download the Hadley interpolated SST records ############
## by Paul Dolder
## 07 / 04 / 2017
#################################################################################

rm(list=ls()) # clear your session

# file with multi-decadal aggregated data 
decades <- c('1870-1900', '1901-1930', '1931-1960', '1961-1990', '1991-2003')

# files with month/yearly data 
yrs     <- c(2004:2016)

# write a list of urls
urls_dec <- paste('http://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST1_SST_',decades,'.txt.gz', sep = '')
urls_yrs <- paste('http://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST1_SST_',yrs,'.txt.gz', sep = '')


## File format for years

# First line in info: month, year etc..
# Then 180 * 360 matrix of number
# -1000 is ICE, -32768 is land
# -36768 is not separated by a space, so we need to change

#Latitudes
minLat <- -89.5 
maxLat <- 89.5

# Longitudes
minLon <- -179.5 
maxLon <-  179.5 

# number of rows to skip for each month records
skips <- seq(1, 1981, l = 12) 


###############################################
### Function to download all the yearly data ##
###############################################

## Each year as a list

#Save to HadDat
HadDat <- lapply(urls_dec, function(url) {
# Download the yearly file and manipulate a bit
tmpFile <- tempfile()
fileName <- gsub(".gz","",basename(url))
download.file(url, tmpFile)
x <- readLines(tmpFile)
y <- gsub('-32768', ' 9999 ', x)
cat(y, file = tmpFile, sep = "\n")

## For each month, save the data in a list
monthLst <- lapply(as.list(skips), function(x) {
DF <- read.table(tmpFile, fill = T, sep = '', skip = x, nrow = 180)
DF <- apply(DF, 1:2, as.numeric) # make sure these are still numbers
colnames(DF) <- seq(minLon, maxLon, 1)
rownames(DF) <-seq(minLat, maxLat, 1)
return(DF) # return the monthly data
 })

# done with the year, unlink 
unlink(tmpFile)
print(praise::praise()) # bit of love

names(monthLst) <- paste('month', 1:12, sep ='_')
return(monthLst) # return the years data

# add a funny comment

})
names(HadDat) <- paste('year', yrs, sep = '_')


## test its ok

class(HadDat)
names(HadDat)

# function to rotate the data
rotate <- function(x) t(apply(x, 2, rev))

par(mfrow = c(3,4))
lapply(1:12, function(x) image(rotate(HadDat[['year_2004']][[x]]), zlim = c(-180, 3500), col = grey(seq(0.1,1,l = 100))))



