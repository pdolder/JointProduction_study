#################################################################################
#########  Code to bulk download the Hadley interpolated SST records ############
## by Paul Dolder
## 07 / 04 / 2017
#################################################################################

rm(list=ls()) # clear your session

#################################################
dat.type <- 'dec.dat'  ## Option for data type ##
#################################################

# file with multi-decadal aggregated data 
decades <- c('1870-1900', '1901-1930', '1931-1960', '1961-1990', '1991-2003')[1]
# files with month/yearly data 
yrs     <- c(2004:2016)

# write a list of urls
if(dat.type == 'dec.dat')  { urls <- paste('http://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST1_SST_',decades,'.txt.gz', sep = '') }
if(dat.type == 'yrs.dat')  { urls <- paste('http://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST1_SST_',yrs,'.txt.gz', sep = '')}

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


###############################################
### Function to download all the yearly data ##
###############################################

if(dat.type == 'yrs.dat') {
# number of rows to skip for each month records
nrows   <- 180
nmonths <- 12
start<-1
skips <- seq(start, nrows*nmonths-nrows+start, l = nmonths) 


## Each year as a list

#Save to HadDat
HadDat <- lapply(urls, function(url) {
# Download the yearly file and manipulate a bit
tmpFile <- tempfile()
fileName <- gsub(".gz","",basename(url))
download.file(url, tmpFile)
x <- readLines(tmpFile)
y <- gsub('-32768', ' -9999 ', x)
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

###############################################

## Save the output as 'flat' data
flat <- lapply(HadDat, function(yr) {
	 res <- lapply(yr, function(mt) {
		DF <- as.data.frame(mt)
		DF$Lat <- rownames(DF)
		DF <- reshape2::melt(DF, id = 'Lat')
		colnames(DF)[2] <- 'Lon'
		DF$Lat <- -as.numeric(as.character(DF$Lat))
		DF$Lon <-  as.numeric(as.character(DF$Lon))
		return(DF)
 })
 names(res) <- 1:12
 res <-  do.call(rbind,res)
 res$month <- sapply(strsplit(rownames(res), '\\.'),'[[',1)
 return(res)

})

flat <- do.call(rbind, flat)
flat$year <- sapply(strsplit(rownames(flat), '\\.'), '[[',1)
flat$year <- sapply(strsplit(flat$year, '\\_'), '[[',2)

library(data.table)

flat <- data.table(flat)
save(flat, file = 'YearlyTempData.RData')

}



###############################################
## Function to download all the decadal data ##
###############################################

if(dat.type == 'dec.dat') {

#urls <-  urls[5] # 1991 - 2003 
urls <-  urls[1] # 1961 - 1990 

# number of rows to skip for each month records
nrows   <- 180
nmonths <- 12
#nyears <-  length(1991:2003) 
nyears <-  length(1961:1990) 
start <-1
skips <- seq(start, nrows * nmonths * nyears - nrows + start, l = nmonths * nyears) 

## Each year as a list

#Save to HadDat
HadDat <- lapply(urls, function(url) {
# Download the yearly file and manipulate a bit
tmpFile <- tempfile()
fileName <- gsub(".gz","",basename(url))
download.file(url, tmpFile)
x <- readLines(tmpFile)
y <- gsub('-32768', ' -9999 ', x)
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

return(monthLst) # return the years data

# add a funny comment

})

yrs <- substr(decades,1,4):c(as.numeric(substr(decades,6,9))-1)
mnths <- 1:12
names(HadDat[[1]]) <- paste(rep(yrs, each = 12), mnths, sep = '_')

###############################################

## Save the output as 'flat' data
flat <- lapply(HadDat[[1]], function(x) {
		DF <- as.data.frame(x)
		DF$Lat <- rownames(DF)
		DF <- reshape2::melt(DF, id = 'Lat')
		colnames(DF)[2] <- 'Lon'
		DF$Lat <- -as.numeric(as.character(DF$Lat))
		DF$Lon <-  as.numeric(as.character(DF$Lon))
		return(DF)

})

flat <- do.call(rbind, flat)

flat$year  <- sapply(strsplit(rownames(flat), '\\_'), '[[',1)
flat$month <- sapply(strsplit(rownames(flat), '\\_'), '[[',2)
flat$month <- sapply(strsplit(flat$month, '\\.'),'[[',1)

library(data.table)

flat <- data.table(flat)
save(flat, file = 'DecTempData1.RData')

}

# Trim the files to something more useable for the Celtic Sea
# -12,  2 Lon
#  48, 52 Lat
library(dplyr)
flat2 <- filter(flat, Lat > -13, Lat < 3, Lon > 47, Lon < 53)

#Dec1 <- flat
#Dec2 <- flat
#Year <- flat

CS <- rbind(Dec1, Dec2, Year)

save(CS, file = 'CelticSeaTempRecords_1961_2016.RData')

## plot a single point

plotDF <- filter(CS, Lat == -7.5, Lon == 51.5, month == 3, year > 2000)

plotDF2 <- plotDF[!plotDF$value %in% c(-9999, -1000),]

plot(plotDF2$year, plotDF2$value/1000)
