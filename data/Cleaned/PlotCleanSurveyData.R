###################################################
# Some plots of the survey data for diagnostics ###
###################################################

rm(list = ls())

library(dplyr)
library(ggplot2)
library(maps)

##################################################

# Load in data
load(file.path('CelticSurveyFormattedSize.RData')) # Datras data by weight
load(file.path('CelticSurvey2FormattedSize.RData')) # Cefas data by weight

DWt <- DF; CWt <- FSS; rm(DF, FSS) # Rename to avoid confusion

load(file.path('DATRAS_No_At_Length.RData')) # Datras data by length 
load(file.path('Cefas_No_At_Length.RData')) # Cefas data by length 

Dln <- DF; Cln <- FSS; rm(DF, FSS)

##################################################

## Combine the datasets

Wt <- data.frame(Survey    = c(DWt$Survey      , as.character(CWt$fldSeriesName)),
		 Year      = c(DWt$Year        , CWt$Year),
		 Month     = c(DWt$Month       , CWt$Month),
		 HaulNo    = c(DWt$HaulNo      , CWt$fldCruiseStationNumber),
		 Lon       = c(DWt$HaulLonMid  , CWt$HaulLonMid),
		 Lat       = c(DWt$HaulLatMid  , CWt$HaulLatMid),
		 HaulDur   = c(DWt$HaulDur     , CWt$fldTowDuration),
		 SweptArea = c(DWt$SweptAreaAdj, CWt$SweptAreaAdj),
		 Species   = c(DWt$Species,      CWt$Species),
		 Kg        = c(DWt$Kg          , CWt$Kg))
rm(DWt, CWt)

Ln <- data.frame(Survey    = c(Dln$Survey      , as.character(Cln$fldSeriesName)),
		 Year      = c(Dln$Year        , Cln$Year),
		 Month     = c(Dln$Month       , Cln$Month),
		 HaulNo    = c(Dln$HaulNo      , Cln$fldCruiseStationNumber),
		 Lon       = c(Dln$HaulLonMid  , Cln$HaulLonMid),
		 Lat       = c(Dln$HaulLatMid  , Cln$HaulLatMid),
		 HaulDur   = c(Dln$HaulDur     , Cln$fldTowDuration),
		 SweptArea = c(Dln$SweptAreaAdj, Cln$SweptAreaAdj),
		 Species   = c(as.character(Dln$SpeciesName) , Cln$Species),
		 Length    = c(Dln$LngtClass   , Cln$Length),
		 Number    = c(Dln$Number      , Cln$Numbers),
		 Kg        = c(Dln$Kg          , Cln$Kg))

rm(Dln, Cln)
gc()

###################################################################


## Write to PDF
pdf(file = 'SurveyDataOverview.pdf', paper = 'a4r', bg = 'white', width = 0, height = 0)

##############################################
####### Plot the station locations ###########
##############################################

Stations <- Wt[!duplicated(paste(Wt$Survey, Wt$Year, Wt$Lon, Wt$Lat)),]

yrs <- sort(unique(Stations$Year))
n.yrs <- length(yrs)

map <- map_data('world', region = c('UK', 'Ireland', 'France'))


print(ggplot() + 
	geom_polygon(data = map, aes(x = long, y = lat, group = group), colour = 'black', fill = 'grey') +
	coord_fixed(xlim = c(-12, 2), ylim = c(48, 52), ratio = 1.3) + 
	geom_point(data  = Stations, aes(x = Lon, y = Lat, colour = Survey), shape = '+') + 
	facet_wrap(~ Year, ncol = 5) +
	theme_classic() + ggtitle('Survey locations by year and survey'))



#############################################
## Plot the temporal extent of the surveys ##
#############################################

surveyyrs <- reshape2::melt(table(Stations$Survey, Stations$Year))

print(ggplot(surveyyrs[surveyyrs$value !=0,], aes(x = Var2, y = Var1)) +
	geom_point(aes(size = value)) + xlab('') + ylab('') +
	theme(legend.title = element_blank()) + geom_vline(xintercept = 1997) + 
	ggtitle('Number of Stations Per Survey Per Year'))

############################################
############# Survey by month ##############
############################################

surveymo <- reshape2::melt(table(Stations$Month, Stations$Year))

print(ggplot(surveymo[surveymo$value !=0,], aes(x = Var2, y = Var1)) +
	geom_point(aes(size = value)) + xlab('') + ylab('') +
	theme(legend.title = element_blank()) + geom_vline(xintercept = 1997) + 
	ggtitle('Number of Stations Per Survey Per Month'))




#############################################
#####    Number of stations per year    #####
#############################################

surveyno <- group_by(Stations, Survey, Year) %>% summarise(n= n())

print(ggplot(surveyno, aes(x = Year, y = n)) + 
	geom_bar(stat = 'identity', aes(fill = Survey), colour = 'black') +
	theme_bw() + theme(axis.text.x = element_text(angle = -90)) +
	ggtitle('No stations per month'))

#############################################
#### Plot the survey spatial coverage    ####
#############################################

Lats_Lons <- group_by(Stations, Year) %>% summarise(minLon = min(Lon), maxLon = max(Lon), meanLon = mean(Lon),
						    minLat = min(Lat), maxLat = max(Lat), meanLat = mean(Lat))
print(ggplot(Lats_Lons, aes(x = Year, y = minLon)) + geom_segment(aes(xend = Year, yend = maxLon), lwd = 2) +
	geom_point(aes(y = meanLon), colour = 'red') +
	theme(axis.text.x = element_text(angle = -90)) + ylim(0, -14) + ylab('') + xlab('') +
	ggtitle('Longitudinal survey coverage: min, max and mean'))

print(ggplot(Lats_Lons, aes(x = Year, y = minLat)) + geom_segment(aes(xend = Year, yend = maxLat), lwd = 2) +
	geom_point(aes(y = meanLat), colour = 'red') +
	theme(axis.text.x = element_text(angle = -90)) + ylim(47, 53) + ylab('') + xlab('') +
	ggtitle('Latitudinal survey coverage: min, max and mean'))


#############################################
## Plot the catches for each species-group ##
#############################################


# put positive catches on same scale

spp <- unique(Wt$Species)

for (s in 1:length(spp)) {


plotDF <- Wt[Wt$Species == spp[s],]

print(ggplot() + 
	geom_polygon(data = map, aes(x = long, y = lat, group = group), colour = 'black', fill = 'grey') +
	coord_fixed(xlim = c(-12, 2), ylim = c(48, 52), ratio = 1.3) + 
	geom_point(data  = plotDF[plotDF$Kg != 0,], aes(x = Lon, y = Lat, size = sqrt(Kg)), colour = 'blue', alpha = 0.5) + 
	scale_size_continuous(limits = range(sqrt(Wt$Kg))) + 
	geom_point(data  = plotDF[plotDF$Kg == 0,], aes(x = Lon, y = Lat), colour = 'red', shape = '+') +
	facet_wrap(~ Year, ncol = 5) + 
	theme_classic() + ggtitle(spp[s]))


}

##############################################

##############################################
############ Plot CPUE by survey #############
##############################################

Wt$HaulDur <- as.numeric(as.character(Wt$HaulDur))

cpue <- group_by(Wt, Survey, Year, Species) %>% summarise(q05 = quantile(Kg/HaulDur * 60, prob = 0.05, na.rm = T),
							  q50 = quantile(Kg/HaulDur * 60, prob = 0.50, na.rm = T),
							  mean = mean(Kg/HaulDur * 60,na.rm = T),
							  q95 = quantile(Kg/HaulDur * 60, prob = 0.95, na.rm = T))

print(ggplot(cpue, aes(x = Year, y = mean)) + geom_line(aes(group = Survey, colour = Survey)) +
	facet_wrap(~Species, ncol = 1, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90)) +
	ylab('Kg per hour tow') + xlab(''))

cpsa <- group_by(Wt, Survey, Year, Species) %>% summarise(q05 = quantile(Kg/SweptArea, prob = 0.05, na.rm = T),
							  q50 = quantile(Kg/SweptArea, prob = 0.50, na.rm = T),
							  mean = mean(Kg/SweptArea,na.rm = T),
							  q95 = quantile(Kg/SweptArea, prob = 0.95, na.rm = T))

print(ggplot(cpue, aes(x = Year, y = mean)) + geom_line(aes(group = Survey, colour = Survey)) +
	facet_wrap(~Species, ncol = 1, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90)) +
	ylab('Density (catch per km2 swept)') + xlab(''))

#################################################
## Plot overall length comp by survey
#################################################

Lcomp <- group_by(Ln, Survey, Species, Length) %>% summarise(n = sum(Number)) %>% as.data.frame()

spp <- unique(Lcomp$Species)

for (s in 1:length(spp)) {

plotDF <- Lcomp[Lcomp$Species == spp[s],]

print(ggplot(plotDF, aes(x = Length, y = n)) + geom_line(aes(colour = Survey)) + facet_wrap( ~Species, scale = 'free_y'))

}

#################################################
## Plot length comp by year and survey
#################################################

Lcomp <- group_by(Ln, Survey, Year, Species, Length) %>% summarise(n = sum(Number)) %>% as.data.frame()

spp <- unique(Lcomp$Species)

for (s in 1:length(spp)) {

plotDF <- Lcomp[Lcomp$Species == spp[s],]

print(ggplot(plotDF, aes(x = Length, y = n)) + geom_line(aes(colour = Survey)) + facet_wrap(Year ~Species, scale = 'free_y'))

}
# close device 
dev.off()
