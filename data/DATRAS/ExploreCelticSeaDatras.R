### Datras download code

library(icesDatras)

getSurveyList()  # returns the acronyms for available surveys.

getSurveyYearList("EVHOE") # returns the years available for a given survey.

getSurveyYearQuarterList(survey = "EVHOE", year = 1998)  # returns the quarters available for a given survey and year.

# Single year
HH <- getHHdata(survey = "EVHOE", year = 1997, quarter = 4)

# Multi-year
HH <- do.call(rbind,
                  lapply(1997:2015,
                         function(year)
                         getHHdata(survey = "EVHOE", year = year, quarter = 4)
                             )
                     )

getSurveyYearList("IE-IGFS") # returns the years available for a given survey.

getSurveyYearQuarterList(survey = "IE-IGFS", year = 2003)  # returns the quarters available for a given survey and year.


HH_IE <- do.call(rbind,
                  lapply(2003:2015,
                         function(year)
                         getHHdata(survey = "IE-IGFS", year = year, quarter = 4)
                             )
                     )



# Add the IGFS
HH <- rbind(HH,HH_IE)

# Plot the shoot -hauls

library(ggmap)
library(ggplot2)

HH$ShootLat  <- as.numeric(HH$ShootLat) 
HH$ShootLong <- as.numeric(HH$ShootLong) 
HH$HaulLat  <- as.numeric(HH$HaulLat) 
HH$HaulLong  <- as.numeric(HH$HaulLong) 


map <- get_map("Celtic Sea", zoom = 6)

png(file = file.path("CelticSurveyStations.png"), width = 800, height = 800)
print(ggmap(map) + geom_segment(data = HH,aes(x = ShootLong, xend = HaulLong, y = ShootLat, yend = HaulLat, colour = Survey)) +
	scale_colour_manual(values = c("red","darkgreen")))
dev.off()



pdf(file.path("CelticSeaSurveys.pdf"), paper = "a4r")

# By Survey
print(ggmap(map) + geom_segment(data = HH,aes(x = ShootLong, xend = HaulLong, y = ShootLat, yend = HaulLat, colour = Survey)) +
	scale_colour_manual(values = c("red","darkgreen")))

# By year
print(ggmap(map) + geom_segment(data = HH,aes(x = ShootLong, xend = HaulLong, y = ShootLat, yend = HaulLat, colour = Year)) +
facet_wrap(~Survey))

dev.off()


## And the length data

# getHLdata() single year

# muti-year
HL <- do.call(rbind,
                  lapply(1997:2015,
                         function(year)
                         getHLdata(survey = "EVHOE", year = year, quarter = 4)
                             )
                     )
# muti-year
HL_IE <- do.call(rbind,
                  lapply(2003:2015,
                         function(year)
                         getHLdata(survey = "IE-IGFS", year = year, quarter = 4)
                             )
                     )

HL <- rbind(HL, HL_IE)

save(HH, HL, file = file.path("CelticSurveyData.RData"))
