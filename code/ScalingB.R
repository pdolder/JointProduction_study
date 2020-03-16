## Calculate the scaling factors for the biomass estimates

library(dplyr)

DF <- read.csv(file.path('..', 'data', 'STECF', 'landings_and_discards_data.csv'))

## Piscatorius and budegassa taken from ICES advice sheets

DF <- rbind(DF, cbind(Measure.Names = rep('landings',26), reg_area_cod = rep('7BCEFGHJK',26), 
		      species = rep(c('PISC', 'BUD'), each = 13), year  = rep(2003:2015, 2),
		      Measure.Values = c(21665,23741,22098,22490,25432,21248,16197,16344,16358,18316,18324,21205,21250,
					 6899, 5767, 5810, 4305, 4690, 5475, 6537, 6994, 6100, 6055, 7669, 6744, 6669)))

DF$Measure.Values <- as.numeric(DF$Measure.Values)

DF <- DF %>% filter(Measure.Names == 'landings', !species %in% c('ANF','MON','MEG')) %>% group_by(species) %>%
	summarise(MaxLandings = max(Measure.Values))

MaxLan <- DF

save(MaxLan, file = 'MaxLandingsPerSpecies.RData')


