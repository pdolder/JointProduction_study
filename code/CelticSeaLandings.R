
library(dplyr)

DF <- read.csv(file.path('..', 'data', 'STECF','landings_and_discards_dataALL.csv'))

DF$InModel <- ifelse(DF$species %in% c("COD","HAD","WHG",
				       "HKE","LEZ","MEG","ANF",
				       "PLE","SOL"), "Yes", "No")

DF2 <- DF %>% filter(Measure.Names == 'landings', year %in% 2011:2015) %>%
	group_by(InModel) %>% summarise(value = sum(Measure.Values, na.rm = T)) %>%
						 as.data.frame()

print(DF2)



