
Index <- read.csv('../../results/2017-02-04_RUN_FULL/Table_for_SS3.csv')

## ICES assessment SSB
CS <-	read.csv(file.path('..','..','..','SFG_talk','CSAssessmentOutput2016.csv'))

library(dplyr) 
CS <- filter(CS, StockDescription %in%  unique(CS$StockDescription)[c(5,7,12)])

CS$Category <-paste(CS$SpeciesName,'Adu', sep = '_')

# scale the data
CS_scaled <- by(CS[,c('High_StockSize','StockSize','Low_StockSize')],
		CS[,'Category'], scale)

CS$scaled <- c(CS_scaled[['Melanogrammus aeglefinus_Adu']][,'StockSize'],
	       CS_scaled[['Merlangius merlangus_Adu']][,'StockSize'],
	       CS_scaled[['Gadus morhua_Adu']][,'StockSize'])

CS$scaledMult <- (CS$High_StockSize - CS$StockSize) / CS$StockSize

Index <- filter(Index, Category %in% unique(grep('Adu', Index$Category, value = T))) 

Index_scaled <- by(Index[,c('Estimate_metric_tons')], Index[,'Category'], scale)

Index$scaled <- c(Index_scaled[['Gadus morhua_Adu']][,1],
		  Index_scaled[['Melanogrammus aeglefinus_Adu']][,1],
				Index_scaled[['Merlangius merlangus_Adu']][,1])

#Index$scaledMult <- (Index$Hi -
#		     Index$Estimate..metric.tonnes.)/Index$Estimate..metric.tonnes.


## Plot indices

Index <- filter(Index, Category %in% c('Gadus morhua_Adu', 'Gadus morhua_Juv',
				      'Melanogrammus aeglefinus_Adu', 'Melanogrammus aeglefinus_Juv',
				      'Merlangius merlangus_Adu','Merlangius merlangus_Juv'),
		Year %in% 1990:2015)

CS <- filter(CS, Category %in% c('Gadus morhua_Adu', 'Gadus morhua_Juv',
				      'Melanogrammus aeglefinus_Adu', 'Melanogrammus aeglefinus_Juv',
				      'Merlangius merlangus_Adu','Merlangius merlangus_Juv'),
	     Year %in% 1990:2015)




library(ggplot2)
print(ggplot(Index, aes(x = Year, y = scaled)) + geom_line() +
#      geom_pointrange(aes(ymin = scaled - scaled * scaledMult, ymax = scaled +  scaled * scaledMult)) +
      facet_wrap(~Category, scale =    'free_y', ncol = 1) + 	expand_limits(y = 0) +
      geom_line(data = CS, aes(x = Year, y = scaled), col = 'blue') +
#      geom_ribbon(data = CS, aes(x =  Year,  ymin =   scaled   - scaled *   scaledMult, ymax =
#    scaled  + scaled *   scaledMult),alpha = 0.5, fill ='blue') +
      ggtitle('Blue = Assessment, black = VAST estimates') + theme_bw())
ggsave('RealativeIndexVRelativeAssessSSB.png')


