
library(reshape2)
library(ggplot2)

load(file.path('..','results','2017-04-08_M2','Save.RData'))
load(file.path('..','results','CovariatesAtKnot.RData'))

DF <- Save$Data_Geostat
DF$Habitat <- Hab$Habitat[match(DF$knot_i, rownames(Hab))]

Depths <- data.frame(Depth = Depths, knot_i = 1:250)
DF$Depth   <- Depths$Depth[match(DF$knot_i, Depths$knot_i)]

ggplot(DF, aes(x = Habitat, y = Catch_KG / AreaSwept_km2)) + geom_boxplot() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))

ggplot(DF, aes(x = abs(Depth), y = Catch_KG / AreaSwept_km2)) + geom_point() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))

ggplot(DF, aes(x = log(abs(Depth)), y = Catch_KG / AreaSwept_km2)) + geom_point() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))
