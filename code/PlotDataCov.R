
library(reshape2)
library(ggplot2)

load(file.path('..','results','2017-04-08_M2','Save.RData'))
load(file.path('..','results','CovariatesAtKnot.RData'))

DF <- Save$Data_Geostat
DF$Habitat <- Hab$Habitat[match(DF$knot_i, rownames(Hab))]

Depths <- data.frame(Depth = Depths, knot_i = 1:250)
DF$Depth   <- Depths$Depth[match(DF$knot_i, Depths$knot_i)]

ggplot(DF, aes(x = Habitat, y = log((Catch_KG / AreaSwept_km2)+1))) + geom_boxplot(notch = T) +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))

ggsave(file.path('..', 'plots', 'HabitatCovDiag.png'), width = 16, height = 16)

ggplot(DF, aes(x = abs(Depth), y = log((Catch_KG / AreaSwept_km2)+1))) + geom_point() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + geom_smooth()
ggsave(file.path('..', 'plots', 'DepthCovDiag.png'), width = 16, height = 16)


ggplot(DF, aes(x = log(abs(Depth)), y = log((Catch_KG / AreaSwept_km2)+1))) + geom_point() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))
