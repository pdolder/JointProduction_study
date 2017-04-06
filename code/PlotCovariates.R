
## Read in Spatial List

load(file.path('..', 'results', 'CovariatesAtKnot.RData'))

HabDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers),
		    Habitat = Hab$Habitat, Depth = Depths)


DF <- data.frame(X = Centers[,'E_km'], Y = Centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29 

LLs <- PBSmapping::convUL(DF)

HabDF$Lon <- LLs$Y
HabDF$Lat <- LLs$X

library(ggplot2); library(cowplot)
library(mapplots)
library(rworldmap)
library(rworldxtra)
library(data.table)
library(broom)

mapdata =ggplot2::fortify(rworldmap::getMap(resolution = 'high'))

xlim = c(-12, -2)
ylim = c(48, 52)
 # if you want to fast subset by lat/long you can do the following:
 # subset to just regions in xlim and ylim, see
# http://stackoverflow.com/a/16574176
 mapdata = data.table(mapdata)
  mapdata = mapdata[mapdata[,.I[any(long %between% xlim) & any(lat %between%
 ylim)], by = list(group)]$V1]

  # make geoms
  coast.poly = geom_polygon(data=mapdata, aes(x=long, y=lat, group=group), colour="#999999", fill="#999999", lwd=0.2)
  coast.outline = geom_path(data=mapdata, aes(x=long, y=lat, group=group), colour="#000000", lwd=0.2)

p1 <- ggplot() + coast.poly + coast.outline + geom_point(data = HabDF, aes(x = lat, y =  lon, colour = Habitat)) +   coord_quickmap(xlim, ylim) +
	theme(legend.position = 'top')
p2 <- ggplot() + coast.poly + coast.outline + geom_point(data = HabDF, aes(x = lat, y = lon, colour = Depth)) + coord_quickmap(xlim, ylim) + 
	theme(legend.position = 'bottom')

plot_grid(p1,p2, ncol = 1)

ggsave(file.path('..', 'plots', 'HabitatCovariates.png'), width = 8, height = 12)
