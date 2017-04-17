
## Read in Spatial List

load(file.path('..', 'results', 'CovariatesAtKnot.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

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
library(RColorBrewer)

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


#####################################
# Plotting all areas

HabDF2 <- MapDetails_List[['PlotDF']]

HabDF2$Habitat <- HabDF$Habitat[match(HabDF2$x2i, HabDF$knot)]
HabDF2$Depth   <- HabDF$Depth[match(HabDF2$x2i, HabDF$knot)]

p1 <- ggplot() + geom_point(data = HabDF2, aes(x = Lon, y =  Lat, colour = Habitat), size = 0.02) + geom_point(data = HabDF, aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
	coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'top', legend.title = element_blank()) + 
	guides(colour = guide_legend(override.aes = list(size = 4, shape = 15)))

cols <- brewer.pal(9, 'Blues')[c(9,5,1)]

p2 <- ggplot() + geom_point(data = HabDF2, aes(x = Lon, y = Lat, colour = Depth), size = 0.02) + geom_point(data = HabDF, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
	coast.poly + coast.outline + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), legend.title = element_blank()) + 
	scale_colour_gradient2(low = cols[1], mid = cols[2], high = cols[3], midpoint = mean(HabDF$Depth), space = "Lab", na.value = "grey50", guide = "colourbar")

plot_grid(p1,p2, labels = c('(a) Substrate', '(b) Depth'), vjust = c(1.5, -1), ncol = 1)

ggsave(file.path('..', 'plots', 'HabitatCovariatesPoly.png'), width = 8, height = 12)


