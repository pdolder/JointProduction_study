####################################
## Plot of an area map and 
## the data and knots ##
#####################################


library(ggplot2); library(cowplot)
library(mapplots)
library(rworldmap)
library(rworldxtra)
library(data.table)
library(broom)
library(RColorBrewer)
library(rgdal)

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

## Read in ICES rectangles
rect <- readOGR(dsn = file.path('..', 'data','ICES_rect'), layer = 'ices_rectangles')
rect@data$id <- rownames(rect@data)
rectangles <- fortify(rect, region = 'id')
rectDF <- merge(rectangles, rect@data, by = "id")

# Read in ICES areas
area <- readOGR(dsn = file.path('..', 'data','ICES_area'), layer = 'ices_areas')
area@data$id <- rownames(area@data)
areas <- fortify(area, region = 'id')
areaDF <- merge(areas, area@data, by = "id")


## The plots

ggplot() + geom_rect(aes(xmin = xlim[1], xmax = xlim[2],  ymin =  ylim[1], ymax = ylim[2]), fill = 'lightblue', colour = 'grey') +
  geom_polygon(data = rectDF, aes(x= long, y = lat, group = group), fill = 'lightblue', colour = 'black',alpha = 0.3) +
  geom_polygon(data = areaDF, aes(x= long, y = lat, group = group), fill = 'lightblue', colour = 'black',alpha = 0.3, size = 2) +
  geom_text(aes(x = -4,  y = 49.5,  label = 'VIIe')) +
  geom_text(aes(x = -11, y = 50,    label = 'VIIj')) +
  geom_text(aes(x = -7,  y = 48.75, label = 'VIIh')) +
  geom_text(aes(x = -7,  y = 51.2,    label = 'VIIg')) +
  geom_text(aes(x = -5.5,  y = 50.5,  label = 'VIIf')) +
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), 
                                                                   legend.title = element_blank()) +
  ggtitle('Area of model application overlayed with ICES areas')
ggsave(file = file.path('..', 'plots', 'AreaMap.png'), width = 12, height = 8)


## Plot data and knots
## Plot map of area to which the model is appled

load(file.path('..', 'results', 'CovariatesAtKnot.RData'))
#load(file.path('..', 'results', '2017-04-08_M2', 'Save.RData'))


SpatialDeltaGLMM::Plot_data_and_knots(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=file.path('..', 'plots'))

# Create a dataframe of the knots
DF <- data.frame(X = Spatial_List$loc_i[,'E_km'], Y = Spatial_List$loc_i[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29

# Get the Lon Lats in decimal degrees
DataP <- PBSmapping::convUL(DF)

##
knots <- data.frame(X = Spatial_List$Kmeans$centers[,1], 
                    Y= Spatial_List$Kmeans$centers[,2])
attr(knots, 'projection') = 'UTM'
attr(knots, "zone") <- 29

knotsP <- PBSmapping::convUL(knots)


p1 <-  ggplot() + geom_rect(aes(xmin = xlim[1], xmax = xlim[2],  ymin =  ylim[1], ymax = ylim[2]), fill = 'lightblue', colour = 'grey') +
  geom_polygon(data = rectDF, aes(x= long, y = lat, group = group), fill = 'lightblue', colour = 'black',alpha = 0.3) +
  geom_polygon(data = areaDF, aes(x= long, y = lat, group = group), fill = 'lightblue', colour = 'black',alpha = 0.3, size = 2) +
  geom_text(aes(x = -4,  y = 49.5,  label = 'VIIe')) +
  geom_text(aes(x = -11, y = 50,    label = 'VIIj')) +
  geom_text(aes(x = -7,  y = 48.75, label = 'VIIh')) +
  geom_text(aes(x = -7,  y = 51.2,    label = 'VIIg')) +
  geom_text(aes(x = -5.5,  y = 50.5,  label = 'VIIf')) +
  coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), 
                                                                   legend.title = element_blank()) +
  geom_point(data = DataP, aes(x = X, y = Y), colour = 'blue', size = 0.1)

p2 <-  ggplot() + geom_rect(aes(xmin = xlim[1], xmax = xlim[2],  ymin =  ylim[1], ymax = ylim[2]), fill = 'lightblue', colour = 'grey') +
  geom_polygon(data = rectDF, aes(x= long, y = lat, group = group), fill = 'lightblue', colour = 'black',alpha = 0.3) +
  geom_polygon(data = areaDF, aes(x= long, y = lat, group = group), fill = 'lightblue', colour = 'black',alpha = 0.3, size = 2) +
  geom_text(aes(x = -4,  y = 49.5,  label = 'VIIe')) +
  geom_text(aes(x = -11, y = 50,    label = 'VIIj')) +
  geom_text(aes(x = -7,  y = 48.75, label = 'VIIh')) +
  geom_text(aes(x = -7,  y = 51.2,    label = 'VIIg')) +
  geom_text(aes(x = -5.5,  y = 50.5,  label = 'VIIf')) +
  coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), 
                                                                   legend.title = element_blank()) +
  geom_point(data = knotsP, aes(x = X, y = Y), colour = 'red', size = 0.8)

  plot_grid(p1,p2, labels = c('(a) data','(b) knots'), ncol = 1, vjust = c(1.5, -1))
  ggsave(file = file.path('..', 'plots', 'SpatialDataAndKnots.png'), width = 8, height = 12)
  
