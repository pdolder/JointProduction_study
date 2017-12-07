##############################################
# Code plotting the standardised 
# spatial densities of all species 
# In a single year
#############################################

rm(list = ls())

run <- '2017-06-16_M1'

load(file.path('..','..', '..', '..','results', 'CovariatesAtKnot.RData'))
load(file.path('..','..', '..' , '..','results', run, 'Save.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

den <- list()

# standardise the densities for a year (2015)
den <- lapply(1:18, function(x) {
		      dens <- Save$Report$Index_xcyl[,x,,]
		      dens <- dens[,26] / sum(dens[,26]) * 100
		      return(dens)
})

# species names
spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')
names(den) <- spp

den <- as.data.frame.table(do.call(rbind, den))
den$Var2 <- rep(1:250, each = 18)
colnames(den) <- c("spp", "x2i", "value")

#############
## the Map ##
#############

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


## The map for all areas

DF2 <- MapDetails_List[['PlotDF']]
DF2 <- DF2[DF2$Include == T,]

PlotDF <- merge(DF2, den)

### The Plot

cols <- brewer.pal(8, 'Set1')

## The plots

for (i in 1:18) {

assign(paste("p",i, sep = ""), ggplot() + geom_point(data = PlotDF[PlotDF$spp==spp[i],], aes(x = Lon, y =  Lat, colour = value), size = 1) + geom_point(data = PlotDF[PlotDF$spp==spp[i],], aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") + facet_wrap(~spp)) 

}

plot_grid(p1,p3,p5,p7,p9,p11,p13,p15,p17, ncol = 3)

plot_grid(p2,p4,p6,p8,p10,p12,p14,p16,p18, ncol = 3)
