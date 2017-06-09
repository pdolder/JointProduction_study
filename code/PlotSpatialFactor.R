library(VAST)

load(file = file.path('..','results', '2017-04-08_M2', 'Save.RData'))

an <- as.numeric

DF <- Save$Data_Geostat

  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))

spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')


category_names = spp 
L_list  <- Var_list <-  vector('list', length = 4)
names(L_list)   <- names(Var_list) <-  c("Omega1", "Epsilon1", "Omega2", "Epsilon2")
Data = Save$TmbData
ParHat = Save$ParHat
Report = Save$Report
for(i in 1:4) {
Par_name = names(L_list)[i] 
if(i %in% c(1,3)) Var_name = paste('Omega','input',substring(Par_name, 6,6),'_sf', sep = '')
if(i %in% c(2,4)) Var_name = paste('Epsilon','input',substring(Par_name, 8,8),'_sft', sep = '')

L_list[[Par_name]] <- calc_cov(L_z = ParHat[[paste0('L_',tolower(Par_name), '_z')]], n_f = Data[['FieldConfig']][[Par_name]], n_c = Data$n_c, returntype = 'loadings_matrix')
rownames(L_list[[Par_name]]) <- category_names
Var_list[[Par_name]]  <- SpatialDFA::Rotate_Fn(L_pj = L_list[[Par_name]], Psi = Report[[Var_name]], RotationMethod = 'PCA', testcutoff = 1e-04)
rownames(Var_list[[Par_name]]$L_pj_rot) <- category_names
}



Mat_sf <- apply(Var_list$"Omega1"$Psi_rot, 1:2, FUN = mean)


## The map 

load(file.path('..', 'results', 'CovariatesAtKnot.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

plotDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers))

DF <- data.frame(X = Centers[,'E_km'], Y = Centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29 

LLs <- PBSmapping::convUL(DF)

plotDF$Lon <- LLs$Y
plotDF$Lat <- LLs$X

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
Mat_sf <- Mat_sf[1:250,]

Mat_sf <- as.data.frame(Mat_sf)
Mat_sf$x2i <- 1:250

DF2 <- MapDetails_List[['PlotDF']]

DF3 <- merge(DF2, Mat_sf)

DF3 <- DF3[DF3$Include==T,]

colnames(DF3)[5:13] <- paste("factor",1:9, sep = '_')

DF4 <- melt(DF3, id = c("x2i", "Lat","Lon","Include")

cols <- brewer.pal(8, 'Set1')

library(dplyr)

p1 <- ggplot() + geom_point(data = filter(DF4, variable %in% paste("factor", 1:3, sep="_")), aes(x = Lon, y =  Lat, colour = value), size = 1) +
			    geom_point(data = filter(DF4, variable %in% paste("factor", 1:3, sep="_")), aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +
  facet_wrap(~variable, ncol = 1)
  #, limits = lim)

p1

ggsave(file = file.path('..', 'plots', 'SpatialFactorLoadingsOmega1.png'))
