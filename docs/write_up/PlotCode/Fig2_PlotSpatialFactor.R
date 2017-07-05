library(VAST)

run <- '2017-06-16_M1'
load(file = file.path('..', '..', '..','results', run, 'Save.RData'))

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


Omega1_sf   <- apply(Var_list$"Omega1"$Psi_rot, 1:2, FUN = mean)
Omega2_sf   <- apply(Var_list$"Omega2"$Psi_rot, 1:2, FUN = mean)
Epsilon1_sf <- Var_list$"Epsilon1"$Psi_rot
Epsilon2_sf <- Var_list$"Epsilon2"$Psi_rot

##############
## The map  ##
##############
load(file.path('..', '..', '..','results', 'CovariatesAtKnot.RData'))

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
library(dplyr) 

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
Omega1_sf <- Omega1_sf[1:250,]
Omega2_sf <- Omega2_sf[1:250,]
Epsilon1_sf <- Epsilon1_sf[1:250,1:9,1:26]
Epsilon2_sf <- Epsilon2_sf[1:250,1:9,1:26]

Omega1_sf <- as.data.frame(Omega1_sf)
Omega2_sf <- as.data.frame(Omega2_sf)
Epsilon1_sf <- as.data.frame.table(Epsilon1_sf)
Epsilon2_sf <- as.data.frame.table(Epsilon2_sf)

Epsilon1_sf$Var1 <- Epsilon2_sf$Var1 <- 1:250
Epsilon1_sf$Var2 <- Epsilon2_sf$Var2 <- rep(paste("factor", 1:9, sep = "_"), each = 250)
Epsilon1_sf$Var3 <- Epsilon2_sf$Var3 <- rep(1990:2015, each = 250 * 9)

colnames(Epsilon1_sf) <- colnames(Epsilon2_sf) <- c('x2i','factor','year','value')

## Just the first 3 factors for the Epsilons and every 5 years, otherwise files are huge
Epsilon1_sf <- filter(Epsilon1_sf, factor %in% paste("factor", 1:3, sep = "_"), year %in% seq(1990,2015, 5))
Epsilon2_sf <- filter(Epsilon2_sf, factor %in% paste("factor", 1:3, sep = "_"), year %in% seq(1990,2015, 5))

Omega1_sf$x2i <- 1:250
Omega2_sf$x2i <- 1:250

DF2 <- MapDetails_List[['PlotDF']]

DF3_O1 <- merge(DF2, Omega1_sf)
DF3_O2 <- merge(DF2, Omega2_sf)
DF3_E1 <- merge(DF2, Epsilon1_sf)
DF3_E2 <- merge(DF2, Epsilon2_sf)

DF3_O1 <- DF3_O1[DF3_O1$Include==T,]
DF3_O2 <- DF3_O2[DF3_O2$Include==T,]
DF3_E1 <- DF3_E1[DF3_E1$Include==T,]
DF3_E2 <- DF3_E2[DF3_E2$Include==T,]

colnames(DF3_O1)[5:13] <- colnames(DF3_O2)[5:13] <- paste("factor",1:9, sep = '_')

DF4_O1 <- melt(DF3_O1, id = c("x2i", "Lat","Lon","Include")); colnames(DF4_O1)[6] <- "Spatial Encounter Prob"
DF4_O2 <- melt(DF3_O2, id = c("x2i", "Lat","Lon","Include")); colnames(DF4_O2)[6] <- "Spatial Density"
DF4_E1 <- DF3_E1; colnames(DF4_E1)[7] <- 'Spatio-temporal Encounter Prob'
DF4_E2 <- DF3_E2; colnames(DF4_E2)[7] <- 'Spatio-temporal Density'

cols <- brewer.pal(8, 'Set1')

#p1 <- ggplot() + geom_point(data = filter(DF4_O1, variable %in% paste("factor", 1:3, sep="_")), aes(x = Lon, y =  Lat, colour = value), size = 1) +
#			    geom_point(data = filter(DF4_O1, variable %in% paste("factor", 1:3, sep="_")), aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
#    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
#  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +
#  facet_grid(variable ~.)

#p2 <- ggplot() + geom_point(data = filter(DF4_O2, variable %in% paste("factor", 1:3, sep="_")), aes(x = Lon, y =  Lat, colour = value), size = 1) +
#			    geom_point(data = filter(DF4_O2, variable %in% paste("factor", 1:3, sep="_")), aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
#    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
#  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +
#  facet_wrap(~variable, ncol = 1)


DF5 <- merge(DF4_O1, DF4_O2)
DF5 <- melt(DF5, id = c("x2i", "Lat", "Lon", "Include","variable"))
colnames(DF5)[c(5,6)] <- c("factor", "parameter")

ggplot() + geom_point(data = filter(DF5, factor %in% paste("factor", 1:3, sep="_")), aes(x = Lon, y =  Lat, colour = value), size = 1) +
			    geom_point(data = filter(DF5, factor %in% paste("factor", 1:3, sep="_")), aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +
  facet_grid(factor~parameter)

ggsave(file = file.path('..', 'figures', 'Figure 2 - SpatialFactorLoadingsOmega1Omega2.png'), width = 8, height = 12)

DF6 <- merge(DF4_E1, DF4_E2)
DF6 <- melt(DF6, id = c("x2i", "Lat", "Lon", "Include","factor","year"))
colnames(DF6)[c(7)] <- c("parameter")


ggplot() + geom_point(data = filter(DF6, parameter == 'Spatio-temporal Encounter Prob', year %in% c(2000,2005,2010,2015)), aes(x = Lon, y =  Lat, colour = value), size = 1) +
			    geom_point(data = filter(DF6, parameter == 'Spatio-temporal Encounter Prob', year %in% c(2000,2005,2010,2015)), aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +   facet_grid(year ~ factor)

ggsave(file = file.path('..', 'figures', 'Suppl - SpatioTempLoadingsEpsilon1Pres.png'), width = 6, height = 6)

ggplot() + geom_point(data = filter(DF6, parameter == 'Spatio-temporal Density'), aes(x = Lon, y =  Lat, colour = value), size = 1) +
			    geom_point(data = filter(DF6, parameter == 'Spatio-temporal Density'), aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") +   facet_grid(year ~ factor)

ggsave(file = file.path('..', 'figures', 'Suppl - SpatioTempLoadingsEpsilon2.png'), width = 8, height = 12)


