##################################################################
##### Combined figure of spatial and species factor loadings #####
##################################################################

library(VAST)
library(ggplot2) 
library(cowplot)
library(ggthemes)
library(mapplots)
library(rworldmap)
library(rworldxtra)
library(data.table)
library(broom)
library(RColorBrewer)
library(dplyr) 


run <- '2017-06-16_M1'
load(file = file.path("..","..",'..', '..', '..','results', run, 'Save.RData'))

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
load(file.path("..","..",'..', '..', '..','results', 'CovariatesAtKnot.RData'))

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
Epsilon1_sf$Var2 <- Epsilon2_sf$Var2 <- rep(paste("Factor", 1:9, sep = " "), each = 250)
Epsilon1_sf$Var3 <- Epsilon2_sf$Var3 <- rep(1990:2015, each = 250 * 9)

colnames(Epsilon1_sf) <- colnames(Epsilon2_sf) <- c('x2i','factor','year','value')

## Save all spatio-temp
Ep1 <- Epsilon1_sf
Ep2 <- Epsilon2_sf

## Just the first 3 factors for the Epsilons and every 5 years, otherwise files are huge
Epsilon1_sf <- filter(Epsilon1_sf, factor %in% paste("Factor", 1:3, sep = " "), year %in% seq(1990,2015, 5))
Epsilon2_sf <- filter(Epsilon2_sf, factor %in% paste("Factor", 1:3, sep = " "), year %in% seq(1990,2015, 5))

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

colnames(DF3_O1)[5:13] <- colnames(DF3_O2)[5:13] <- paste("Factor",1:9, sep = ' ')

DF4_O1 <- melt(DF3_O1, id = c("x2i", "Lat","Lon","Include")); colnames(DF4_O1)[6] <- "(a) average spatial encounter probability"
DF4_O2 <- melt(DF3_O2, id = c("x2i", "Lat","Lon","Include")); colnames(DF4_O2)[6] <- "(b) average spatial density"
DF4_E1 <- DF3_E1; colnames(DF4_E1)[7] <- '(a) Spatiotemporal encounter probability'
DF4_E2 <- DF3_E2; colnames(DF4_E2)[7] <- '(b) Spatiotemporal density'

cols <- brewer.pal(8, 'Set1')

DF5 <- merge(DF4_O1, DF4_O2)
DF5 <- melt(DF5, id = c("x2i", "Lat", "Lon", "Include","variable"))
colnames(DF5)[c(5,6)] <- c("factor", "parameter")

####################

### Plot the temporal trends

Ep1 <- filter(Ep1, factor %in% c(paste('Factor', 1:3, sep = " ")))

Ep1med <- Ep1 %>% group_by(factor, year) %>% summarise(value = median(value))

ggplot(Ep1,aes(x = year, y = value)) + geom_line(colour = 'grey90', aes(group = x2i)) +
	facet_wrap(~factor, ncol = 1) + geom_point(data = Ep1med, aes(x = year, y = value)) +
	geom_line(data = Ep1med, aes(x = year, y = value))


###############################
### Species factor loadings ###
###############################

PCA.DFs <- lapply(Var_list, function(x) {
			  DF <- as.data.frame(x$L_pj_rot)
			  colnames(DF) <- paste('Factor',1:9, sep = ' ')
			  return(DF)
	     })

## explained variance
## sum(L_pj_rot[f]^2)/sum(L_pj_rot^2)

params <- names(Var_list)

var.DF <- lapply(params, function(y) {
DF <- data.frame(variable = y, 
		 factor = paste("Factor", 1:9, sep = " "),
		 values = sapply(1:9, function(x) {round(100 * sum(Var_list[[y]]$L_pj_rot[,x]^2) / sum(Var_list[[y]]$L_pj_rot^2),1)})
)
		 return(DF)
})

var.DF <- do.call(rbind, var.DF)

## Omega 1

O1 <- data.frame(species = rownames(PCA.DFs[["Omega1"]]), reshape2::melt(PCA.DFs[["Omega1"]]))
colnames(O1)[2] <- "factor"

# Colour for the labels and bars
O1$col <- ifelse(O1$value > 0, "blue", "red")

## Omega 2

O2 <- data.frame(species = rownames(PCA.DFs[["Omega2"]]), reshape2::melt(PCA.DFs[["Omega2"]]))
colnames(O2)[2] <- "factor"

# Colour for the labels and bars
O2$col <- ifelse(O2$value > 0, "blue", "red")

## Epsilon 1

E1 <- data.frame(species = rownames(PCA.DFs[["Epsilon1"]]), reshape2::melt(PCA.DFs[["Epsilon1"]]))
colnames(E1)[2] <- "factor"

# Colour for the labels and bars
E1$col <- ifelse(E1$value > 0, "blue", "red")

## Epsilon 2

E2 <- data.frame(species = rownames(PCA.DFs[["Epsilon2"]]), reshape2::melt(PCA.DFs[["Epsilon2"]]))
colnames(E2)[2] <- "factor"

# Colour for the labels and bars
E2$col <- ifelse(E2$value > 0, "blue", "red")



#############################
## Plot with all combined ###
#############################

## Spatial plots 

# 3 factors, 2 variables

var <- unique(DF5$parameter)

####################################################
## Loop through to assign all plots to a variable ##
####################################################

for (i in 1:2) {

	for (j in 1:3) {
assign(paste("sp", i, "f", j, sep = ""),

ggplot() + geom_point(data = filter(DF5, parameter == var[i], factor == paste("Factor", j, sep = " ")), aes(x = Lon, y = Lat, colour = value), size = 1) + geom_point(data = filter(DF5, parameter == var[i], factor == paste("Factor", j, sep = " ")) , aes(x = Lon, y = Lat, colour = value), size = 0.5, shape = 3) +
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = 'none',plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') + theme(axis.text = element_blank()) +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") + theme_void() + theme(legend.position=  "none") + geom_text(aes(x = seq(-10,-2.5,2.5), y = 48), label = seq(-10,-2.5,2.5), size = 4, colour = "black") + geom_text(aes(x = -12.5+0.6, y = seq(48,52,1)), label = c("48/-12.5", seq(49,52,1)), size = 4, colour = "black")


)	}

}

###################
## Species plots ##
###################

mins <- round(min(rbind(O1,O2)$value)-0.1,1)
maxs <- round(max(rbind(O1,O2)$value)+0.1,1)

for (i in 1:3) {

	if(i %in% 1) {
assign(paste("Spp",1,"f",i,sep = ""),
ggplot(filter(O1, factor == paste("Factor", i, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col, alpha = abs(value)), colour = "black") + ylim(mins, maxs) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none", plot.margin=unit(c(0,0,0,0),"mm")) +
	 scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") +# coord_flip() +
	ggtitle(paste("Variance explained = ", filter(var.DF, variable == "Omega1", factor == paste("Factor", i, sep = " "))$values, "%")) + theme_void() + theme(legend.position=  "none", plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
)

	}

	if(i %in% 2) {
assign(paste("Spp",1,"f",i,sep = ""),
ggplot(filter(O1, factor == paste("Factor", i, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col, alpha = abs(value)), colour = "black") + ylim(mins, maxs) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none", plot.margin=unit(c(0,0,0,0),"mm")) +
	 scale_fill_manual(values = c("blue")) + ylab("") + xlab("") +# coord_flip() +
	ggtitle(paste("Variance explained = ", filter(var.DF, variable == "Omega1", factor == paste("Factor", i, sep = " "))$values, "%")) + theme_void() + theme(legend.position=  "none", plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
)

	}


	if(i == 3) {
	
	assign(paste("Spp",1,"f",i,sep = ""),
ggplot(filter(O1, factor == paste("Factor", i, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col, alpha = abs(value)), colour = "black") + ylim(mins, maxs) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none", plot.margin=unit(c(0,0,0,0),"mm")) +
	 scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") +# coord_flip() +
	ggtitle(paste("Variance explained = ", filter(var.DF, variable == "Omega1", factor == paste("Factor", i, sep = " "))$values, "%")) + theme(legend.position=  "none", plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 14, face = "bold"))
)

	}

}

for (i in 1:3) {

	if(i %in% 1:2) {

assign(paste("Spp",2,"f",i,sep = ""),
ggplot(filter(O2, factor == paste("Factor", i, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col, alpha = abs(value)), colour = "black") + ylim(mins, maxs) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none", plot.margin = unit(c(0,0,0,0),"mm")) +
	 scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") + #coord_flip() +
	 ggtitle(paste("Variance explained = ", filter(var.DF, variable == "Omega2", factor == paste("Factor", i, sep = " "))$values, "%")) + theme_void() + theme(legend.position=  "none", plot.title = element_text(size = 12, face = "bold", hjust = 0.5))


)


	}


if(i == 3) {

assign(paste("Spp",2,"f",i,sep = ""),
ggplot(filter(O2, factor == paste("Factor", i, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col, alpha = abs(value)), colour = "black") + ylim(mins, maxs) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none", plot.margin = unit(c(0,0,0,0),"mm")) +
	 scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") + #coord_flip() +
	 ggtitle(paste("Variance explained = ", filter(var.DF, variable == "Omega2", factor == paste("Factor", i, sep = " "))$values, "%")) + theme(legend.position=  "none", plot.title = element_text(size = 14, face = "bold", hjust = 0.5), axis.ticks = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 14, face = "bold"))


)


}



}


### Now to combine it all

## The labels 

l1 <- ggplot() + geom_text(aes(1.06,1), label = "(A) average", size = 10) + theme_void()+ 
	xlim(0.9,1.1) + ylim(0.9,1.1)

l2 <- ggplot() + geom_text(aes(0.93,1), label = "encounter", size = 10) + theme_void()+ 
	xlim(0.9,1.1) + ylim(0.9,1.1)

l3 <- ggplot() + geom_text(aes(1.06,1), label = "(B) average", size = 10) + theme_void()+ 
	xlim(0.9,1.1) + ylim(0.9,1.1)

l4 <- ggplot() + geom_text(aes(0.92,1), label = "density", size = 10) + theme_void()+ 
	xlim(0.9,1.1) + ylim(0.9,1.1)

l5 <- ggplot() + geom_text(aes(1,1), label = "") + theme_void()

l6 <- ggplot() + geom_text(aes(1,1), label = "Factor 1", angle = -90, size = 10) + theme_void() + 
	xlim(0.9,1.1) + ylim(0.9,1.1)

l7 <- ggplot() + geom_text(aes(1,1), label = "Factor 2", angle = -90, size = 10) + theme_void()
l8 <- ggplot() + geom_text(aes(1,1), label = "Factor 3", angle = -90, size = 10) + theme_void()

# 5 X 4
## Order:

comb_plot <- plot_grid(l1, l2, l3, l4, l5,
	  Spp1f1, sp1f1, sp2f1, Spp2f1, l6, 
	  Spp1f2, sp1f2, sp2f2, Spp2f2, l7,
	  Spp1f3, sp1f3, sp2f3, Spp2f3, l8,
	  ncol = 5, nrow = 4,
	  rel_heights = c(0.2,1,1,1),
	  rel_widths = c(1,1,1,1,0.2
			 ))

#save_plot("Fig1.eps", comb_plot, ncol = 5, nrow = 3, device = cairo_ps)
save_plot("Fig1.png", comb_plot, ncol = 5, nrow = 3)

