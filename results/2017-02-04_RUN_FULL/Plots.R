
# Load libraries
library(TMB)
library(ThorsonUtilities)
library(VAST)

# This is where all runs will be located
DateFile <- file.path(paste(getwd(),'/',sep=''))

load('Save.RData')

###############
# Settings
###############

#########################
### VAST CPP version ###
  Version = "VAST_v2_1_0"
########################
## Spatial settings ###
########################
  Method = c("Grid", "Mesh")[2]
  #grid_size_km = 20 
  n_x = c(50, 100, 250, 500, 1000, 2000)[3] # Number of stations
  Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
 strata.limits <- data.frame('STRATA'="All_areas") # Decide on strata for use when calculating indices
  Region = "Celtic_Sea"# Determine region
  Catch_units <- 'Kg'
  max_dist <- 50
########################
#### Model settings ####
########################
  FieldConfig = c("Omega1"=10, "Epsilon1"=10, "Omega2"=10, "Epsilon2"=10) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  ObsModel = c(2,0)  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
  OverdispersionConfig = c("eta1" = 0,"eta2" = 0) # 0 - number of factors
  BiasCorr = FALSE 

##############################
diag.plots <- TRUE      ## Do you want to plot diagnostics ??
######################
#### Prepare data ####
######################

  # Read or simulate trawl data
  load(file.path('..','..','data', 'Cleaned','CelticSurveyFormattedSize.RData')) ## EVHOE and IE-IGFS
  load(file.path('..','..','data', 'Cleaned','CelticSurvey2FormattedSize.RData')) ## Various Cefas surveys

  # Combine the survey data
  DF2 <- DF

  ac <- as.character
  DF <- data.frame(Survey        = c(DF2$Ship,        ac(FSS$fldSeriesName)),
		   Year          = c(DF2$Year,        ac(FSS$Year)),
		   Station       = c(DF2$StNo,        FSS$fldCruiseStationNumber),
		   Lat           = c(DF2$HaulLatMid,  FSS$HaulLatMid),
		   Lon           = c(DF2$HaulLonMid,  FSS$HaulLonMid),
		   AreaSwept_km2 = c(DF2$SweptArea,   FSS$SweptArea),
		   spp           = c(DF2$SpeciesName, ac(FSS$Species)),
		   Kg            = c(DF2$Kg,          FSS$Kg))

 table(DF$Survey, DF$Year)		   

 ## Subset years to best data - based on data exp. doc
 DF <- DF[DF$Year %in% c(1990:2015),]

 sort(unique(DF$spp)) 

  DF$SpeciesName <- factor(DF$spp) # drop empty factors
  DF$Ship        <- factor(DF$Survey)
  DF$Year        <- factor(DF$Year)

  an <- as.numeric
  Data_Geostat = data.frame("spp"=DF[,"SpeciesName"], 
		       "Year"=DF[,"Year"], 
		       "Catch_KG"=DF[,"Kg"], 
		       "AreaSwept_km2"=DF[,'AreaSwept_km2'], 
		       "Vessel"= DF[,'Ship'] ,
		       "Lat"=DF[,"Lat"], 
		       "Lon"=DF[,"Lon"] )

## Prepare the fixed vessel covariates, Q_ik
Vess_Cov <- vector_to_design_matrix(Data_Geostat[,'Vessel'])
Vess_Cov <- Vess_Cov[,-7] # Relative to WCGFS

# Read in the habitat covariate function to generate X_xj
source(file.path('..','..' ,'data', 'Covariates', 'HabitatCovariateFunc.R'))
source(file.path('..','..' ,'data', 'Covariates', 'DepthCovariateFunc.R'))


##############################
##### Extrapolation grid #####
##############################
  # Get extrapolation data
 Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample = max_dist)

  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn(n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )

#### Prep data
Data_Geostat = cbind(Data_Geostat, Spatial_List$loc_i, "knot_i"=Spatial_List$knot_i)

#### Assign habitat class
#Hab <- HabAssignFunc(Kmeans = Spatial_List$Kmeans, zone = 30, locationHabMap = file.path('..','..' ,'data', 'Covariates', '201208_EUSeaMap_Atlantic_Habitats'), nameHabMap = '201208_EUSeaMap_Atlantic_Habitats')  
Hab <- HabAssignFunc(Kmeans = Spatial_List$Kmeans, zone = 30, locationHabMap = file.path('..', '..', 'data', 'Covariates', 'CelticSeaMap'), nameHabMap = 'CelticSeaMap') 
 
Hab2 <- vector_to_design_matrix(Hab$Habitat)
Hab2 <- Hab2[,-6] # Relative to ???

## Assign the depth                                                                                                 
Depths <- DepthAssignFunc(Kmeans = Spatial_List$Kmeans, zone = 30, locationDepths = file.path('..', '..','data', 'Covariates', 'Bathy.RData'))
## Combine the covariates
Qs <- cbind(Hab2, Depths)                                                                                               

if(diag.plots == T) {
  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))

  ## Plot encounter probability 
  Enc_prob <- SpatialDeltaGLMM::Check_encounter_prob(Report = Save$Report, Data_Geostat = Data_Geostat, DirName = DateFile)

 # Positive catch rate Q-Q plot
  Q = SpatialDeltaGLMM::QQ_Fn( TmbData=Save$TmbData, Report=Save$Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

# Get region-specific settings for plots
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn("Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List" = Extrapolation_List)

  ## Plot residuals
  SpatialDeltaGLMM:::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=Save$TmbData, Report=Save$Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

  # Plot Anisotropy  
  SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Save$Report, TmbData=Save$TmbData )

  # Plot covariances
  Cov_List = Summarize_Covariance( Report=Save$Report, ParHat=Save$ParHat, Data=Save$TmbData, SD=Save$Opt$SD, plot_cor=TRUE, category_names=levels(DF[,'SpeciesName']), figname=paste0("Spatio-temporal_covariances"), plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )

  TmbData <- Save$TmbData
  # Plot overdispersion
  Plot_Overdispersion( filename1=paste0(DateFile,"Overdispersion"), filename2=paste0(DateFile,"Overdispersion--panel"), Data=Save$TmbData, ParHat=Save$ParHat, Report=Save$Report, ControlList1=list("Width"=5, "Height"=10, "Res"=600, "Units"='in'), ControlList2=list("Width"=TmbData$n_c, "Height"=TmbData$n_c, "Res"=200, "Units"='in'))

  # Plot index
  SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=Save$TmbData, Sdreport=Save$Opt$SD, Year_Set=seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year'])))), strata_names=strata.limits[,1], category_names=levels(DF[,'SpeciesName']), use_biascorr=BiasCorr, cex = 0.3)


 # Plot surface - this is for all spp
  Dim = c( "Nrow"=ceiling(sqrt(length(Years2Include))), "Ncol"=ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))))
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1:3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Save$Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateFile,"Field_"), Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], category_names=levels(DF[,'SpeciesName']), mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8, Legend=MapDetails_List[["Legend"]])


## Plot factors

fac <- Plot_factors(Report = Save$Report, ParHat = Save$ParHat, Data = Save$TmbData, SD = Save$Opt$SD, mapdetails_list = MapDetails_List, Year_Set = Year_Set, category_names = levels(DF[,'SpeciesName']), plotdir = DateFile)

PCAstyle <- TRUE

if(PCAstyle == TRUE) {
category_names = levels(DF[,'SpeciesName'])
L_list = vector('list', length = 4)
names(L_list) = c("Omega1", "Epsilon1", "Omega2", "Epsilon2")
Data = Save$TmbData
ParHat = Save$ParHat
Report = Save$Report
i <- 'Epsilon2'
Par_name = 'Epsilon2'
Var_name = 'Epsiloninput2_sft'
L_list[[i]] <- calc_cov(L_z = ParHat[[paste0('L_',tolower(Par_name), '_z')]], n_f = Data[['FieldConfig']][[Par_name]], n_c = Data$n_c, returntype = 'loadings_matrix')
rownames(L_list[[i]]) <- category_names
Var_rot = SpatialDFA::Rotate_Fn(L_pj = L_list[[i]], Psi = Report[[Var_name]], RotationMethod = 'PCA', testcutoff = 1e-04)

rownames(Var_rot$L_pj_rot) <- paste(rep(c('cod','meg','bud','pis','had','whg','hke','ple','sol'), each = 2), c('adu','juv'), sep = "_")

plot(Var_rot$L_pj_rot[,1], Var_rot$L_pj_rot[,2], xlim = c(-1,1), ylim = c(-1,1), pch = 16, xlab = 'Factor 1', ylab = 'Factor 2')
abline(v = 0)
abline(h = 0)
arrows(0,0,Var_rot$Eigen$vectors[1,1], Var_rot$Eigen$vectors[1,2], col = 'red', length = 0.1)
text(Var_rot$Eigen$vectors[1,1], Var_rot$Eigen$vectors[1,2], label = 'Loading 1', col = 'red', cex = 0.8)
arrows(0,0,Var_rot$Eigen$vectors[2,1], Var_rot$Eigen$vectors[2,2], col = 'red', length = 0.1)
text(Var_rot$Eigen$vectors[2,1], Var_rot$Eigen$vectors[2,2], label = 'Loading 2', col = 'red', cex = 0.8)
text(Var_rot$L_pj_rot[,1], Var_rot$L_pj_rot[,2], label = rownames(Var_rot$L_pj_rot))

}


plot.cov <- TRUE


if(plot.cov == TRUE) {
HabDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers),
		    Habitat = Hab$Habitat, Depth = Depths)

library(ggplot2); library(cowplot)
p1 <- ggplot(HabDF, aes(x = x, y = y)) + geom_point(aes(colour = Habitat)) + theme_classic() + theme(legend.position = 'top')
p2 <- ggplot(HabDF, aes(x = x, y = y))  + geom_point(aes(colour = Depth)) + theme(legend.position = 'top') +
       scale_colour_gradient2(low = "black", mid = "grey90",   high = "blue", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")

plot_grid(p1,p2)
ggsave(file = 'CovariatesPredUpdated.png', width = 10, height = 6)

# Count of vessel covariates
VessDF <- data.frame(group = paste(Data_Geostat[,'Vessel'],Data_Geostat[,'spp'], sep = '_'), value = Data_Geostat[,'Catch_KG'])
VessDF$Zeros <- ifelse(VessDF$value == 0, 'Zero', 'Positive')

VessDF <- as.matrix(table(VessDF$group, VessDF$Zeros))
VessDF <- data.frame(group = rownames(VessDF), Positive = VessDF[,'Positive'], Zeros = VessDF[,'Zero'], Tot = rowSums(VessDF[,1:2]))
VessDF[,c('Positive','Zeros')] <- VessDF[,c('Positive','Zeros')]/VessDF[,c('Tot')]
VessDF <- reshape2::melt(VessDF, id = c('group'))

ggplot(VessDF[VessDF$variable != 'Tot',], aes(x = group, y = value)) + geom_bar(stat = 'identity', aes(fill = variable)) + coord_flip() +
	theme(axis.text = element_text(size = 4)) + ylab('Percentage zeros')
ggsave(file = 'CovariatesQ.png', width = 14, height = 6)

### Plot habitat covariates vs data

## add hab to data

HabCov <- Data_Geostat
HabCov$Habitat <- HabDF$Habitat[match(HabCov$knot_i, HabDF$knot)]
HabCov$Habitat <- factor(HabCov$Habitat)
HabCov$Survey_spp <- paste(HabCov$Vessel, HabCov$spp)

library(dplyr) ; library(tidyr)

HabCovsppAgg <- group_by(HabCov, spp, Habitat) %>% summarise(Catch_KG = sum(Catch_KG)) %>% complete(spp, Habitat) %>% as.data.frame()
HabCovsppAgg[is.na(HabCovsppAgg$Catch_KG),]

## No missing spp - habitat combinations

HabCovsppAgg <- group_by(HabCov, Vessel, Habitat) %>% summarise(Catch_KG = sum(Catch_KG)) %>% complete(Vessel, Habitat) %>% as.data.frame()
HabCovsppAgg[is.na(HabCovsppAgg$Catch_KG),]

## Carlhelmar - seabed missing
## CEXP - rock or other hard substrate missing


ggplot(HabCov, aes(x = spp, y = Catch_KG)) + geom_boxplot() + facet_wrap(~Habitat) + theme(axis.text.x = element_text(angle = -90, hjust = 0))
ggsave(file = 'Covariates_Hab_sppUpdated.png', width = 14, height = 6)
ggplot(HabCov, aes(x = Vessel, y = Catch_KG)) + geom_boxplot() + facet_wrap(~Habitat) + theme(axis.text.x = element_text(angle = -90, hjust = 0))
ggsave(file = 'Covariates_Hab_vessUpdated.png', width = 14, height = 6)

}









}



}
###############




