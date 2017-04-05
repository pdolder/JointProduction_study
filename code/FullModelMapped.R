
# Install TMB
# Must be installed from: https://github.com/kaskr/adcomp

# Install INLA
# Must be installed from: http://www.r-inla.org/download

#If Install geostatistical delta-GLMM package
if(!"VAST" %in% installed.packages()[,1]) devtools::install_github("james-thorson/VAST")
if(!"ThorsonUtilities" %in% installed.packages()[,1]) devtools::install_github("james-thorson/utilities")

# Load libraries
library(TMB)
library(ThorsonUtilities)
library(VAST)

library(INLA)
#INLA:::inla.dynload.workaround() 

run <- 'RUN_DIAG_MAPPED_HABFIX'

# This is where all runs will be located
DateFile  <- file.path('..','results',paste(Sys.Date(),'_',run,'/', sep = ""))

dir.create(DateFile)

###############
# Settings
###############

#########################
### VAST CPP version ###
  Version = "VAST_v2_4_0"
########################
## Spatial settings ###
########################
  Method = c("Grid", "Mesh")[2]
  #grid_size_km = 20 
  n_x = c(10, 50, 100, 250, 500, 1000, 2000)[4] # Number of stations
  Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
 strata.limits <- data.frame('STRATA'="All_areas") # Decide on strata for use when calculating indices
  Region = "Celtic_Sea"# Determine region
  Catch_units <- 'Kg'
  max_dist <- 50
########################
#### Model settings ####
########################

  FieldConfig = c("Omega1"= 9, "Epsilon1"= 9, "Omega2"= 9, "Epsilon2"= 9) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  ObsModel = c(2,0)  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
  OverdispersionConfig = c("eta1" = 0,"eta2" = 0) # 0 - number of factors
  BiasCorr = FALSE 
#######################
##### Save options ####
#######################
  # Save options for future records
  Record = ThorsonUtilities::bundlelist( c("Version","Method","n_x","FieldConfig","RhoConfig", "ObsModel", "OverdispersionConfig", "Kmeans_Config","Catch_units","BiasCorr","Region","strata.limits") )
  capture.output( Record, file=paste0(DateFile,"Record.txt"))
  save(Record, file=paste0(DateFile,"Record.RData"))

  diag.plots <- FALSE  ## Do you want to plot diagnostics ??
######################
#### Prepare data ####
######################

  # Read or simulate trawl data
  load(file.path('..','data', 'Cleaned','CelticSurveyFormattedSize.RData')) ## EVHOE and IE-IGFS
  load(file.path('..','data', 'Cleaned','CelticSurvey2FormattedSize.RData')) ## Various Cefas surveys

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

 species <- sort(unique(DF$spp)) 

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

## try species only 
#spp <- sapply(strsplit(as.character(Data_Geostat[,'spp']), '\\_'), '[',1)

#sz  <- sapply(strsplit(as.character(Data_Geostat[,'spp']), '\\_'), '[',2)

## 0. Vessel only 
#Vess_Cov <- vector_to_design_matrix(paste(Data_Geostat[,'Vessel']))
 
## 1. Vessel and species concatenated 
Vess_Cov <- vector_to_design_matrix(paste(Data_Geostat[,'Vessel'],Data_Geostat[,'spp'], sep = '_'))

# 1a. Drop only single vessel-species combo
# Vess_Cov  <- Vess_Cov[,-1]

# 1b. Drop set of vessel-species combos
#Vess_Cov <- Vess_Cov[,-c(1:18)] # - All relative to Carlhelmar for same spp
Vess_Cov <- Vess_Cov[,-grep('CEXP', colnames(Vess_Cov))]  ## Relative to the Celtic Explorer for each spp

# 2. Vessel and species separately
#Vess_Cov <- vector_to_design_matrix(Data_Geostat[,'Vessel'])
#spp_Cov <- vector_to_design_matrix(Data_Geostat[,'spp'])
#Vess_Cov <- cbind(Vess_Cov[,-1],spp_Cov[,-1])

# 2b. Vessel and species separated, interacting
# use the design formula function in R

#vess <- factor(Data_Geostat[,'Vessel'])
#sp   <- factor(Data_Geostat[,'spp'])
#Vess_Cov <- model.matrix(~ vess)

# remove the column name
#colnames(Vess_Cov) <- sub('vess', '', colnames(Vess_Cov))
#colnames(Vess_Cov) <- sub('sp', '', colnames(Vess_Cov))


#Vess_Cov <- Vess_Cov[,-1]  ## Drop the intercept


# 3. Vessel / species excluding size
#Vess_Cov <- vector_to_design_matrix(Data_Geostat[,'Vessel'])
#spp_Cov <- vector_to_design_matrix(spp)
#Vess_Cov <- cbind(Vess_Cov[,-1],spp_Cov[,-1])

# Read in the habitat covariate function to generate X_xj
source(file.path('..', 'data', 'Covariates', 'HabitatCovariateFunc.R'))
source(file.path('..', 'data', 'Covariates', 'DepthCovariateFunc.R'))

##############################
##### Extrapolation grid #####
##############################
  # Get extrapolation data
 Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample = max_dist)

  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn(n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )

#### Prep data
Data_Geostat = cbind(Data_Geostat, Spatial_List$loc_i, "knot_i"=Spatial_List$knot_i)

#### Assign habitat class

# Fix for point supposedly on land
Centers <- Spatial_List$Kmeans$centers
Centers[191,'N_km'] <- Centers[191,'N_km'] - 20

Hab <- HabAssignFunc(KmeansCenters = Centers, zone = 29, locationHabMap = file.path('..', 'data', 'Covariates', 'CelticSeaMap'), nameHabMap = 'CelticSeaMap')  
# 
Hab2 <- vector_to_design_matrix(Hab$Habitat)
Hab2 <- Hab2[,-6] # Relative to ???

## Assign the depth                                                                                                 
Depths <- unlist(DepthAssignFunc(KmeansCenters = Centers, zone = 29, locationDepths = file.path('..', 'data', 'Covariates', 'Bathy.RData')))
### Combine the covariates
StatCov <- cbind(Hab2, Depths)                                                                                               
###
HabDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers),
		    Habitat = Hab$Habitat, Depth = Depths)

library(ggplot2); library(cowplot)
p1 <- ggplot(HabDF, aes(x = x, y = y)) + geom_point(aes(colour = Habitat)) + theme_classic() + theme(legend.position = 'top')
p2 <- ggplot(HabDF, aes(x = x, y = y))  + geom_point(aes(colour = Depth)) + theme(legend.position = 'top') +
       scale_colour_gradient2(low = "black", mid = "grey90",   high = "blue", midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")

plot_grid(p1,p2)

################################
#### Make and Run TMB model ####
################################
  # Make TMB data list

## Note - removed the habitat covariates for testing !!!
  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "Q_ik" = Vess_Cov,"s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=as.numeric(Data_Geostat[,'Year']), "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method)

  # Make TMB object
#  dyn.unload( paste0(DateFile,"/",dynlib(TMB:::getUserDLL())))
  setwd(DateFile) # so executables go to the right place...

  #####################################
  ## Trying externally estimated Qs ###
  #####################################

  ## loading in values
  load(file.path('../','QEstimates.RData'))

  # Define parameter list
  Params <- Param_Fn(Version = Version, DataList = TmbData, RhoConfig = RhoConfig) # - list of parameters

  # extract and assign lambda k values
  lam1_k <- unlist(lapply(Qs, function(x) x$par.fixed[names(x$par.fixed) == 'lambda1_k']))
  lam2_k <- unlist(lapply(Qs, function(x) x$par.fixed[names(x$par.fixed) == 'lambda2_k']))

  Params$lambda1_k <- lam1_k
  Params$lambda2_k <- lam2_k

  # cod juv, bud juv should be fixed at their adult values 
  pos_juv <- which(colnames(Vess_Cov) %in% c('CARLHELMAR_Gadus morhua_Juv','CARLHELMAR_Lophius budegassa_Juv'))
  pos_adu <- which(colnames(Vess_Cov) %in% c('CARLHELMAR_Gadus morhua_Adu','CARLHELMAR_Lophius budegassa_Adu'))
 
  Params$lambda1_k[pos_juv] <- Params$lambda1_k[pos_adu]
  Params$lambda2_k[pos_juv] <- Params$lambda2_k[pos_adu]


  map <- Make_Map(TmbData = TmbData, TmbParams = Params, CovConfig = TRUE, Q_Config = TRUE,RhoConfig = RhoConfig) # make the map
 
 # Replace the values for carlhelmar cod/bud juv  
  map$lambda1_k[] <- NA 
  map$lambda1_k <- factor(map$lambda1_k)
  map$lambda2_k[] <- NA
  map$lambda2_k <- factor(map$lambda2_k)

  # Relevel
  levels(map$lambda1_k) <- map$lambda1_k
  levels(map$lambda2_k) <- map$lambda2_k

  # Parameters with fixed gear estimates and a map
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_xi,
  "Parameters" = Params, "Map" = map)
  Obj = TmbList[["Obj"]]

  # Run model
  Opt = TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0('../',DateFile), bias.correct=BiasCorr) 
  Report = Obj$report()

#################################  
######## Save outputs ###########
#################################
Save = list(Obj = Obj,"Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData, "Data_Geostat" = Data_Geostat)
save(Save, file=paste0('../',DateFile,"Save.RData"))

##########################
# Make diagnostic plots ##
##########################

if(diag.plots == T) {
  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))

  ## Plot encounter probability 
  Enc_prob <- SpatialDeltaGLMM::Check_encounter_prob(Report = Save$Report, Data_Geostat = Data_Geostat, DirName = file.path('../',DateFile))

 # Positive catch rate Q-Q plot
  Q = SpatialDeltaGLMM::QQ_Fn( TmbData=Save$TmbData, Report=Save$Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

# Get region-specific settings for plots
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList)

  ## Plot residuals
  SpatialDeltaGLMM:::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

  # Plot Anisotropy  
  SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0('../',DateFile,"Aniso.png"), Report=Save$Report, TmbData=Save$TmbData )

  # Plot covariances
  Cov_List = Summarize_Covariance( Report=Save$Report, ParHat=Save$ParHat, Data=Save$TmbData, SD=Save$Opt$SD, plot_cor=TRUE, category_names=levels(DF[,'SpeciesName']), figname=paste0("Spatio-temporal_covariances"), plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )

  TmbData <- Save$TmbData
  # Plot overdispersion
  Plot_Overdispersion( filename1=paste0(DateFile,"Overdispersion"), filename2=paste0(DateFile,"Overdispersion--panel"), Data=Save$TmbData, ParHat=Save$ParHat, Report=Save$Report, ControlList1=list("Width"=5, "Height"=10, "Res"=600, "Units"='in'), ControlList2=list("Width"=TmbData$n_c, "Height"=TmbData$n_c, "Res"=200, "Units"='in'))

  # Plot index
  SpatialDeltaGLMM::PlotIndex_Fn( DirName=file.path('../', DateFile), TmbData=Save$TmbData, Sdreport=Save$Opt$SD, Year_Set=seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year'])))), strata_names=strata.limits[,1], category_names=levels(DF[,'SpeciesName']), use_biascorr=BiasCorr, cex = 0.3)


 # Plot surface - this is for all spp
  Dim = c( "Nrow"=ceiling(sqrt(length(Years2Include))), "Ncol"=ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))))
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1:3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Save$Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0('../',DateFile,"Field_"), Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], category_names=levels(DF[,'SpeciesName']), mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8, Legend=MapDetails_List[["Legend"]])


## Plot factors

Plot_factors(Report = Save$Report, ParHat = Save$ParHat, Data = Save$TmbData, SD = Save$Opt$SD, mapdetails_list = MapDetails_List, Year_Set = Year_Set, 
	     category_names = levels(DF[,'Species']), plotdir = file.path('../',DateFile))
}
###############

