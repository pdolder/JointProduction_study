
# Install TMB
# Must be installed from: https://github.com/kaskr/adcomp

# Install INLA
# Must be installed from: http://www.r-inla.org/download

# Install geostatistical delta-GLMM package
if(!"VAST" %in% installed.packages()[,1]) devtools::install_github("james-thorson/VAST")
if(!"ThorsonUtilities" %in% installed.packages()[,1]) devtools::install_github("james-thorson/utilities")
# if(!"FishData" %in% installed.packages()[,1]) devtools::install_github("james-thorson/FishData")

# setwd("C:/Users/James.Thorson/Desktop/Project_git/VAST/examples/")

# Load libraries
library(TMB)
library(ThorsonUtilities)
library(VAST)

run <- 'RUN_SIMPLE_BIG'

# This is where all runs will be located
DateFile  <- file.path('..','results',paste(Sys.Date(),'_',run,'/', sep = ""))

dir.create(DateFile)

###############
# Settings
###############

  Version = "VAST_v1_8_0"
  Method = c("Grid", "Mesh")[2]
  n_x = c(100, 250, 500, 1000, 2000)[1] # Number of stations
  FieldConfig = c("Omega1"=6, "Epsilon1"=6, "Omega2"=6, "Epsilon2"=6) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  ObsModel = c(2,0)  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
  OverdispersionConfig = c("eta1" = 6,"eta2" = 6) # 0 - number of factors
  Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
  BiasCorr = FALSE 

  # Determine region
  Region = "Celtic_Sea"
  Catch_units <- 'Kg'
  max_dist <- 50


  # Decide on strata for use when calculating indices
  strata.limits <- data.frame('STRATA'="All_areas")

  # Save options for future records
  Record = ThorsonUtilities::bundlelist( c("Version","Method","n_x","FieldConfig","RhoConfig", "ObsModel", "OverdispersionConfig", "Kmeans_Config","Catch_units","BiasCorr","Region","strata.limits") )
  capture.output( Record, file=paste0(DateFile,"Record.txt"))
  save(Record, file=paste0(DateFile,"Record.RData"))

################
# Prepare data
# (THIS WILL VARY FOR DIFFERENT DATA SETS) 
################

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

 sort(unique(DF$spp)) 

## At size
 DF <- DF[(DF$spp %in% c('Gadus morhua_Juv',
			 'Gadus morhua_Adu',
			 'Melanogrammus aeglefinus_Juv',
			 'Melanogrammus aeglefinus_Adu',
		         'Merlangius merlangus_Juv',		 
			 'Merlangius merlangus_Adu')),]

# Trim years and surveys
  DF <- DF[(DF$Year %in% c(1992:2015)),]
#  DF <- DF[(!DF$Survey %in% c('CARLHELMAR', 'NWGFS')),]
  
  ## Add missing zeros
#  DF <- reshape2::dcast(DF, Survey + Year + Station + Lat + Lon + AreaSwept_km2 ~ spp, value.var = 'Kg', fill = 0)
#  DF <- reshape2::melt(DF, id = c('Survey','Year','Station','Lat','Lon','AreaSwept_km2'), value.name = 'Kg')
  colnames(DF)[7] <- 'spp'

  DF$SpeciesName <- factor(DF$spp) # drop empty factors
  DF$Ship        <- factor(DF$Survey)
  DF$Year        <- factor(DF$Year)

  an <- as.numeric
  Data_Geostat = cbind("spp"=DF[,"SpeciesName"], 
		       "Year"=DF[,"Year"], 
		       "Catch_KG"=DF[,"Kg"], 
		       "AreaSwept_km2"=DF[,'AreaSwept_km2'], 
		       "Vessel"= DF[,'Ship'] ,
		       "Lat"=DF[,"Lat"], 
		       "Lon"=DF[,"Lon"] )

## Prepare the fixed vessel covariates, Q_ik
#Vess_Cov <- vector_to_design_matrix(Data_Geostat[,'Vessel'])
#Vess_Cov <- Vess_Cov[,-1]

# Read in the habitat covariate function to generate X_xj
#source(file.path('..', 'data', 'Covariates', 'HabitatCovariateFunc.R'))

 
  # Get extrapolation data
 Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample = max_dist)

  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn(n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
  Data_Geostat = cbind( Data_Geostat, Spatial_List$loc_UTM, "knot_i"=Spatial_List$knot_i )

#Hab <- HabAssignFunc(Kmeans = Spatial_List$Kmeans, zone = 29, locationHabMap = file.path('..', 'data', 'Covariates', '201208_EUSeaMap_Atlantic_Habitats'), nameHabMap = '201208_EUSeaMap_Atlantic_Habitats')  
  
#Hab2 <- vector_to_design_matrix(Hab$Habitat)
#Hab2 <- Hab2[,-6]

################
# Make and Run TMB model
# (THIS WILL BE SIMILAR FOR EVERY DATA SET) 
################

  # Make TMB data list
  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1 , "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method )

  # Make TMB object
  #dyn.unload( paste0(DateFile,"/",dynlib(TMB:::getUserDLL())) )
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"= getwd(), "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)
  Obj = TmbList[["Obj"]]

  # Run model
  Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=BiasCorr )
  Report = Obj$report()

  # Save stuff
Save = list(Obj = Obj,"Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData, "Data_Geostat" = Data_Geostat)
 save(Save, file=paste0(DateFile,"Save.RData"))


