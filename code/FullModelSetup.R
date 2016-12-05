
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

run <- 'RUN_FULL'

# This is where all runs will be located
DateFile <- file.path(getwd(),paste(Sys.Date(),'_',run,'/', sep = ""))
dir.create(DateFile)

###############
# Settings
###############

  Version = "VAST_v1_8_0"
  Method = c("Grid", "Mesh")[2]
  #grid_size_km = 20 
  n_x = c(100, 250, 500, 1000, 2000)[1] # Number of stations
  FieldConfig = c("Omega1"=4, "Epsilon1"=5, "Omega2"=5, "Epsilon2"=6) # 1=Presence-absence; 2=Density given presence; #Epsilon=Spatio-temporal; #Omega=Spatial
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0) # Structure for beta or epsilon over time: 0=None (default); 1=WhiteNoise; 2=RandomWalk; 3=Constant
  ObsModel = c(2,0)  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
  OverdispersionConfig = c("eta1" = 0,"eta2" = 0) # 0 - number of factors
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
  load(file.path('..','data','CelticSurveyFormattedSize.RData')) ## EVHOE and IE-IGFS
  load(file.path('..','data','CelticSurvey2FormattedSize.RData')) ## Various Cefas surveys

  # Combine the survey data
  DF$SpeciesName <- toupper(DF$SpeciesName)
  FSS$fldScientificName <- toupper(FSS$fldScientificName)

  DF2 <- DF

  ac <- as.character
  DF <- data.frame(Survey        = c(DF2$Ship,        ac(FSS$fldSeriesName)),
		   Year          = c(DF2$Year,        ac(FSS$Year)),
		   Station       = c(DF2$StNo,        FSS$fldCruiseStationNumber),
		   Lat           = c(DF2$HaulLatMid,  FSS$HaulLatMid),
		   Lon           = c(DF2$HaulLonMid,  FSS$HaulLonMid),
		   AreaSwept_km2 = c(DF2$SweptArea,   FSS$SweptArea),
		   spp           = c(DF2$SpeciesName, ac(FSS$fldScientificName)),
		   Kg            = c(DF2$Kg,          FSS$Kg))

 table(DF$Survey, DF$Year)		   

 sort(unique(DF$spp)) 

#  DF <- DF[(DF$spp %in% c('GADUS MORHUA',
#			  'MELANOGRAMMUS AEGLEFINUS',
#			  'MERLANGIUS MERLANGUS')),]
			 # 'MERLUCCIUS MERLUCCIUS',
			 # 'LOPHIUS PISCATORIUS',
			 # 'LEPIDORHOMBUS WHIFFIAGONIS')),]

## At size
 DF <- DF[(DF$spp %in% c('GADUS MORHUA_JUV',
			 'GADUS MORHUA_ADU',
			 'MELANOGRAMMUS AEGLEFINUS_JUV',
			 'MELANOGRAMMUS AEGLEFINUS_ADU',
		         'MERLANGIUS MERLANGUS_JUV',		 
			 'MERLANGIUS MERLANGUS_ADU',
			 'MERLUCCIUS MERLUCCIUS_ADU',
			 'MERLUCCIUS MERLUCCIUS_JUV',
			 'LOPHIUS PISCATORIUS_ALL',
			 'LEPIDORHOMBUS WHIFFIAGONIS_ADU',
			 'LEPIDORHOMBUS WHIFFIAGONIS_JUV')),]

# Trim years
  DF <- DF[(DF$Year %in% c(1990:2015)),]
  
  ## Add missing zeros
  DF <- reshape2::dcast(DF, Survey + Year + Station + Lat + Lon + AreaSwept_km2 ~ spp, value.var = 'Kg', fill = 0)
  DF <- reshape2::melt(DF, id = c('Survey','Year','Station','Lat','Lon','AreaSwept_km2'), value.name = 'Kg')
  colnames(DF)[7] <- 'spp'

  DF$SpeciesName <- factor(DF$spp) # drop empty factors
  DF$Ship        <- factor(DF$Survey)
  DF$Year        <- factor(DF$Year)

  # Remove some surveys
  # DF <- DF[DF$Survey %in% c('CEXP','THA2','WCGFS','Q1SWIBTS','Q1SWBEAM'),]

  an <- as.numeric
  Data_Geostat = cbind("spp"=DF[,"SpeciesName"], 
		       "Year"=DF[,"Year"], 
		       "Catch_KG"=DF[,"Kg"], 
		       "AreaSwept_km2"=DF[,'AreaSwept_km2'], 
		       "Vessel"= DF[,'Ship'] ,
		       "Lat"=DF[,"Lat"], 
		       "Lon"=DF[,"Lon"] )

## Prepare the fixed vessel covariates, Q_ik
Vess_Cov <- as.matrix(Data_Geostat[,'Vessel']-1)

# Read in the habitat covariates X_xj
load(file.path('..','data','KmeansHab.RData'))
KmeanHab$Habitat <- factor(KmeanHab$Habitat)
Hab <- as.matrix(as.numeric(KmeanHab$Habitat)-1)
 
  # Get extrapolation data
  Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn( Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample = max_dist)

  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn(n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )
  Data_Geostat = cbind( Data_Geostat, Spatial_List$loc_UTM, "knot_i"=Spatial_List$knot_i )

################
# Make and Run TMB model
# (THIS WILL BE SIMILAR FOR EVERY DATA SET) 
################

  # Make TMB data list
  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], Q_ik = Vess_Cov, "X_xj" = Hab,"s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method )

  # Make TMB object
  #dyn.unload( paste0(DateFile,"/",dynlib(TMB:::getUserDLL())) )
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x)
  Obj = TmbList[["Obj"]]

  # Run model
  Opt = TMBhelper::Optimize( obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=BiasCorr )
  Report = Obj$report()

  # Save stuff
Save = list(Obj = Obj,"Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData, "Data_Geostat" = Data_Geostat)
 save(Save, file=paste0(DateFile,"Save.RData"))

################
# Make diagnostic plots
################

  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))


  # Plot Anisotropy  
  SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )

  # Plot covariances
  Cov_List = Summarize_Covariance( report=Report, parhat=Save$ParHat, tmbdata=TmbData, sd_report=Opt$SD, plot_cor=TRUE, names_set=levels(DF[,'SpeciesName']), figname=paste0(DateFile,"Spatio-temporal_covariances"), plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )

  # Plot overdispersion
  Plot_Overdispersion( filename1=paste0(DateFile,"Overdispersion"), filename2=paste0(DateFile,"Overdispersion--panel"), Data=TmbData, ParHat=Save$ParHat, Report=Report, ControlList1=list("Width"=5, "Height"=10, "Res"=600, "Units"='in'), ControlList2=list("Width"=TmbData$n_c, "Height"=TmbData$n_c, "Res"=200, "Units"='in'))

  # Plot index
  SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Opt$SD, Year_Set=seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year'])))), strata_names=strata.limits[,1], category_names=levels(DF[,'SpeciesName']), use_biascorr=BiasCorr, cex = 0.3)


 # Positive catch rate Q-Q plot
  Q = SpatialDeltaGLMM::QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

 # Plot surface - this is for all spp
  Dim = c( "Nrow"=ceiling(sqrt(length(Years2Include))), "Ncol"=ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))))
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=1:9, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateFile,"Field_"), Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], category_names=levels(DF[,'SpeciesName']), mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8, Legend=MapDetails_List[["Legend"]])

#SpatialDeltaGLMM::Plot_data_and_knots(Data_Extrap = Spatial_List$NN_Extrap, 
#				      Extrap_Area_km2 = Spatial_List$a_xl, 
#				      loc_x = Spatial_List$loc_x,
#				      Spatial_List$loc_x, 
#				      Data_Geostat = Data_Geostat,
#				      Plot_name = 'Data_and_knots.png',
#				      Data_name = 'Data_by_year.png')

#### END ###

