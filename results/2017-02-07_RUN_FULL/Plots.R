
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
  n_x = c(50, 100, 250, 500, 1000, 2000)[2] # Number of stations
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
 DF <- DF[DF$Year %in% c(2000:2015),]

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

##############################
##### Extrapolation grid #####
##############################
  # Get extrapolation data
  Extrapolation_List <- Save$Extrapolation_List
  Spatial_List <- Save$Spatial_List


#### Prep data
Data_Geostat = cbind(Data_Geostat, Spatial_List$loc_i, "knot_i"=Spatial_List$knot_i)

if(diag.plots == T) {
  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))

  ## Plot encounter probability 
  Enc_prob <- SpatialDeltaGLMM::Check_encounter_prob(Report = Save$Report, Data_Geostat = Data_Geostat, DirName = DateFile)

 # Positive catch rate Q-Q plot
  Q = SpatialDeltaGLMM::QQ_Fn( TmbData=Save$TmbData, Report=Save$Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

# Get region-specific settings for plots
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn("Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap,
						    "Extrapolation_List" = Extrapolation_List)

  ## Plot residuals
  SpatialDeltaGLMM:::plot_residuals(Lat_i=Data_Geostat[,'Lat'], Lon_i=Data_Geostat[,'Lon'], TmbData=TmbData, Report=Report, Q=Q, savedir=DateFile, MappingDetails=MapDetails_List[["MappingDetails"]], PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=DateFile, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8)

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

Plot_factors(Report = Save$Report, ParHat = Save$ParHat, Data = Save$TmbData, SD = Save$Opt$SD, mapdetails_list = MapDetails_List, Year_Set = Year_Set, category_names = levels(DF[,'SpeciesName']), plotdir = DateFile)
}
###############


cov <- data.frame(cov = rep(colnames(Save$Vessel_Covariate), times = 2),
		  MLE_enc = params$MLE[params$Param == 'lambda1_k'],
		  MLE_pos = params$MLE[params$Param == 'lambda2_k'])

cov$Survey  <- sapply(strsplit(as.character(cov$cov), '_'), '[', 1) 
cov$Species <- sapply(strsplit(as.character(cov$cov), '_'), '[', 2) 
cov$Stage <- sapply(strsplit(as.character(cov$cov), '_'), '[', 3) 

sd_report <- Save$Opt$SD

library(ggplot2)

ggplot(cov, aes(x = Survey, y = MLE_enc)) + geom_point(aes(colour = Stage)) +
	facet_wrap(~Species) + theme(axis.text.x = element_text(angle = -90,hjust = 0))

ggplot(cov, aes(x = Survey, y = exp(MLE_pos))) + geom_point(aes(colour = Stage)) +
	facet_wrap(~Species) + theme(axis.text.x = element_text(angle = -90,hjust = 0))


