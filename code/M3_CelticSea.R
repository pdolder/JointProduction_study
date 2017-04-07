#############################################################
#############################################################
## VAST implementation for 18 species groups in the Celtic ##
## Sea. 
## 3 model set ups:
## M0 = Vessel as random effect, no Q or static covariates
## M1 = Vessel as fixed effect (Q covariate), no static cov
## M2 = Vessel as fixed effect, habitat and depth as static cov
## M3 = Vessel as fixed effect - estimated together, habitat depth as stat cov
############################################################
## 06 April 2017
## Paul Dolder
##############################################################

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
INLA:::inla.dynload.workaround() ## Needed on older linux machines 

# setup run
mod <- 'M3'
run <- mod 

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


###############################
### Catchability covariates ##
## Prepare the fixed vessel covariates, Q_ik

if(mod %in% c('M1', 'M2', 'M3')) {
## Vessel and species concatenated 
Vess_Cov <- vector_to_design_matrix(paste(Data_Geostat[,'Vessel'],Data_Geostat[,'spp'], sep = '_'))

# Drop set of vessel-species combos
Vess_Cov <- Vess_Cov[,-grep('CEXP', colnames(Vess_Cov))]  ## All spp relative to the Celtic Explorer 

}

##############################

##############################
##### Extrapolation grid #####
##############################
  # Get extrapolation data
 Extrapolation_List = SpatialDeltaGLMM::Prepare_Extrapolation_Data_Fn(Region=Region, strata.limits=strata.limits, observations_LL=Data_Geostat[,c('Lat','Lon')], maximum_distance_from_sample = max_dist)

  # Calculate spatial information for SPDE mesh, strata areas, and AR1 process
  Spatial_List = SpatialDeltaGLMM::Spatial_Information_Fn(n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile )

#### Prep data
Data_Geostat = cbind(Data_Geostat, Spatial_List$loc_i, "knot_i"=Spatial_List$knot_i)


#############################
#### Static covariates ######
#############################

if (mod %in% c('M2','M3')) {
#### Assign habitat class
# Read in the habitat covariate function to generate X_xj
source('HabitatCovariateFunc.R')
source('DepthCovariateFunc.R')

# Fix for point supposedly on land - move it down 20km
Centers <- Spatial_List$Kmeans$centers
Centers[191,'N_km'] <- Centers[191,'N_km'] - 20

Hab <- HabAssignFunc(KmeansCenters = Centers, zone = 29, locationHabMap = file.path('..', 'data', 'Covariates', 'CelticSeaMap'), nameHabMap = 'CelticSeaMap')  
# 
Hab2 <- vector_to_design_matrix(Hab$Habitat)
Hab2 <- Hab2[,-grep('Seabed', colnames(Hab2))] # Relative to seabed 

## Assign the depth                                                                                                 
Depths <- unlist(DepthAssignFunc(KmeansCenters = Centers, zone = 29, locationDepths = file.path('..', 'data', 'Covariates', 'Bathy.RData')))
### Combine the covariates
StatCov <- cbind(Hab2, Depths)                                                                                               

}

################################
#### Make and Run TMB model ####
################################
  # Make TMB data list

## Note - removed the habitat covariates for testing !!!

# M0 - no covariates, random vessel
if (mod == 'M0') {
  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=as.numeric(Data_Geostat[,'Year']), "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method)
}

# M1 - vessel as fixed effect
if (mod == 'M1') {
  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "Q_ik" = Vess_Cov,"s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=as.numeric(Data_Geostat[,'Year']), "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method)
}

# M2 - vessel as fixed effect + habitat covariates
if (mod %in% c('M2','M3')) {
 TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "Q_ik" = Vess_Cov,"X_xj" = StatCov, "s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=as.numeric(Data_Geostat[,'Year']), "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method)

}

  # Make TMB object
#  dyn.unload( paste0(DateFile,"/",dynlib(TMB:::getUserDLL())))
  setwd(DateFile) # so executables go to the right place...



##########################
## Random vessel effect, #
## generate map internally
##########################

if (mod == 'M0') {
# Parameters with fixed gear estimates and a map
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_xi)
  Obj = TmbList[["Obj"]]
}

  ################################################
  ## When using catchability covariates, these ###
  ## are estimated from another fit ##############
  ################################################

if (mod %in% c('M1', 'M2')) {
  ## loading in values
  load(file.path('../','QEstimates.RData')) ## results folder

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

}

if (mod == 'M3') {

  ## loading in values
 load(file.path('../','QEstimates.RData')) ## results folder - to fix the problem gear-species combinations

 # Define parameter list
  Params <- Param_Fn(Version = Version, DataList = TmbData, RhoConfig = RhoConfig) # - list of parameters

  # extract and assign lambda k values - is this case, use these as the
# starting pos
  lam1_k <- unlist(lapply(Qs, function(x) x$par.fixed[names(x$par.fixed) == 'lambda1_k']))
  lam2_k <- unlist(lapply(Qs, function(x) x$par.fixed[names(x$par.fixed) == 'lambda2_k']))

  # cod juv, bud juv should be fixed at their adult values 
  pos_juv <- which(colnames(Vess_Cov) %in% c('CARLHELMAR_Gadus morhua_Juv','CARLHELMAR_Lophius budegassa_Juv'))
  pos_adu <- which(colnames(Vess_Cov) %in% c('CARLHELMAR_Gadus morhua_Adu','CARLHELMAR_Lophius budegassa_Adu'))
 
  Params$lambda1_k[pos_juv] <- lam1_k[pos_adu]
  Params$lambda2_k[pos_juv] <- lam2_k[pos_adu]

  map <- Make_Map(TmbData = TmbData, TmbParams = Params, CovConfig = TRUE, Q_Config = TRUE,RhoConfig = RhoConfig) # make the map
 
 # Replace the values for carlhelmar cod/bud juv  - only fix these parameters,
  # allow rest to be estimated
  map$lambda1_k[pos_juv] <- NA 
  map$lambda2_k[pos_juv] <- NA

  # Relevel
  levels(map$lambda1_k) <- map$lambda1_k
  levels(map$lambda2_k) <- map$lambda2_k

  # Parameters with fixed gear estimates and a map
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_xi, "Parameters" = Params, "Map" = map)
  Obj = TmbList[["Obj"]]

}



  ############
  # Run model
  ###########
  Opt = TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0('../',DateFile), bias.correct=BiasCorr) 
  Report = Obj$report()

#################################  
######## Save outputs ###########
#################################
Save = list(Obj = Obj,"Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData, "Data_Geostat" = Data_Geostat)
save(Save, file=paste0('../',DateFile,"Save.RData"))


