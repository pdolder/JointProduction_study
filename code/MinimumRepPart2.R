
library(TMB)
library(VAST)

load(file = 'MinExampleHess.RData')

  TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "OverdispersionConfig" = OverdispersionConfig, "c_i"=as.numeric(Data_Geostat[,'spp'])-1, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "Q_ik" = Vess_Cov,"s_i"=Data_Geostat[,'knot_i']-1, "t_iz"=as.numeric(Data_Geostat[,'Year']), "a_xl"=Spatial_List$a_xl, "MeshList"=Spatial_List$MeshList, "GridList"=Spatial_List$GridList, "Method"=Spatial_List$Method)

  # Make TMB object
#  dyn.unload( paste0(DateFile,"/",dynlib(TMB:::getUserDLL())))
  setwd(DateFile) # so executables go to the right place...

  # Parameters with fixed gear estimates and a map
  TmbList = Build_TMB_Fn("TmbData"=TmbData, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_xi)
  Obj = TmbList[["Obj"]]

  # Run model
  Opt = TMBhelper::Optimize(obj=Obj, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, savedir=paste0('../',DateFile), bias.correct=BiasCorr) 
  Report = Obj$report()


  Opt$SD  
  
  ## hessian not positive definitive, yet if I run the same model on each
  ## species individually (estimating gear-specific catchability covariates) it
  ## is fine


