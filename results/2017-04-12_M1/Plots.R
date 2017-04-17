
# Load libraries
library(TMB)
library(ThorsonUtilities)
library(VAST)

# This is where all runs will be located
DateFile <- file.path(paste(getwd(),'/',sep=''))

load('Save.RData')
load(file.path('..', 'CovariatesAtKnot.RData'))
###############
# Settings
###############

################
# Make diagnostic plots
################
an <- as.numeric
strata.limits <- data.frame('STRATA'="All_areas")
BiasCorr <- FALSE

DF <- Save$Data_Geostat

  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))

spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')



###################


  ## Plot encounter probability 
  Enc_prob <- SpatialDeltaGLMM::Check_encounter_prob(Report = Save$Report, Data_Geostat = Save$Data_Geostat, DirName = DateFile)

  # Plot Anisotropy  
  SpatialDeltaGLMM::PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Save$Report, TmbData=Save$TmbData )

  # Plot covariances
  Cov_List = Summarize_Covariance( Report=Save$Report, ParHat=Save$ParHat, Data=Save$TmbData, SD=Save$Opt$SD, plot_cor=TRUE, category_names=spp, figname=paste0("Spatio-temporal_covariances"), plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )

  # Plot overdispersion
#  Plot_Overdispersion( filename1=paste0(DateFile,"Overdispersion"), filename2=paste0(DateFile,"Overdispersion--panel"), Data=Save$TmbData, ParHat=Save$ParHat, Report=Save$Report, ControlList1=list("Width"=5, "Height"=10, "Res"=600, "Units"='in'), ControlList2=list("Width"=TmbData$n_c, "Height"=TmbData$n_c, "Res"=200, "Units"='in'))

  # Plot index
  SpatialDeltaGLMM::PlotIndex_Fn( DirName=DateFile, TmbData=Save$TmbData, Sdreport=Save$Opt$SD, Year_Set=seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year'])))), strata_names=strata.limits[,1], category_names=spp, use_biascorr=BiasCorr, cex = 0.3)


 # Positive catch rate Q-Q plot
  Q = SpatialDeltaGLMM::QQ_Fn( TmbData=Save$TmbData, Report=Save$Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

 # Plot surface - this is for all spp
  Dim = c( "Nrow"=ceiling(sqrt(length(Years2Include))), "Ncol"=ceiling(length(Years2Include)/ceiling(sqrt(length(Years2Include)))))
  MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
SpatialDeltaGLMM::PlotResultsOnMap_Fn(plot_set=3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Save$Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateFile,"Field_"), Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], category_names=levels(DF[,'spp']), mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]], cex=1.8, Legend=MapDetails_List[["Legend"]])


## Plot factors

Plot_factors(Report = Save$Report, ParHat = Save$ParHat, Data = Save$TmbData, SD = Save$Opt$SD, mapdetails_list = MapDetails_List, Year_Set = Year_Set, 
	     category_names = spp, plotdir = DateFile)

PCAstyle <- TRUE

if(PCAstyle == TRUE) {
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

## Plots
library(ggplot2)
library(ggthemes)
library(cowplot)

dem <- paste(rep(c('cod', 'had', 'whg'), each = 2),rep(c('adu','juv'),times = 3), sep =  '_')
fla <- paste(rep(c('ple', 'sol'), each = 2),rep(c('adu','juv'),times = 2), sep =  '_')
dep <- paste(rep(c('hke', 'bud', 'pis', 'meg'), each = 2),rep(c('adu','juv'),times = 4), sep =  '_')


PCA.DFs <- lapply(Var_list, function(x) {
			  DF <- as.data.frame(x$L_pj_rot)
			  colnames(DF) <- paste('Factor',1:9, sep = '_')
			  DF$group <- ifelse(rownames(DF) %in% dem, 'roundfish', ifelse(rownames(DF) %in% fla, 'flatfish', 'deep'))
			  return(DF)
	     })

grps <- data.frame(Factor_1 = rep(-1.1,3), Factor_2 = c(-0.9, -0.95, -1), group = c('roundfish', 'flatfish', 'deep'))

p1 <- ggplot(PCA.DFs[['Epsilon1']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
	expand_limits(y = c(-1,0.5), x = c(-1,0.5)) + 	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw() + 
	theme(legend.position = 'none', axis.text = element_text(face = 'bold')) + xlab('Factor 1 (39.8 % variance explained)') + ylab('Factor 2 (16.7 % variance explained)') +
	geom_rect(aes(xmin  = -1.25, xmax = -0.95, ymin = -1.05, ymax = -0.85), fill= 'white', colour = 'grey') +
	geom_label(data = grps, aes(label = group, fill = group), alpha = 0.2) 

p2 <- ggplot(PCA.DFs[['Epsilon2']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
	expand_limits(y = c(-1,0.5), x = c(-1,0.5)) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw()  +
	theme(legend.position = 'none', axis.text = element_text(face = 'bold'))+ xlab('Factor 1 (43.3 % variance explained)') + ylab('Factor 2 (13.9 % variance explained)')


plot_grid(p1,p2, labels = c('(a) Factor loading for spatio-temporal \n        encounter probability',
			    '(b) Factor loading for spatio-temporal \n            catch rates'), hjust = -0.2, vjust = 1.5)
ggsave(file = 'PCAstyle_Plots_Spatiotemporal.png', width = 16, height = 8)

## Spatial 

p3 <- ggplot(PCA.DFs[['Omega1']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
	expand_limits(y = c(-2.6,0.1), x = c(-2.5,3.5)) +	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw() + 
	theme(legend.position = 'none', axis.text = element_text(face = 'bold')) + xlab('Factor 1 (41 % variance explained)') + ylab('Factor 2 (31.3 % variance explained)') +
	geom_rect(aes(xmin  = -1.25, xmax = -0.95, ymin = -1.05, ymax = -0.85), fill= 'white', colour = 'grey') +
	geom_label(data = grps, aes(label = group, fill = group), alpha = 0.2) 

p4 <- ggplot(PCA.DFs[['Omega2']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw()  +
	theme(legend.position = 'none', axis.text = element_text(face = 'bold'))+ xlab('Factor 1 (35.6 % variance explained)') + ylab('Factor 2 (18.1 % variance explained)')

plot_grid(p3,p4, labels = c('(a) Factor loading for spatial \n        encounter probability',
			    '(b) Factor loading for spatial \n            catch rates'), hjust = -0.2, vjust = 1.5)
ggsave(file = 'PCAstyle_Plots_Spatial.png', width = 16, height = 8)

}


#SpatialDeltaGLMM::Plot_data_and_knots(Data_Extrap = Spatial_List$NN_Extrap, 
#				      Extrap_Area_km2 = Spatial_List$a_xl, 
#				      loc_x = Spatial_List$loc_x,
#				      Spatial_List$loc_x, 
#				      Data_Geostat = Data_Geostat,
#				      Plot_name = 'Data_and_knots.png',
#				      Data_name = 'Data_by_year.png')

#### END ###

