rm(list = ls())

load(file.path('..', 'results', 'CovariatesAtKnot.RData'))
load(file.path('..', 'results', '2017-04-08_M2', 'Save.RData'))
load(file.path('..','results','Qs_with_labels.RData'))


MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
Report <- Save$Report
###########################################
## Function for getting the catch per gear
###########################################

inv.logit <-function(x) {
  exp(x) / (1+exp(x))
  
}

GetPrediction <- function(enc_or_catch, spp.n , gear, location, year, area = 0.5, Report, Qs_all) {
  
  n <- spp.n # the number of the species-group
  loc <- location # knot to associate
  gr  <- gear # name of gear, where we look up the lambda_Ks
  yr <- year # year as numeric
  a  <- area # area for tow, use a standardised 0.5km
  
  # Encounter probability
  if(enc_or_catch == 'encounter') {
    zeta1 <- Qs_all[Qs_all$Survey == gr & Qs_all$spp == names(n) & Qs_all$Param == 'lambda1_k','Estimate']
    P <- Report$Omega1_sc[loc,n] + Report$eta1_x[loc] + zeta1 + Report$beta1_ct[n, yr] + Report$Epsilon1_sct[loc,n,yr]
  }
   
  # Positive catch rates
  if(enc_or_catch == 'catch') {
    zeta2 <- Qs_all[Qs_all$Survey == gr & Qs_all$spp == names(n) & Qs_all$Param == 'lambda2_k','Estimate']
    P <- Report$Omega2_sc[loc,n] + Report$eta2_x[loc] + zeta2 + Report$beta2_ct[n, yr] + Report$Epsilon2_sct[loc,n,yr]
  } 
    
    
    if(enc_or_catch == 'encounter') {
            R <- inv.logit(P)
    }
    if(enc_or_catch == 'catch') {
      R <- a * exp(P)
    }
    
    return(R)
      
  }


spp <- c('cod_adu' = 1, 'cod_juv' = 2,"meg_adu" = 3, "meg_juv" = 4,"bud_adu" = 5,"bud_juv" = 6,"pisc_adu" = 7,"pisc_juv" = 8,
         "had_adu" = 9,"had_juv" = 10, "whg_adu" = 11,"whg_juv" = 12,"hke_adu" = 13,"hke_juv" = 14,"ple_adu" = 15,
         "ple_juv" = 16,"sol_adu" = 17,"sol_juv" = 18)


yr <- 1:length(1990:2015)


# Now run for locations 66 and 79
PredDF <- expand.grid(location = 1:250, spp = c("cod_adu", "had_adu"),
            gear = c('THA2','NWGFS'), year = 22:26, EncProb = NA, PosCatch = NA, PredCatchKg = NA)


# Look through locations, species and gears
for (loc in unique(PredDF$location)) {
    for (sp in unique(PredDF$spp)) {
      for (gr in unique(PredDF$gear)) {
        for (y in unique(PredDF$year)) {
        
        # Predicted encounter prob
        PredDF$EncProb[PredDF$location == loc & PredDF$spp == sp & PredDF$gear == gr & PredDF$year == y] <- 
          GetPrediction(enc_or_catch = 'encounter', spp.n = spp[sp], gear = gr, location = loc,
                        year = y, Report = Report, Qs_all = Qs_all)
        
        # Predicted catch rate
        PredDF$PosCatch[PredDF$location == loc & PredDF$spp == sp & PredDF$gear == gr & PredDF$year == y] <- 
          GetPrediction(enc_or_catch = 'catch', spp.n = spp[sp], gear = gr, location = loc,
                        year = y, Report = Report, Qs_all = Qs_all)
      
          }
            }

        }
}

PredDF$PredCatchKg <- PredDF$EncProb * PredDF$PosCatch

# Want to plot as scatter plot, spp1 v spp2
library(reshape2); library(dplyr)

TechPlotDF <- dcast(PredDF, location + gear + year ~ spp, value.var = "PredCatchKg")

## Getting the convex hull
# https://stats.stackexchange.com/questions/22805/how-to-draw-neat-polygons-around-scatterplot-regions-in-ggplot2
library(plyr)
find_hull <- function(TechPlotDF) TechPlotDF[chull(TechPlotDF$cod_adu, TechPlotDF$had_adu), ]
hulls <- ddply(TechPlotDF, c("year","gear"), find_hull)

library(ggplot2)

ggplot(filter(TechPlotDF, gear == 'THA2'), aes(x = cod_adu, y = had_adu)) + geom_point() +
	facet_wrap(~year) + geom_polygon(data = filter(hulls,gear=="THA2"), alpha = 0.2) +
	theme_classic()

ggplot(filter(TechPlotDF, gear == 'NWGFS'), aes(x = cod_adu, y = had_adu)) + geom_point() +
	facet_wrap(~year) + geom_polygon(data = filter(hulls,gear=="NWGFS"), alpha = 0.2) +
	theme_classic()

ggplot(TechPlotDF, aes(x = cod_adu, y = had_adu)) + geom_point(aes(colour = factor(year))) +
	facet_wrap(~gear) + geom_polygon(data = hulls, alpha = 0.2, aes(colour = factor(year), fill = factor(year))) +
	theme_classic()

