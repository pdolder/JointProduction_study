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
PredDF <- expand.grid(location = 1:250, spp = c("cod_adu", "had_adu", "whg_adu", "ple_adu", "sol_adu"),
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

## First pass guestimate of raising factor
# maximum biomass in global index
# maximum landings in STECF data

Raise <- TRUE
if(Raise == T) {
RaisingFactor <- data.frame(species = c('cod_adu','had_adu','whg_adu','ple_adu','sol_adu','hke_adu', 'bud_adu', 'pisc_adu', 'meg_adu'), 
	   MaxB = c(2000, 21000, 20000, 5000, 400, 10000, 0.4, 2, 1500),
	   MaxL = c(7000, 17000, 13000, 2000, 2300, 43000, 7700, 25500, 15000))

# Assume biomass 2 X Landings
RaisingFactor$RF <- (RaisingFactor$MaxL * 2) / RaisingFactor$MaxB

PredDF$RF <- RaisingFactor$RF[match(PredDF$spp, RaisingFactor$species)]
PredDF$PredCatchKg <- PredDF$PredCatchKg * PredDF$RF
}

# Want to plot as scatter plot, spp1 v spp2
library(reshape2); library(dplyr)

TechPlotDF <- dcast(PredDF, location + gear + year ~ spp, value.var = "PredCatchKg")

TechPlotDF$year[TechPlotDF$year==22] <- 2011
TechPlotDF$year[TechPlotDF$year==23] <- 2012
TechPlotDF$year[TechPlotDF$year==24] <- 2013
TechPlotDF$year[TechPlotDF$year==25] <- 2014
TechPlotDF$year[TechPlotDF$year==26] <- 2015
TechPlotDF$year  <- as.factor(TechPlotDF$year)
## Getting the convex hull
# https://stats.stackexchange.com/questions/22805/how-to-draw-neat-polygons-around-scatterplot-regions-in-ggplot2
library(plyr)
library(ggplot2)

## Cod:Haddock
find_hull <- function(TechPlotDF) TechPlotDF[chull(TechPlotDF$cod_adu, TechPlotDF$had_adu), ]
hulls <- ddply(TechPlotDF, c("year","gear"), find_hull)

TechPlot_Otter <- filter(TechPlotDF, gear == 'THA2')

p1 <- ggplot(TechPlot_Otter, aes(x = cod_adu, y = had_adu)) + geom_point(aes(colour = year)) +
	geom_polygon(data = filter(hulls,gear=="THA2"), aes(colour = year, fill = year), alpha = 0.1) +
	theme_bw() + expand_limits(x = c(0,max(c(TechPlot_Otter$cod_adu, TechPlot_Otter$had_adu))), 
				   y = c(0,max(c(TechPlot_Otter$cod_adu, TechPlot_Otter$had_adu)))) +
xlab("cod") + ylab("haddock")


find_hull <- function(TechPlotDF) TechPlotDF[chull(TechPlotDF$ple_adu, TechPlotDF$sol_adu), ]
hulls <- ddply(TechPlotDF, c("year","gear"), find_hull)


TechPlot_Beam <- filter(TechPlotDF, gear == 'NWGFS')
p2 <- ggplot(TechPlot_Beam, aes(x = ple_adu, y = sol_adu)) + geom_point(aes(colour = year)) +
	geom_polygon(data = filter(hulls,gear=="NWGFS"), aes(colour = year, fill = year), alpha = 0.1) +
	theme_bw() + expand_limits(x = c(0,max(c(TechPlot_Beam$ple_adu, TechPlot_Beam$sol_adu))), 
				   y = c(0,max(c(TechPlot_Beam$ple_adu, TechPlot_Beam$sol_adu)))) +
xlab("plaice") + ylab("sole")


library(cowplot)
plot_grid(p1,p2, labels = c("(a) cod v haddock using otter trawl",
			    "(b) plaice v sole using beam trawl"), hjust = -0.25, vjust = 2)

ggsave(file = file.path('..','results','2017-04-08_M2','TechnicalEfficiency.png'), width = 12, height = 6)


