#####################################################################
### Code for plotting differences in density for cod, had and whg ###
#####################################################################

rm(list = ls())

run <- '2017-06-16_M1'

load(file.path('..', '..', '..','results', 'CovariatesAtKnot.RData'))
load(file.path('..', '..' , '..','results', run, 'Save.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

cod <- cod_or <-  Save$Report$Index_xcyl[,1,,]
had <- had_or <- Save$Report$Index_xcyl[,9,,]
whg <- whg_or <- Save$Report$Index_xcyl[,11,,]
ple <- ple_or <- Save$Report$Index_xcyl[,15,,]
sol <- sol_or <- Save$Report$Index_xcyl[,17,,]
hke <- hke_or <- Save$Report$Index_xcyl[,13,,]
pis <- pis_or <- Save$Report$Index_xcyl[,7,,]
meg <- meg_or <- Save$Report$Index_xcyl[,3,,]

# Check overall index
#par(mfrow = c(1,3))
#plot(colSums(cod), type = 'b')
#plot(colSums(had), type = 'b') 
#plot(colSums(whg), type = 'b')

# Now standardise the 2015 values
cod <- cod[,26] / sum(cod[,26]) * 100
had <- had[,26] / sum(had[,26]) * 100 
whg <- whg[,26] / sum(whg[,26]) * 100
ple <- ple[,26] / sum(ple[,26]) * 100 
sol <- sol[,26] / sum(sol[,26]) * 100
hke <- hke[,26] / sum(hke[,26]) * 100
pis <- pis[,26] / sum(pis[,26]) * 100
meg <- meg[,26] / sum(meg[,26]) * 100

# Make the plot dataframe

plotDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers))

DF <- data.frame(X = Centers[,'E_km'], Y = Centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29 

LLs <- PBSmapping::convUL(DF)

plotDF$Lon <- LLs$Y
plotDF$Lat <- LLs$X

## Add the densities
plotDF$cod <- cod
plotDF$had <- had
plotDF$whg <- whg
plotDF$ple <- ple 
plotDF$sol <- sol 
plotDF$hke <- hke
plotDF$pis <- pis
plotDF$meg <- meg

## Area densities relative to... 
plotDF$cod_had  <- plotDF$cod - plotDF$had
plotDF$cod_whg  <- plotDF$cod - plotDF$whg
plotDF$had_whg  <- plotDF$had - plotDF$whg

plotDF$ple_sol  <- plotDF$ple - plotDF$sol

plotDF$hke_pis  <- plotDF$hke - plotDF$pis
plotDF$hke_meg  <- plotDF$hke - plotDF$meg
plotDF$pis_meg  <- plotDF$pis - plotDF$meg

#############
## the Map ##
#############

library(ggplot2); library(cowplot)
library(mapplots)
library(rworldmap)
library(rworldxtra)
library(data.table)
library(broom)
library(RColorBrewer)

mapdata =ggplot2::fortify(rworldmap::getMap(resolution = 'high'))

xlim = c(-12, -2)
ylim = c(48, 52)
 # if you want to fast subset by lat/long you can do the following:
 # subset to just regions in xlim and ylim, see
# http://stackoverflow.com/a/16574176
 mapdata = data.table(mapdata)
  mapdata = mapdata[mapdata[,.I[any(long %between% xlim) & any(lat %between%
 ylim)], by = list(group)]$V1]

  # make geoms
  coast.poly = geom_polygon(data=mapdata, aes(x=long, y=lat, group=group), colour="#999999", fill="#999999", lwd=0.2)
  coast.outline = geom_path(data=mapdata, aes(x=long, y=lat, group=group), colour="#000000", lwd=0.2)


## The map for all areas

DF2 <- MapDetails_List[['PlotDF']]

DF2$cod_had <- plotDF$cod_had[match(DF2$x2i, plotDF$knot)]
DF2$cod_whg <- plotDF$cod_whg[match(DF2$x2i, plotDF$knot)]
DF2$had_whg <- plotDF$had_whg[match(DF2$x2i, plotDF$knot)]

DF2$ple_sol <- plotDF$ple_sol[match(DF2$x2i, plotDF$knot)]

DF2$hke_pis <- plotDF$hke_pis[match(DF2$x2i, plotDF$knot)]
DF2$hke_meg <- plotDF$hke_meg[match(DF2$x2i, plotDF$knot)]
DF2$pis_meg <- plotDF$pis_meg[match(DF2$x2i, plotDF$knot)]



##########################
### Catch compositions ###
##########################

load(file.path('..', '..', '..','results','Qs_with_labels.RData'))


##
library(dplyr)


#####################################################################################################
##### Caclculating predicted catch for a gear, year, cell and species - Jim's email of 15/5/2017 ####
## Rather than grabbing the predictions, which may not exist for a particular gear at a location, ###
## we need to do longhand
#####################################################################################################
##### apply over location
ggplot() + coast.poly + coast.outline + coord_quickmap(xlim, ylim) + theme(legend.position = 'bottom', legend.text = element_text(angle = -90), legend.title = element_blank())   +
  geom_label(data = LLs, aes(x = X, y = Y, label = 1:250), colour = 'black', size = 2)

## Choose locations

Report <- Save$Report

Data_Geostat[(Data_Geostat$spp == 'Gadus morhua_Adu' & 
        Data_Geostat$Year == 2015 & 
         Data_Geostat$knot_i %in% 24       ),]

i_Ref <- which(Data_Geostat$spp == 'Gadus morhua_Adu' & 
                         Data_Geostat$Year == 2015 & 
                         Data_Geostat$Vessel == 'THA2' &
                         Data_Geostat$knot_i == 24)

Report$R1_i[i_Ref] 
Report$R2_i[i_Ref]

Report$R1_i[i_Ref] * Report$R2_i[i_Ref]

## And for an 'out of bag' prediction
# Linear predictor for a single observation is:

# First part
#P1_i = Omega1_sc(s,c) + eta_x(s) + zeta1_i(i) + eta1_vc 

# Omega1_sc(s,c) = Spatial intercept, location s, species c
# eta_x(s) = spatial variation
# zeta1_i = gear effect [repeats same value], i.e. is lambda1_k for all observations of that gear/species
# eta_vc is ZERO - this is where there is some structure (e.g. AR)

# Second part
#P1_i + beta1_ct + Epsilon1_sct * exp(log_sigmaratio1_z) + eta1_xct

# beta1_ct(s,t) is spatio-temporal variation, dim = spp, year
# Epsilon1_sct is location (266 values ??), species, time
# exp(log_sigmaratio1_z) = log ratio of variance of t_iz (0 ??)
# eta1_xct(s,c,t) -habitat covariates 250, 18, 26 - ZEROS

## Example: first data point is hake (spp 13), THA2 (gear 5) at knot 154 in 1997 (year 8)

P1 <- Report$Omega1_sc[154,13] + Report$eta1_x[154] + Report$zeta1_i[1] + 0
P1 <- P1 + Report$beta1_ct[13,8] + Report$Epsilon1_sct[154,13,8] + exp(log(0)) + 0


# Third part: inverse logit =
inv.logit <-function(x) {
  exp(x) / (1+exp(x))
  
}

print(R1 <- inv.logit(P1))
identical(R1, Report$R1_i[1])


## Lets convert this to a function which works both on the encounter prob and positive catch rates

spp <- c('cod_adu' = 1, 'cod_juv' = 2,"meg_adu" = 3, "meg_juv" = 4,"bud_adu" = 5,"bud_juv" = 6,"pisc_adu" = 7,"pisc_juv" = 8,
         "had_adu" = 9,"had_juv" = 10, "whg_adu" = 11,"whg_juv" = 12,"hke_adu" = 13,"hke_juv" = 14,"ple_adu" = 15,
         "ple_juv" = 16,"sol_adu" = 17,"sol_juv" = 18)


yr <- 1:length(1990:2015)


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



# test the function
GetPrediction(enc_or_catch = 'catch', spp.n = spp["hke_adu"], gear = 'THA2', location = 154, 
              year = 26, Report = Report, Qs_all = Qs_all)


# Now run for locations 66 and 79
PredDF <- expand.grid(location = c(66,79,216), spp = c('cod_adu','had_adu','whg_adu',
                                                   'ple_adu','sol_adu','hke_adu','pisc_adu','meg_adu'),
            gear = c('THA2','NWGFS'), year = 26, EncProb = NA, PosCatch = NA, PredCatchKg = NA)


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

# As a catch composition
library(dplyr)
TotCatch <- PredDF %>% 
  group_by(location, gear, year) %>% 
  summarise(TotalCatch = sum(PredCatchKg))

PredDF$TotalCatch <- TotCatch$TotalCatch[match(paste(PredDF$location, PredDF$gear, PredDF$year),
                                               paste(TotCatch$location, TotCatch$gear, TotCatch$year))]

PredDF$perc <- PredDF$PredCatchKg / PredDF$TotalCatch * 100


###################
## ~~~ Plot ~~~ ###
###################

cols <- brewer.pal(8, 'Set1')
lim <- range(DF2[,5:11])
## The plots
p1 <-  ggplot() + geom_point(data = filter(DF2, Include == T), aes(x = Lon, y =  Lat, colour = cod_had), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) + 
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal", plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[2], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") 


p2 <- ggplot() + geom_point(data = filter(DF2, Include == T), aes(x = Lon, y = Lat, colour = cod_whg), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal",plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[3], mid = 'white', high = cols[1], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")#, limits = lim)  

p3 <- ggplot() + geom_point(data = filter(DF2, Include ==T), aes(x = Lon, y = Lat, colour = had_whg), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal",plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[3], mid = 'white', high = cols[2], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")#, limits = lim)

p4 <- ggplot() + geom_point(data = filter(DF2, Include ==T), aes(x = Lon, y = Lat, colour = hke_pis), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal",plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[7], mid = 'white', high = cols[6], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")#, limits = lim)

p5 <- ggplot() + geom_point(data = filter(DF2, Include ==T), aes(x = Lon, y = Lat, colour = hke_meg), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
    coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal",plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[8], mid = 'white', high = cols[6], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")#, limits = lim)

p6 <- ggplot() + geom_point(data = filter(DF2, Include ==T), aes(x = Lon, y = Lat, colour = pis_meg), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
  coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal",plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[8], mid = 'white', high = cols[7], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar")#, limits = lim)


p7 <- ggplot() + geom_point(data = filter(DF2, Include ==T), aes(x = Lon, y = Lat, colour = ple_sol), size = 1) + geom_point(data = DF2, aes(x = Lat, y = Lon), size = 0.5, shape = 3) +
  coast.poly + coast.outline  + coord_quickmap(xlim, ylim) + theme(legend.position = c(0.1,1), legend.justification = c(0,0), legend.title = element_blank(), legend.text = element_text(size = 8), legend.direction = "horizontal",plot.margin=unit(c(0,0,0,0),"mm")) + xlab('') + ylab('') +
  scale_colour_gradient2(low = cols[5], mid = 'white', high = cols[4], midpoint = 0, space = "Lab", na.value = "grey50", guide = "colourbar") + #, limits = lim) +
  geom_text(data = LLs[c(66,79,216),], aes(x = X, y = Y), label = c('A','B','C'), size = 3) + geom_point(data = LLs[c(66,79, 216),], aes(x = X, y = Y), shape = 'o', colour = 'red', size = 14)

## Change the gear names

levels(PredDF$gear)[levels(PredDF$gear) == "THA2"]  <- "Otter"
levels(PredDF$gear)[levels(PredDF$gear) == "NWGFS"] <- "Beam"

## convert the catch composition locations

PredDF$location[PredDF$location == 66]  <- 'A'
PredDF$location[PredDF$location == 79]  <- 'B'
PredDF$location[PredDF$location == 216] <- 'C'

p8  <- ggplot(filter(PredDF,gear %in% c('Otter','Beam')), 
            aes(x = factor(gear), y = perc)) + geom_bar(stat= 'identity',aes(fill = spp), colour = 'black') + facet_wrap(~location) +
  scale_fill_manual(values = cols, labels = c('cod','haddock','whiting','plaice','sole','hake','anglerfishes','megrim')) + xlab('') + ylab('') + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.title = element_blank(),plot.margin=unit(c(0,0,0,0),"mm"))


comb_plot <- ggdraw() + draw_plot(p1, x = 0, y = 0.6, width = 0.3, height = 0.3) + # cod:haddock
  draw_plot(p2, x = 0, y = 0.3, width = 0.3, height = 0.3) + # cod:whiting
  draw_plot(p3, x = 0, y =  0, width = 0.3, height = 0.3) + # haddock:whiting
  draw_plot(p4, x = 0.3, y = 0.6, width = 0.3, height = 0.3)  + # hake:monk
  draw_plot(p5, x = 0.3, y = 0.3, width = 0.3, height = 0.3) + # hake:meg
  draw_plot(p6, x = 0.3, y = 0, width = 0.3, height = 0.3) + # monk:meg
  draw_plot(p7, x = 0.6, y = 0.6, width = 0.3, height = 0.3) + # # plaice:sole
  
  draw_plot(p8, x = 0.65, y = 0, width = 0.3, height = 0.55) + # catch composition
  draw_plot_label(c("(A) haddock : cod", "(B) whiting : cod", "(C) whiting : haddock", '(D) anglerfishes : hake','(E) megrim : hake',
                    '(F) megrim : anglerfishes','(G) sole : plaice','(H) Catch composition'), 
                  x = c(0.05, 0.05, 0.05, 0.35, 0.35, 0.35, 0.65, 0.65), y = c(0.92, 0.62, 0.32,0.92,0.62,0.32, 0.92, 0.62), 
                  size = 15, hjust = 0)

save_plot(file.path('..','figures','Figure 4 - DensityDifferencesFigureswithCC.png'), comb_plot,ncol = 3, nrow = 3)



