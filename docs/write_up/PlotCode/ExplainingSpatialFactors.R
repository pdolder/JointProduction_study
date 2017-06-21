######################################################################
## Spatial factor loading variance explanation by covariates
######################################################################

###########################
### The factor loadings ###
###########################
library(VAST)

run <- '2017-06-16_M1'
load(file = file.path('..', '..', '..','results', run, 'Save.RData'))

an <- as.numeric
DF <- Save$Data_Geostat

  Year_Set = seq(min(an(as.character(DF[,'Year']))),max(an(as.character(DF[,'Year']))))
  Years2Include = which(Year_Set %in% sort(unique(an(as.character(DF[,'Year'])))))

spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')


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

Omega1_sf   <- apply(Var_list$"Omega1"$Psi_rot, 1:2, FUN = mean)
Omega2_sf   <- apply(Var_list$"Omega2"$Psi_rot, 1:2, FUN = mean)

## The map for all areas
Omega1_sf <- Omega1_sf[1:250,]
Omega2_sf <- Omega2_sf[1:250,]

Omega1_sf <- as.data.frame(Omega1_sf)
Omega2_sf <- as.data.frame(Omega2_sf)

Omega1_sf$knot <- 1:250
Omega2_sf$knot <- 1:250

colnames(Omega1_sf)[1:9] <- colnames(Omega2_sf)[1:9] <- paste("factor",1:9, sep = '_')

######################
### The covariates ###
######################

## Read in Spatial List
load(file.path('..', '..', '..','results', 'CovariatesAtKnot.RData'))

MapDetails_List = SpatialDeltaGLMM::MapDetails_Fn( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )

HabDF <- data.frame(x = Spatial_List$Kmeans$centers[,1], 
		    y= Spatial_List$Kmeans$centers[,2],
		    knot = 1:nrow(Spatial_List$Kmeans$centers),
		    Habitat = Hab$Habitat, Depth = Depths)

DF <- data.frame(X = Centers[,'E_km'], Y = Centers[,'N_km'])
attr(DF, 'projection') = 'UTM'
attr(DF, "zone") <- 29 

LLs <- PBSmapping::convUL(DF)

HabDF$Lon <- LLs$Y
HabDF$Lat <- LLs$X

################
### Combined ###
################

## Covariates and Omega1
O1 <- merge(HabDF, Omega1_sf)
O2 <- merge(HabDF, Omega2_sf)


##############################
## Relationship with log Depth
##############################
plot(O1$factor_1 ~ log(abs(O1$Depth)))
summary(lm(O1$factor_1 ~ log(abs(O1$Depth))))
cor(O1$factor_1, log(abs(O1$Depth)))
cor.test(O1$factor_1, log(abs(O1$Depth)))
cor(O1$factor_2, log(abs(O1$Depth)))
cor(O1$factor_3, log(abs(O1$Depth)))

# Relationship with habitat type
library(ggplot2); library(cowplot)
p1 <- ggplot(data = O1, aes(x = Habitat, y = factor_1)) + geom_boxplot()
p2 <- ggplot(data = O1, aes(x = Habitat, y = log(abs(Depth)))) + geom_boxplot()

plot_grid(p1, p2, ncol = 1)

summary(aov(O1$factor_1 ~ O1$Habitat))
coefficients(aov(O1$factor_1 ~ O1$Habitat))

## Hows about a Random Forest classification tree to look at covariate
## contribution to variance

library(randomForest)
fit1 <- randomForest(O1$factor_1 ~ log(abs(O1$Depth)) + O1$Habitat)
print(fit1) # view results 
importance(fit1)
fit2 <- randomForest(O1$factor_2 ~ log(abs(O1$Depth)) + O1$Habitat)
print(fit2) # view results 
importance(fit2)
fit3 <- randomForest(O1$factor_3 ~ log(abs(O1$Depth)) + O1$Habitat)
print(fit3) # view results 
importance(fit3)

# Covariates and Omega2
plot(O2$factor_1 ~ log(abs(O2$Depth)))
summary(lm(O2$factor_1 ~ log(abs(O2$Depth))))
cor(O2$factor_1, log(abs(O2$Depth)))
cor.test(O2$factor_1, log(abs(O2$Depth)))
cor(O2$factor_2, log(abs(O2$Depth)))
cor(O2$factor_3, log(abs(O2$Depth)))

fit1 <- randomForest(O2$factor_1 ~ log(abs(O2$Depth)) + O2$Habitat)
print(fit1) # view results 
importance(fit1)
fit2 <- randomForest(O2$factor_2 ~ log(abs(O2$Depth)) + O2$Habitat)
print(fit2) # view results 
importance(fit2)
fit3 <- randomForest(O2$factor_3 ~ log(abs(O2$Depth)) + O2$Habitat)
print(fit3) # view results 
importance(fit3)



