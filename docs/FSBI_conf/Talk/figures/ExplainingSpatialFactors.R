######################################################################
## Spatial factor loading variance explanation by covariates
######################################################################

###########################
### The factor loadings ###
###########################
library(VAST)

run <- '2017-06-16_M1'
load(file = file.path('..','..', '..', '..','results', run, 'Save.RData'))

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
load(file.path('..','..', '..', '..','results', 'CovariatesAtKnot.RData'))

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
abline(coef(lm(O1$factor_1 ~ log(abs(O1$Depth)))))
summary(lm(O1$factor_1 ~ log(abs(O1$Depth))))
cor(O1$factor_1, log(abs(O1$Depth)))
cor.test(O1$factor_1, log(abs(O1$Depth)))
cor(O1$factor_2, log(abs(O1$Depth)))
cor(O1$factor_3, log(abs(O1$Depth)))

# Relationship with habitat type
library(ggplot2); library(cowplot)
library(ggthemes)
p1 <- ggplot(data = O1, aes(x = log(abs(Depth)), y = factor_1)) + geom_point() + 
	geom_smooth(method = 'lm') + xlab("log(Depth)") + ylab("factor 1 score") +
	ggtitle('Spatial encounter probability adjusted R2 = 0.72') 
p2 <- ggplot(data = O1, aes(x = Habitat, y = factor_1)) + geom_boxplot() + 
	xlab("Substrate type") + ylab("factor 1 score") + theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
	ggtitle('Spatial encounter probability') 

print(p1)
ggsave(file = 'Factor1_DepthO1.png', width = 12, height = 8)

print(p2)
ggsave(file = 'Factor1_HabitatO1.png', width = 12, height = 8)


#plot_grid(p1, p2, ncol = 1)

#summary(aov(O1$factor_1 ~ O1$Habitat))
#coefficients(aov(O1$factor_1 ~ O1$Habitat))

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

# Relationship with habitat type
p1 <- ggplot(data = O2, aes(x = log(abs(Depth)), y = factor_1)) + geom_point() + 
	geom_smooth(method = 'lm') + xlab("log(Depth)") + ylab("factor 1 score") +
	ggtitle("Spatial densit adjusted R2 = 0.51") 
p2 <- ggplot(data = O2, aes(x = Habitat, y = factor_1)) + geom_boxplot() + 
	xlab("Substrate type") + ylab("factor 1 score") + theme(axis.text.x = element_text(angle = -45, hjust = 0)) +
	ggtitle("Spatial density") 

print(p1)
ggsave(file = 'Factor1_DepthO2.png', width = 12, height = 8)

print(p2)
ggsave(file = 'Factor1_HabitatO2.png', width = 12, height = 8)



fit1 <- randomForest(O2$factor_1 ~ log(abs(O2$Depth)) + O2$Habitat)
print(fit1) # view results 
importance(fit1)
fit2 <- randomForest(O2$factor_2 ~ log(abs(O2$Depth)) + O2$Habitat)
print(fit2) # view results 
importance(fit2)
fit3 <- randomForest(O2$factor_3 ~ log(abs(O2$Depth)) + O2$Habitat)
print(fit3) # view results 
importance(fit3)


###########################################################

## Spatio-temporal factor loadings

Epsilon1_sf   <- as.data.frame.table(Var_list$"Epsilon1"$Psi_rot)
Epsilon1_sf$Var1 <- 1:266
Epsilon1_sf$Var2 <- rep(1:9, each = 266)
Epsilon1_sf$Var3 <- rep(1990:2015, each  = 266 * 9)
colnames(Epsilon1_sf) <- c("knot", "factor", "year", "value")

Epsilon1_sf <- Epsilon1_sf[Epsilon1_sf$knot %in% 1:250,]

## Load in the temperature data
load(file.path('..','..', '..', '..','data', 'Covariates' , 'YearlyMeanandCumSumSSTatKnot.RData'))


print(ggplot(Temps, aes(x = Year, y = TempMean)) + geom_line(aes(colour = factor(knot)))  +
	theme(legend.position = "none") + ylab("T degrees C") + xlab(""))
ggsave('Temp.png', width = 12, height = 4)

Epsilon1_sf$TempMean <- Temps$TempMean[match(paste(Epsilon1_sf$knot, Epsilon1_sf$year),
					     paste(Temps$knot, Temps$Year))]

Epsilon1_sf$TempCumSum <- Temps$TempCumSum[match(paste(Epsilon1_sf$knot, Epsilon1_sf$year),
					     paste(Temps$knot, Temps$Year))]

library(ggplot2)

p3 <- ggplot(Epsilon1_sf[Epsilon1_sf$factor == 1,], aes(x = TempMean, y = value)) + geom_point() +
	geom_smooth(method = 'lm') + xlab("Mean temperature") + ylab("factor 1 score")

print(p3)
ggsave(file = 'Factor1_Temp.eps', height = 8, width = 12)


EpF1 <- Epsilon1_sf[Epsilon1_sf$factor == 1,]
EpF1 <- EpF1[!is.na(EpF1$TempMean),]

fit1 <- randomForest(EpF1$value ~ log(EpF1$TempCumSum))
print(fit1) # view results 
importance(fit1)



