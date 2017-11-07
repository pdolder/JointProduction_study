
library(reshape2)
library(ggplot2)

load(file.path('..','results','2017-04-08_M2','Save.RData'))
load(file.path('..','results','CovariatesAtKnot.RData'))

DF <- Save$Data_Geostat
DF$Habitat <- Hab$Habitat[match(DF$knot_i, rownames(Hab))]

Depths <- data.frame(Depth = Depths, knot_i = 1:250)
DF$Depth   <- Depths$Depth[match(DF$knot_i, Depths$knot_i)]

ggplot(DF, aes(x = Habitat, y = log((Catch_KG / AreaSwept_km2)+1))) + geom_boxplot(notch = T) +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))

ggsave(file.path('..', 'plots', 'HabitatCovDiag.png'), width = 16, height = 16)

ggplot(DF, aes(x = abs(Depth), y = log((Catch_KG / AreaSwept_km2)+1))) + geom_point() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + geom_smooth()
ggsave(file.path('..', 'plots', 'DepthCovDiag.png'), width = 16, height = 16)


DF$Depth2 <- DF$Depth^2

m_1 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ spp , data = DF)
m0 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ abs(Depth) , data = DF)
m1 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ abs(Depth) + spp , data = DF)
m2 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ abs(Depth) + abs(Depth2) + spp , data = DF)
m3 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ abs(Depth) + abs(Depth2) + Habitat + spp , data = DF)
m4 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ abs(Depth) + abs(Depth2) + Habitat * spp , data = DF)
m5 <- glm((log(Catch_KG / AreaSwept_km2+1)) ~ abs(Depth) * abs(Depth2) * Habitat * spp , data = DF)

summary(m_1)
summary(m0)
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

AIC(m0,m1,m2,m3,m4,m5)
BIC(m0,m1,m2,m3,m4,m5)


ggplot(DF, aes(x = log(abs(Depth)), y = log((Catch_KG / AreaSwept_km2)+1))) + geom_point() +
	facet_wrap(~spp, scale = 'free_y') + theme(axis.text.x = element_text(angle = -90, hjust = 0))
