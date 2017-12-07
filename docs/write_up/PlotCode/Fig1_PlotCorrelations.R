######################################################
### | Code to Plot the inter-species correlation | ###
######################################################
library(VAST)

run <- '2017-06-16_M1'

load(file.path('..', '..', '..', 'results', run, 'Save.RData'))


## Extract parameter estimates 
sdreport <- summary(Save$Opt$SD)

spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')

Cov_List = Summarize_Covariance( Report=Save$Report, ParHat=Save$ParHat, Data=Save$TmbData, SD=Save$Opt$SD, plot_cor=TRUE, category_names=spp, figname=paste0("Spatio-temporal_covariances"), plotTF=c("Omega1"=TRUE,"Epsilon1"=TRUE,"Omega2"=TRUE,"Epsilon2"=TRUE), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,5,2,2) )


COR_O1  <- Cov_List[["Cor_omega1"]][,,"Estimate"] 
SE_O1   <- Cov_List[["Cor_omega1"]][,,"Std.Error"] ; SE_O1[is.na(SE_O1)] <- 0
Up_O1   <- COR_O1 + 1.96 * SE_O1
Lo_O1   <- COR_O1 - 1.96 * SE_O1
# a psuedo-pvalue matrix for corrplot
P_O1 <- Lo_O1 < 0 & Up_O1 > 0
P_O1 <- ifelse(P_O1 == TRUE, 0.5, 0.01)


COR_O2  <- Cov_List[["Cor_omega2"]][,,"Estimate"] 
SE_O2   <- Cov_List[["Cor_omega2"]][,,"Std.Error"] ; SE_O2[is.na(SE_O2)] <- 0
Up_O2   <- COR_O2 + 1.96 * SE_O2
Lo_O2   <- COR_O2 - 1.96 * SE_O2
# a psuedo-pvalue matrix for corrplot
P_O2 <- Lo_O2 < 0 & Up_O2 > 0
P_O2 <- ifelse(P_O2 == TRUE, 0.5, 0.01)

COR_E1  <-Cov_List[["Cor_epsilon1"]][,,"Estimate"] 
SE_E1   <- Cov_List[["Cor_epsilon1"]][,,"Std.Error"] ; SE_E1[is.na(SE_E1)] <- 0
Up_E1   <- COR_E1 + 1.96 * SE_E1
Lo_E1   <- COR_E1 - 1.96 * SE_E1
# a psuedo-pvalue matrix for corrplot
P_E1 <- Lo_E1 < 0 & Up_E1 > 0
P_E1 <- ifelse(P_E1 == TRUE, 0.5, 0.01)

COR_E2  <-Cov_List[["Cor_epsilon2"]][,,"Estimate"] 
SE_E2   <- Cov_List[["Cor_epsilon2"]][,,"Std.Error"] ; SE_E2[is.na(SE_E2)] <- 0
Up_E2   <- COR_E2 + 1.96 * SE_E2
Lo_E2   <- COR_E2 - 1.96 * SE_E2
# a psuedo-pvalue matrix for corrplot
P_E2 <- Lo_E2 < 0 & Up_E2 > 0
P_E2 <- ifelse(P_E2 == TRUE, 0.5, 0.01)

library(corrplot)
library(RColorBrewer)

spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')

# Omega 1
png(file = file.path("..", 'plots', "Omega1_Correlations.png"), width = 800, height = 800)
corrplot(COR_O1,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black', p.mat = P_O1)
dev.off()

# Omega 2
png(file = file.path("..", 'figures', "Omega2_Correlations.png"), width = 800, height = 800)
corrplot(COR_O2,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black', p.mat = P_O2)
dev.off()

# Epsilon 1
png(file = file.path("..", 'plots', "Epsilon1_Correlations.png"), width = 800, height = 800)
corrplot(COR_E1,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black', p.mat = P_E1)
dev.off()

# Epsilon 2
png(file = file.path("..", 'plots', "Epsilon2_Correlations.png"), width = 800, height = 800)
corrplot(COR_E2,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black', p.mat = P_E2)
dev.off()

# can feed p.mat and significance level to cross out non-significant
# correlations, e.g. insig = 'blank'

## As a pair of plots - omegas
png(file = file.path("..", 'figures', "Figure 1 - Omega1Omega2_Correlations.png"), width = 1600, height = 800)
par(mfrow= c(1,2))
corrplot(COR_O1,  order="hclust",addCoef.col = 'grey20',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = "grey90" , p.mat = P_O1, 
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(a) Spatial Encounter probability')
corrplot(COR_O2,  order="hclust",addCoef.col = 'grey20',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = 'grey90' , p.mat = P_O2, 
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(b) Spatial Density')
dev.off()

## Omegas with CIs

png(file = file.path("..", 'figures', "Figure 1 - Omega1Omega2_Correlations_Conf.png"), width = 1600, height = 800)
par(mfrow= c(1,2))
corrplot(COR_O1,  order="hclust", 
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = "grey90" , p.mat = P_O1, low = Lo_O1, upp = Up_O1, insig = 'blank',
		  cl.pos = 'n', plotC = "rect",addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(a) Spatial Encounter probability')
corrplot(COR_O2,  order="hclust",
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = 'grey90' , p.mat = P_O2, low = Lo_O2, upp = Up_O2, insig = 'blank', 
		  cl.pos = 'n', plotC = "rect", addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(b) Spatial Density')
dev.off()

## Omegas with non-sig blank

## Hack: need to pass a vector of colours, identifying the non-sig as the same colour
## as the background
cols_o1 <- corrplot(COR_O1, order = 'hclust')
# need the significant values in the plot order
for (r in unique(rownames(cols_o1))) {
for(c in unique(rownames(cols_o1))) {
cols_o1[rownames(cols_o1) == r, colnames(cols_o1)== c]  <- P_O1[rownames(P_O1) == r, colnames(P_O1) == c] 
}}
cols_o1 <- ifelse(cols_o1 < 0.05, "grey20", "grey90")
cols_o1 <- as.vector(cols_o1)

cols_o2 <- corrplot(COR_O2, order = 'hclust')
# need the significant values in the plot order
for (r in unique(rownames(cols_o2))) {
for(c in unique(rownames(cols_o2))) {
cols_o2[rownames(cols_o2) == r, colnames(cols_o2)== c]  <- P_O2[rownames(P_O2) == r, colnames(P_O2) == c] 
}}
cols_o2 <- ifelse(cols_o2 < 0.05, "grey20", "grey90")
cols_o2 <- as.vector(cols_o2)


png(file = file.path("..", 'figures', "Figure 1 - Omega1Omega2_Correlations_blank.png"), width = 1600, height = 800)
par(mfrow= c(1,2))
corrplot(COR_O1,  order="hclust" ,addCoef.col = cols_o1,
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = "grey90" , p.mat = P_O1, insig = 'blank',
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(a) Spatial Encounter probability', number.cex = 1.5)
corrplot(COR_O2,  order="hclust", addCoef.col = cols_o2,
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = 'grey90' , p.mat = P_O2, insig = 'blank', 
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(b) Spatial Density', number.cex = 1.5)
dev.off()

## As a pair of plots - epsilons
png(file = file.path("..", 'figures', "Epsilon1Epsilon2_Correlations.png"), width = 1600, height = 800)
par(mfrow= c(1,2))
corrplot(COR_E1,  order="hclust",addCoef.col = 'grey20',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = "grey90" , p.mat = P_E1, insig = 'blank',
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(a) Spatio-temporal Encounter probability')
corrplot(COR_E2,  order="hclust",addCoef.col = 'grey20',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = 'grey90' , p.mat = P_E1, insig = 'blank',
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(b) Spatio-temporal Density')
dev.off()



## Hack: need to pass a vector of colours, identifying the non-sig as the same colour
## as the background
cols_e1 <- corrplot(COR_E1, order = 'hclust')
# need the significant values in the plot order
for (r in unique(rownames(cols_e1))) {
for(c in unique(rownames(cols_e1))) {
cols_e1[rownames(cols_e1) == r, colnames(cols_e1)== c]  <- P_E1[rownames(P_E1) == r, colnames(P_E1) == c] 
}}
cols_e1 <- ifelse(cols_e1 < 0.05, "grey20", "grey90")
cols_e1 <- as.vector(cols_e1)

cols_e2 <- corrplot(COR_E2, order = 'hclust')
# need the significant values in the plot order
for (r in unique(rownames(cols_e2))) {
for(c in unique(rownames(cols_e2))) {
cols_e2[rownames(cols_e2) == r, colnames(cols_e2)== c]  <- P_E2[rownames(P_E2) == r, colnames(P_E2) == c] 
}}
cols_e2 <- ifelse(cols_e2 < 0.05, "grey20", "grey90")
cols_e2 <- as.vector(cols_e2)


png(file = file.path("..", 'figures', "Epsilon1Epsilon2_Correlations_blank.png"), width = 1600, height = 800)
par(mfrow= c(1,2))
corrplot(COR_E1,  order="hclust" ,addCoef.col = cols_e1,
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = "grey90" , p.mat = P_E1, insig = 'blank',
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(a) Spatio-temporal Encounter probability')
corrplot(COR_E2,  order="hclust", addCoef.col = cols_e2,
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = 'grey90' , p.mat = P_E2, insig = 'blank', 
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(b) Spatio-temporal Density')
dev.off()

#####################################

## Comparing the average vs the spatio-temporal correlations

##

## Encounter probability::
O1E1 <- reshape2::melt(COR_O1)
colnames(O1E1) <- c("spp1", "spp2", "O1")
E1 <- reshape2::melt(COR_E1)

O1E1$E1 <- E1$value

plot(O1E1$O1, O1E1$E1)
summary(lm(O1E1$E1 ~ O1E1$O1))
abline(lm(O1E1$E1 ~ O1E1$O1))

cor.test(O1E1$O1, O1E1$E1)

## Density::
O2E2 <- reshape2::melt(COR_O2)
colnames(O2E2) <- c("spp1", "spp2", "O2")
E2 <- reshape2::melt(COR_E2)

O2E2$E2 <- E2$value

plot(O2E2$O2, O2E2$E2)
summary(lm(O2E2$E2 ~ O2E2$O2))
anbline(lm(O2E2$E2 ~ O2E2$O2))

cor.test(O2E2$O2, O2E2$E2)

