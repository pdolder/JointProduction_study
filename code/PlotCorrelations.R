## Loading matrix function
loadings_matrix <- function(L_val, n_rows, n_cols) {
	Dat <- matrix(NA, nrow = n_rows, ncol = n_cols)
	count <- 1
	for (r in 1:n_rows) {
		for (c in 1:n_cols) {
			if(r>=c) {
			Dat[r,c] <- L_val[count]
			count <- 1+count
			}
			else
			Dat[r,c] <- 0.0
		}
	}
	return(Dat)
}

## Convert upper cov to cor function
convert_upper_cov_to_cor <- function(cov) {
nrows <- nrow(cov)
ncols <- ncol(cov)
cor <- matrix(NA, ncol = ncols, nrow = nrows)
for(i in 1:nrows) {
	for(j in 1:ncols) {
	cor[i,j] <- cov[i,j] / cov[i,i]^0.5 / cov[j,j]^0.5
		}
}
return(cor)
}


## correlation function
return_cor <- function(L_omega1_z, n_c, FieldConfig = 9) {
L1_omega_cf <- loadings_matrix(L_omega1_z, n_c, FieldConfig)
lowercov_uppercor_omega1  <-  L1_omega_cf %*% t(L1_omega_cf)
lowercov_uppercor_omega1  <- convert_upper_cov_to_cor(lowercov_uppercor_omega1)
return(lowercov_uppercor_omega1)
}


load(file.path('..','results', '2017-04-08_M2','Save.RData'))

# Omega 1
L_omega1_z <- Save$Opt$SD$par.fixed[names(Save$Opt$SD$par.fixed) == 'L_omega1_z']
L_omega2_z <- Save$Opt$SD$par.fixed[names(Save$Opt$SD$par.fixed) == 'L_omega2_z']
L_epsilon1_z <- Save$Opt$SD$par.fixed[names(Save$Opt$SD$par.fixed) == 'L_epsilon1_z']
L_epsilon2_z <- Save$Opt$SD$par.fixed[names(Save$Opt$SD$par.fixed) == 'L_epsilon2_z']

COR_O1  <- return_cor(L_omega1_z, n_c = 18, 9)
COR_O2  <- return_cor(L_omega2_z, n_c = 18, 9)
COR_E1  <- return_cor(L_epsilon1_z, n_c = 18, 9)
COR_E2  <- return_cor(L_epsilon2_z, n_c = 18, 9)

library(corrplot)
library(RColorBrewer)

spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')

colnames(COR_O1)  <- rownames(COR_O1) <- colnames(COR_O2)  <- rownames(COR_O2) <- colnames(COR_E1)  <- rownames(COR_E1) <- colnames(COR_E2)  <- rownames(COR_E2) <- spp 

# Omega 1
png(file = file.path("..", 'plots', "Omega1_Correlations.png"), width = 800, height = 800)
corrplot(COR_O1,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black')
dev.off()

# Omega 2
png(file = file.path("..", 'plots', "Omega2_Correlations.png"), width = 800, height = 800)
corrplot(COR_O2,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black')
dev.off()

# Epsilon 1
png(file = file.path("..", 'plots', "Epsilon1_Correlations.png"), width = 800, height = 800)
corrplot(COR_E1,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black')
dev.off()

# Epsilon 2
png(file = file.path("..", 'plots', "Epsilon2_Correlations.png"), width = 800, height = 800)
corrplot(COR_E2,  order="hclust",addCoef.col = 'grey90',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),
		  cl.pos = 'n', addrect = 3, tl.col = 'black')
dev.off()

# can feed p.mat and significance level to cross out non-significant
# correlations, e.g. insig = 'blank'

## As a pair of plots
png(file = file.path("..", 'plots', "Omega1Omega2_Correlations.png"), width = 1600, height = 800)
par(mfrow= c(1,2))
corrplot(COR_O1,  order="hclust",addCoef.col = 'grey20',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = "grey90" ,
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(a) Spatial Encounter probability')
corrplot(COR_O2,  order="hclust",addCoef.col = 'grey20',
	          col=rev(brewer.pal(n=8, name="RdYlBu")),mar=c(0,0,1,0), bg = 'grey90' ,
		  cl.pos = 'n', addrect = 3, tl.col = 'black', tl.cex = 1.5, title = '(b) Spatial Density')
dev.off()
