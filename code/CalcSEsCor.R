## Trying to recover SEs for correlations without rerunning VAST

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

L_omega1_z <- Save$Opt$SD$par.fixed[names(Save$Opt$SD$par.fixed) == 'L_omega1_z']

library(TMB)
sd.sum <- summary(Save$Opt$SD)
sd.sum <- sd.sum[rownames(sd.sum)=='L_omega1_z',]

COR  <- return_cor(L_omega1_z, n_c = 18, 9)
COR_up <- return_cor(c(L_omega1_z + sd.sum[,2]), n_c = 18, 9)
COR_lo <- return_cor(c(L_omega1_z - sd.sum[,2]), n_c = 18, 9)


correlations <- data.frame(Estimate = as.vector(COR),
			   Lower    = as.vector(COR_lo),
			   Upper    = as.vector(COR_up))

