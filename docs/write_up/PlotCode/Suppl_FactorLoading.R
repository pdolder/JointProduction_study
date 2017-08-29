
# Load libraries
library(VAST)

# This is where all runs will be located
run <- '2017-06-16_M1'

load(file.path('..', '..', '..', 'results', run, 'Save.RData'))
load(file.path('..', '..', '..', 'results', 'CovariatesAtKnot.RData'))

DF <- Save$Data_Geostat
spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')

###################

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
}

PCA.DFs <- lapply(Var_list, function(x) {
			  DF <- as.data.frame(x$L_pj_rot)
			  colnames(DF) <- paste('Factor',1:9, sep = ' ')
			  return(DF)
	     })

## explained variance
## sum(L_pj_rot[f]^2)/sum(L_pj_rot^2)

params <- names(Var_list)

var.DF <- lapply(params, function(y) {
DF <- data.frame(variable = y, 
		 factor = paste("Factor", 1:9, sep = " "),
		 values = sapply(1:9, function(x) {round(100 * sum(Var_list[[y]]$L_pj_rot[,x]^2) / sum(Var_list[[y]]$L_pj_rot^2),1)})
)
		 return(DF)
})

var.DF <- do.call(rbind, var.DF)

## Plots
library(ggplot2)
library(dplyr)
library(ggthemes)
## Plot the factor loadings (first 3)

## Omega 1

O1 <- data.frame(species = rownames(PCA.DFs[["Omega1"]]), reshape2::melt(PCA.DFs[["Omega1"]]))
colnames(O1)[2] <- "factor"

# Colour for the labels and bars
O1$col <- ifelse(O1$value > 0, "blue", "red")

ggplot(filter(O1, factor %in% paste("Factor", 1:3, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col), colour = "black") + facet_wrap(~ factor) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none") +
	ggtitle("average spatial encounter probability") + scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") + 
	geom_text(data = filter(var.DF, variable == "Omega1", factor %in% paste("Factor", 1:3, sep = " ")), aes(x = 7, y = 3,label = paste("Variance explained = ", values, "%")), inherit.aes  =F, size = 3)

ggsave(file = file.path('..', 'figures', 'Suppl_FactorLoading_Omega1.png'), width = 12, height = 4)



## Omega 2

O2 <- data.frame(species = rownames(PCA.DFs[["Omega2"]]), reshape2::melt(PCA.DFs[["Omega2"]]))
colnames(O2)[2] <- "factor"

# Colour for the labels and bars
O2$col <- ifelse(O2$value > 0, "blue", "red")

ggplot(filter(O2, factor %in% paste("Factor", 1:3, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col), colour = "black") + facet_wrap(~ factor ) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none") + 
	ggtitle("average spatial density") + scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") + 
	geom_text(data = filter(var.DF, variable == "Omega2", factor %in% paste("Factor", 1:3, sep = " ")), aes(x = 7, y = 3,label = paste("Variance explained = ", values, "%")), inherit.aes  =F, size = 3)


ggsave(file = file.path('..', 'figures', 'Suppl_FactorLoading_Omega2.png'), width = 12, height = 4)


## Epsilon 1

E1 <- data.frame(species = rownames(PCA.DFs[["Epsilon1"]]), reshape2::melt(PCA.DFs[["Epsilon1"]]))
colnames(E1)[2] <- "factor"

# Colour for the labels and bars
E1$col <- ifelse(E1$value > 0, "blue", "red")

ggplot(filter(E1, factor %in% paste("Factor", 1:3, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col), colour = "black") + facet_wrap(~ factor ) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none") +
	ggtitle("spatiotemporal encounter probability") + scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") + 
	geom_text(data = filter(var.DF, variable == "Epsilon1", factor %in% paste("Factor", 1:3, sep = " ")), aes(x = 7, y = 3,label = paste("Variance explained = ", values, "%")), inherit.aes  =F, size = 3)


Ggsave(file = file.path('..', 'figures', 'Suppl_FactorLoading_Epsilon1.png'), width = 12, height = 4)



## Epsilon 2

E2 <- data.frame(species = rownames(PCA.DFs[["Epsilon2"]]), reshape2::melt(PCA.DFs[["Epsilon2"]]))
colnames(E2)[2] <- "factor"

# Colour for the labels and bars
E2$col <- ifelse(E2$value > 0, "blue", "red")

ggplot(filter(E2, factor %in% paste("Factor", 1:3, sep = " ")), aes(x = species, y = value)) +
	geom_bar(stat = 'identity', aes(fill = col), colour = "black") + facet_wrap(~ factor ) +
 	theme_tufte() + theme(axis.text.x = element_text(angle = -90, hjust = 0), legend.position = "none") + 
	ggtitle("spatiotemporal density") + scale_fill_manual(values = c("red", "blue")) + ylab("") + xlab("") + 
	geom_text(data = filter(var.DF, variable == "Epsilon2", factor %in% paste("Factor", 1:3, sep = " ")), aes(x = 7, y = 3,label = paste("Variance explained = ", values, "%")), inherit.aes  =F, size = 3)


ggsave(file = file.path('..', 'figures', 'Suppl_FactorLoading_Epsilon2.png'), width = 12, height = 4)



