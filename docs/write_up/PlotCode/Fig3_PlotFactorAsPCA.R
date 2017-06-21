
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
## Plots
library(ggplot2)
library(ggthemes)
library(cowplot)

dem <- paste(rep(c('cod', 'had', 'whg'), each = 2),rep(c('adu','juv'),times = 3), sep =  '_')
fla <- paste(rep(c('ple', 'sol'), each = 2),rep(c('adu','juv'),times = 2), sep =  '_')
dep <- paste(rep(c('hke', 'bud', 'pis', 'meg'), each = 2),rep(c('adu','juv'),times = 4), sep =  '_')


PCA.DFs <- lapply(Var_list, function(x) {
			  DF <- as.data.frame(x$L_pj_rot)
			  colnames(DF) <- paste('Factor',1:9, sep = '_')
			  DF$group <- ifelse(rownames(DF) %in% dem, 'roundfish', ifelse(rownames(DF) %in% fla, 'flatfish', 'deep'))
			  return(DF)
	     })

grps <- data.frame(Factor_1 = rep(-1.1,3), Factor_2 = c(-0.9, -0.95, -1), group = c('roundfish', 'flatfish', 'deep'))

#p1 <- ggplot(PCA.DFs[['Epsilon1']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
#	expand_limits(y = c(-1,0.5), x = c(-1,0.5)) + 	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw() + 
#	theme(legend.position = 'none', axis.text = element_text(face = 'bold')) + xlab('Factor 1 (39.8 % variance explained)') + ylab('Factor 2 (16.7 % variance explained)') +
#	geom_rect(aes(xmin  = -1.25, xmax = -0.95, ymin = -1.05, ymax = -0.85), fill= 'white', colour = 'grey') +
#	geom_label(data = grps, aes(label = group, fill = group), alpha = 0.2) 
#
#p2 <- ggplot(PCA.DFs[['Epsilon2']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
#	expand_limits(y = c(-1,0.5), x = c(-1,0.5)) + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw()  +
#	theme(legend.position = 'none', axis.text = element_text(face = 'bold'))+ xlab('Factor 1 (43.3 % variance explained)') + ylab('Factor 2 (13.9 % variance explained)')
#
#
#plot_grid(p1,p2, labels = c('(a) Factor loading for spatio-temporal \n        encounter probability',
#			    '(b) Factor loading for spatio-temporal \n            catch rates'), hjust = -0.2, vjust = 1.5)
#ggsave(file = 'PCAstyle_Plots_Spatiotemporal.png', width = 16, height = 8)

## Spatial 

#p3 <- ggplot(PCA.DFs[['Omega1']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
#	expand_limits(y = c(-2.6,0.1), x = c(-2.5,3.5)) +	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw() + 
#	theme(legend.position = 'none', axis.text = element_text(face = 'bold')) + xlab('Factor 1 (41 % variance explained)') + ylab('Factor 2 (31.3 % variance explained)') +
#	geom_rect(aes(xmin  = -1.25, xmax = -0.95, ymin = -1.05, ymax = -0.85), fill= 'white', colour = 'grey') +
#	geom_label(data = grps, aes(label = group, fill = group), alpha = 0.2) 
#
#p4 <- ggplot(PCA.DFs[['Omega2']], aes(x = Factor_1, y = Factor_2)) + geom_point() + geom_label(aes(label = rownames(PCA.DFs[[1]]), fill = group), nudge_y = 0.03, nudge_x = 0.06, alpha = 0.2) +
#	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw()  +
#	theme(legend.position = 'none', axis.text = element_text(face = 'bold'))+ xlab('Factor 1 (35.6 % variance explained)') + ylab('Factor 2 (18.1 % variance explained)')
#
#plot_grid(p3,p4, labels = c('(a) Factor loading for spatial \n        encounter probability',
#			    '(b) Factor loading for spatial \n            catch rates'), hjust = -0.2, vjust = 1.5)
#ggsave(file = 'PCAstyle_Plots_Spatial.png', width = 16, height = 8)

#}


## Now try with images
library(png)
library(gridGraphics)

fish_pics <- file.path('..','..','fish_pics')

img1 <- readPNG(file.path(fish_pics, "plaice-col.png"))
ple <- rasterGrob(img1, interpolate=T, width = 1)

img2 <- readPNG(file.path(fish_pics,"cod.png"))
cod <- rasterGrob(img2, interpolate=T, width = 1)

img3 <- readPNG(file.path(fish_pics,"haddock.png"))
had  <- rasterGrob(img3, interpolate=T, width = 1)

img4 <- readPNG(file.path(fish_pics,"whiting.png"))
whg  <- rasterGrob(img4, interpolate=T, width = 1)

img5 <- readPNG(file.path(fish_pics,"sole.png"))
sol  <- rasterGrob(img5, interpolate=T, width = 1)

img6 <- readPNG(file.path(fish_pics,"hake.png"))
hke  <- rasterGrob(img6, interpolate=T, width = 1)

img7 <- readPNG(file.path(fish_pics,"megrim.png"))
meg  <- rasterGrob(img7, interpolate=T, width = 1)

img8 <- readPNG(file.path(fish_pics,"bud.png"))
bud  <- rasterGrob(img8, interpolate=T, width = 1)

img9 <- readPNG(file.path(fish_pics,"pisc.png"))
pisc  <- rasterGrob(img9, interpolate=T, width = 1)


#b1 +  annotation_custom(g1, xmin=-0.7, xmax=-0.6, ymin=0.2, ymax=0.3) 

ysize_adu <- dim(img1)[1]/5000
xsize_adu <- dim(img1)[2]/5000

ysize_juv <- dim(img1)[1]/10000
xsize_juv <- dim(img1)[2]/10000

#################################

# g1: spatio-temp encounter

################################

# dataframe with positions
LinesDFEp1 <- data.frame(x1 = PCA.DFs[["Epsilon1"]][["Factor_1"]], 
			 y1 = PCA.DFs[["Epsilon1"]][["Factor_2"]],
			 x2 = PCA.DFs[["Epsilon1"]][["Factor_1"]],
			 y2 = PCA.DFs[["Epsilon1"]][["Factor_2"]])

rownames(LinesDFEp1) <- rownames(PCA.DFs[[1]])

LinesDFEp1$x2 <- ifelse(LinesDFEp1$x2 > 0, LinesDFEp1$x2 + 0.1 * LinesDFEp1$x2,
			                   LinesDFEp1$x2 - 0.1 * LinesDFEp1$x2)

LinesDFEp1$y2 <- ifelse(LinesDFEp1$y2 > 0, LinesDFEp1$y2 + 0.1 * LinesDFEp1$y2,
			                   LinesDFEp1$y2 - 0.1 * LinesDFEp1$y2)

# flip some to avoid overlaps
LinesDFEp1$y2[rownames(LinesDFEp1) %in% c("meg_adu","pisc_juv","hke_juv")] <-
 LinesDFEp1$y2[rownames(LinesDFEp1) %in% c("meg_adu","pisc_juv","hke_juv")] + -0.1

LinesDFEp1$x2[rownames(LinesDFEp1) %in% c("meg_adu")] <- LinesDFEp1$x2[rownames(LinesDFEp1) %in% c("meg_adu")] + -0.1
LinesDFEp1$x2[rownames(LinesDFEp1) %in% c("hke_juv")] <- LinesDFEp1$x2[rownames(LinesDFEp1) %in% c("hke_juv")] + 0.1
LinesDFEp1$x2[rownames(LinesDFEp1) %in% c("bud_adu")] <- LinesDFEp1$x2[rownames(LinesDFEp1) %in% c("bud_adu")] + 0.1

LinesDFEp1$y2[rownames(LinesDFEp1) %in% c("whg_adu")]  <- LinesDFEp1$y2[rownames(LinesDFEp1) %in% c("whg_adu")] +0.1
LinesDFEp1$y2[rownames(LinesDFEp1) %in% c("bud_adu")] <- LinesDFEp1$y2[rownames(LinesDFEp1) %in% c("bud_adu")] - 0.1


# Base plot with points
b1 <- ggplot(PCA.DFs[['Epsilon1']], aes(x = Factor_1, y = Factor_2)) + geom_point()  +
	expand_limits(y = c(-1,0.5), x = c(-1,0.5)) + 	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw() + 
	xlab('Factor 1 (39.8 % variance explained)') + ylab('Factor 2 (16.7 % variance explained)') 

g1 <- b1  + geom_segment(data = LinesDFEp1, aes(x = x1, xend = x2, y = y1, yend = y2)) +
#	geom_label(aes(label = rownames(LinesDFEp1)))# +
# cod adult
annotation_custom(cod, xmin = -0.32, xmax = -0.32 + xsize_adu, 
		      ymin = 0.1, ymax = 0.1 + ysize_adu) + 
# cod juvenile
annotation_custom(cod, xmin = -0.9, xmax = -0.9 + xsize_juv, 
		      ymin = 0.38, ymax = 0.38 + ysize_juv) +
# haddock adult
annotation_custom(had, xmin = -0.7, xmax = -0.7 + xsize_adu, 
		      ymin = -0.2, ymax = -0.2 + ysize_adu) + 
# haddock juvenile
annotation_custom(had, xmin = -1.05, xmax=-1.05 + xsize_juv, 
		      ymin = 0.25, ymax=0.25 + ysize_juv) + 
# whiting adult 
annotation_custom(whg, xmin = -0.85, xmax = -0.85 + xsize_adu, 
		      ymin = 0.08, ymax = 0.08 + ysize_adu) +
# whiting juvenile
annotation_custom(whg, xmin = -0.95, xmax = -0.95 + xsize_juv, 
		      ymin =0.22, ymax = 0.22 + ysize_juv) +
# hake adult
annotation_custom(hke, xmin = -0.28, xmax = -0.28 + xsize_adu, 
		      ymin = -0.54, ymax = -0.54 + ysize_adu) +
# hake juvenile
annotation_custom(hke, xmin = -0.25, xmax = -0.25 + xsize_juv,
		      ymin = -0.67, ymax = -0.67 + ysize_juv) +
# megrim adult
annotation_custom(meg, xmin = -0.65, xmax = -0.65 + xsize_adu, 
		      ymin = -0.70, ymax = -0.70 + ysize_adu) +
# megrim juvenile
annotation_custom(meg, xmin = -0.26, xmax= -0.26 + xsize_juv, 
		      ymin = -0.17, ymax= -0.17 + ysize_juv) +
# bud adult
annotation_custom(bud, xmin = -0.1, xmax = -0.1 + xsize_adu, 
		      ymin = -0.65, ymax = -0.65 + ysize_adu) +
# bud juvenile
annotation_custom(bud, xmin = -0.09, xmax = -0.09 + xsize_juv, 
		      ymin = -0.7, ymax = -0.7 + ysize_juv) +
# pis adult
annotation_custom(pisc, xmin = -0.4, xmax = -0.4 + xsize_adu,
		      ymin = -0.45, ymax = -0.45 + ysize_adu) +
# pis juvenile
annotation_custom(pisc, xmin = -0.33, xmax = -0.33 + xsize_juv, 
		      ymin = -0.71, ymax = -0.71 + ysize_juv) +
# plaice adult
annotation_custom(ple, xmin = -0.53, xmax = -0.53 + xsize_adu, 
		      ymin = -0.01, ymax = -0.01 + ysize_adu) + 
# plaice juvenile
annotation_custom(ple, xmin = -0.56, xmax = -0.56 + xsize_juv, 
		      ymin = 0.07, ymax = 0.07 + ysize_juv) + 
# sole adult
annotation_custom(sol, xmin = -0.47, xmax = -0.47 + xsize_adu, 
		      ymin = -0.15, ymax = -0.15 + ysize_adu) + 
# sole juvenile
annotation_custom(sol, xmin = -0.4, xmax = -0.4 + xsize_juv, 
		      ymin = 0.22, ymax = 0.22 + ysize_juv)  


#################################

# g2: spatio-temp positive 

################################

# dataframe with positions
LinesDFEp2 <- data.frame(x1 = PCA.DFs[["Epsilon2"]][["Factor_1"]], 
			 y1 = PCA.DFs[["Epsilon2"]][["Factor_2"]],
			 x2 = PCA.DFs[["Epsilon2"]][["Factor_1"]],
			 y2 = PCA.DFs[["Epsilon2"]][["Factor_2"]])

rownames(LinesDFEp2) <- rownames(PCA.DFs[[2]])

LinesDFEp2$x2 <- ifelse(LinesDFEp2$x2 > 0, LinesDFEp2$x2 + 0.1 * LinesDFEp2$x2,
			                   LinesDFEp2$x2 - 0.1 * LinesDFEp2$x2)

LinesDFEp2$y2 <- ifelse(LinesDFEp2$y2 > 0, LinesDFEp2$y2 + 0.1 * LinesDFEp2$y2,
			                   LinesDFEp2$y2 - 0.1 * LinesDFEp2$y2)

# avoid overlaps
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("bud_juv")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("bud_juv")] + 0.05
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("bud_adu")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("bud_adu")] + 0.05
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("pisc_adu")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("pisc_adu")] + 0.05
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("pisc_adu")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("pisc_adu")] - 0.05
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("meg_adu")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("meg_adu")] + 0.1
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("meg_adu")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("meg_adu")] + 0.05
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("hke_adu")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("hke_adu")] - 0.02
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("pisc_juv")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("pisc_juv")] - 0.1
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("pisc_juv")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("pisc_juv")] - 0.02
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("ple_adu")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("ple_adu")] - 0.05
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("ple_juv")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("ple_juv")] - 0.1
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("ple_juv")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("ple_juv")] + 0.15
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("sol_adu")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("sol_adu")] - 0.1
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("sol_adu")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("sol_adu")] - 0.1
LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("sol_juv")] <- LinesDFEp2$x2[rownames(LinesDFEp2) %in% c("sol_juv")] - 0.13
LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("sol_juv")] <- LinesDFEp2$y2[rownames(LinesDFEp2) %in% c("sol_juv")] - 0.14


# Base plot with points
b2 <- ggplot(PCA.DFs[['Epsilon2']], aes(x = Factor_1, y = Factor_2)) + geom_point()  +
	expand_limits(y = c(-1,0.5), x = c(-1,0.5)) + 	geom_vline(xintercept = 0) + geom_hline(yintercept = 0) + theme_bw() + 
	xlab('Factor 1 (43.3 % variance explained)') + ylab('Factor 2 (13.9 % variance explained)')

g2 <-  b2  + geom_segment(data = LinesDFEp2, aes(x = x1, xend = x2, y = y1, yend = y2)) + 
#	geom_text(aes(label = rownames(LinesDFEp1))) #+
# cod adult
annotation_custom(cod, xmin = -0.18, xmax = -0.18  + xsize_adu, 
		      ymin = -0.15, ymax = -0.15 + ysize_adu) + 
# cod juvenile
annotation_custom(cod, xmin = -0.26, xmax = -0.26 + xsize_juv, 
		      ymin = -0.16, ymax = -0.16 + ysize_juv) +
# haddock adult
annotation_custom(had, xmin = -0.6, xmax = -0.6 + xsize_adu, 
		      ymin = -0.28, ymax = -0.28 + ysize_adu) + 
# haddock juvenile
annotation_custom(had, xmin = -0.82, xmax= -0.82 + xsize_juv, 
		      ymin = -0.05, ymax= -0.05 + ysize_juv) + 
# whiting adult 
annotation_custom(whg, xmin = -0.8 , xmax = -0.8 + xsize_adu, 
		      ymin = -0.2, ymax = -0.2 + ysize_adu) +
# whiting juvenile
annotation_custom(whg, xmin = -0.7 , xmax = -0.7 + xsize_juv, 
		      ymin = 0.03, ymax = 0.03 + ysize_juv) +
# hake adult
annotation_custom(hke, xmin = -0.12, xmax = -0.12 + xsize_adu, 
		      ymin = -0.02, ymax = -0.02 + ysize_adu) +
# hake juvenile
annotation_custom(hke, xmin = -0.09, xmax = -0.09 + xsize_juv,
		      ymin = 0.07, ymax = 0.07 + ysize_juv) +
# megrim adult
annotation_custom(meg, xmin = -0.08, xmax = -0.08 + xsize_adu, 
		      ymin = 0.3, ymax = 0.3 + ysize_adu) +
# megrim juvenile
annotation_custom(meg, xmin = -0.1, xmax=  -0.1 + xsize_juv, 
		      ymin = 0.4 , ymax= 0.4 + ysize_juv) +
# bud adult
annotation_custom(bud, xmin = 0.09, xmax = 0.09 + xsize_adu, 
		      ymin = 0.15, ymax = 0.15 + ysize_adu) +
# bud juvenile
annotation_custom(bud, xmin = 0.07, xmax = 0.07 + xsize_juv, 
		      ymin = 0.28, ymax = 0.28 + ysize_juv) +
# pis adult
annotation_custom(pisc, xmin = 0.05, xmax = 0.05 + xsize_adu,
		      ymin = 0.07, ymax = 0.07 + ysize_adu) +
# pis juvenile
annotation_custom(pisc, xmin = -0.12, xmax = -0.12 + xsize_juv, 
		      ymin = 0.12, ymax = 0.12 + ysize_juv) +
# plaice adult
annotation_custom(ple, xmin = -0.4, xmax = -0.4  + xsize_adu, 
		      ymin = 0.3, ymax = 0.3 + ysize_adu) + 
# plaice juvenile
annotation_custom(ple, xmin = -0.5, xmax = -0.5 + xsize_juv, 
		      ymin = 0.45, ymax = 0.4 + ysize_juv) + 
# sole adult
annotation_custom(sol, xmin = -0.4, xmax = -0.4 + xsize_adu, 
		      ymin = 0.1, ymax = 0.1 + ysize_adu) + 
# sole juvenile
annotation_custom(sol, xmin = -0.3, xmax = -0.3 + xsize_juv, 
		      ymin = 0.08, ymax = 0.08 + ysize_juv)  

plot_grid(g1,g2, labels = c('(a) Factor loading for spatio-temporal encounter probability',
			    '(b) Factor loading for spatio-temporal catch rates'), hjust = -0.15, vjust = 2)
#ggsave(file = 'PCAstyle_Plots_Spatial_FishPics.png', width = 16, height = 8)

## Provide a legend

leg <- ggplot(data.frame(x = rep(0,1), y = rep(0,1)), aes(x= x, y = y))  + xlim(-0.22,0.3) + ylim(0,1) +  geom_point(colour = "white") +
	annotation_custom(cod, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.9, ymax = 0.9 + ysize_adu) + geom_text(x = 0.2, y = 1, label = "cod", size = 4) + 
	annotation_custom(had, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.8, ymax = 0.8 + ysize_adu) + geom_text(x = 0.2, y = 0.9, label = "haddock", size = 4) +
	annotation_custom(whg, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.7, ymax = 0.7 + ysize_adu) + geom_text(x = 0.2, y = 0.8, label = "whiting", size = 4) +
	annotation_custom(hke, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.6, ymax = 0.6 + ysize_adu) + geom_text(x = 0.2, y = 0.7, label = "hake", size = 4) +
	annotation_custom(meg, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.5, ymax = 0.5 + ysize_adu) + geom_text(x = 0.2, y = 0.6, label = "megrim", size = 4) +
	annotation_custom(pisc, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.37, ymax = 0.37 + ysize_adu) + geom_text(x = 0.2, y = 0.45, label = "monkfish \n (piscatorius)", size = 4) +
	annotation_custom(bud, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.25, ymax = 0.25 + ysize_adu) + geom_text(x = 0.2, y = 0.32, label = "monkfish \n (budegassa)", size = 4) +
	annotation_custom(ple, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.12, ymax = 0.12 + ysize_adu) + geom_text(x = 0.2, y = 0.2, label = "plaice", size = 4) +
	annotation_custom(sol, xmin = -0.2, xmax = -0.2 + xsize_adu, 
			  ymin = 0.01, ymax = 0.01 + ysize_adu) + geom_text(x = 0.2, y = 0.1, label = "sole", size = 4)  +
	      xlab("") + ylab("") + theme_void()

plot_grid(g1,g2,leg, ncol = 3, rel_widths =c(8,8,3),  labels = c('(a) Factor loading for spatio-temporal encounter probability',
			    '(b) Factor loading for spatio-temporal catch rates',''), hjust = -0.15, vjust = 2)
ggsave(file =  file.path('..', 'figures', 'Figure 3 - PCAstyle_Plots_SpatioTemp.png'), width = 16, height = 8)


#### END ###

