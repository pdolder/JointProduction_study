############################################################
## Code for plotting the Catchability parameter estimates ##
############################################################

library(TMB)
library(RColorBrewer)
library(ggthemes)

#display.brewer.all(colorblindFriendly = T)
pal <- brewer.pal(6, 'Dark2')

#load(file.path('..', 'QEstimates.RData')) # Extract the 2 fixed parameter estimates - NO!! These are fixed to their adult equivilent

#QEst <- Qs

# List of survey and species
Survey <- c('CARLHELMAR','NWGFS','Q1SWBEAM','Q4SWIBTS','THA2', 'WCGFS')
spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')

#Qs_list <- lapply(Qs, summary, p.value = T)

#Qs <- as.data.frame(do.call(rbind, Qs_list))

#Qs$spp <- rep(spp, each = nrow(Qs_list[[1]]))
#Qs$Survey <- Survey

#Qs$Param <- rownames(Qs)

#Qs <- dplyr::filter(Qs, Param %in% c('lambda1_k', 'lambda2_k'),Survey == 'CARLHELMAR', spp %in% c('cod_juv', 'bud_juv'))

## Now extract the model estimates
load('Save.RData')

sd.sum <- summary(Save$Opt$SD, p.value = T)
sd.sum <- as.data.frame(sd.sum)
sd.sum$Param <- rownames(sd.sum) 
sd.sum <- dplyr::filter(sd.sum, Param %in% c('lambda1_k', 'lambda2_k'))

# The carlhelmar cod_juv and bud_juv positions are 2 and 6 for lambda1_k, as indicated in
# colnames(Save$TmbData$Q_ik), and 108 and 112 for lambda2_k

#Qs_all <- rbind(sd.sum[1,c(1:4)],  # 1
#		Qs[Qs$Survey == 'CARLHELMAR' & Qs$spp == 'cod_juv' & Qs$Param == 'lambda1_k',c(1:4)], # 2
#		sd.sum[2:4,c(1:4)], # 3:5
#		Qs[Qs$Survey == 'CARLHELMAR' & Qs$spp == 'bud_juv' & Qs$Param == 'lambda1_k',c(1:4)], #6
#		sd.sum[5:107,c(1:4)], # 7:109
#		Qs[Qs$Survey == 'CARLHELMAR' & Qs$spp == 'cod_juv' & Qs$Param == 'lambda2_k',c(1:4)], # 110
#		sd.sum[108:110,c(1:4)],  # 111:113
#		Qs[Qs$Survey == 'CARLHELMAR' & Qs$spp == 'bud_juv' & Qs$Param == 'lambda2_k',c(1:4)], # 114
#		sd.sum[111:212,c(1:4)]) 

Qs_all <- rbind(sd.sum[1,c(1:4)],  # 1
		sd.sum[1,c(1:4)],  # cod juv lambda 1 
		sd.sum[2:4,c(1:4)], # 3:5
		sd.sum[4,c(1:4)], # bud juv lambda 1
		sd.sum[5:107,c(1:4)], # 7:109
		sd.sum[107,c(1:4)], # cod juv lambda 2
		sd.sum[108:110,c(1:4)],  # 111:113
		sd.sum[110,c(1:4)],  # bud juv lambda 2
		sd.sum[111:212,c(1:4)]) 


Qs_all$spp <- rep(spp, each = 1)
Qs_all$Survey <- rep(Survey, each = 18)
Qs_all$Param <- rep(c('lambda1_k', 'lambda2_k'), each = 108)

Qs_all$concat <- colnames(Save$TmbData$Q_ik) ## make sure we've put labels in the right places...

# cod juv and bud juv don't have std.errors
Qs_all[Qs_all$Survey == 'CARLHELMAR' & Qs_all$spp %in% c('cod_juv','bud_juv'),2]  <- NA # ['std. Error']

dplyr::filter(Qs_all, Survey == 'CARLHELMAR', spp %in% c('cod_juv','cod_adu', 'bud_juv','bud_adu'))

head(Qs_all)  ## check params in right positions 
## YES



library(ggplot2)
library(dplyr)

colnames(Qs_all)[2] <- 'Std.Error'

ggplot(filter(Qs_all, Param %in% c('lambda1_k','lambda2_k')), aes(x = spp, y = Estimate)) +
	coord_flip() + facet_wrap(~Param) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	geom_point(data = filter(Qs_all, Survey =='CARLHELMAR', spp %in% c('cod_juv','bud_juv')), aes(y = Estimate), colour = 'red', shape = 2) +
	theme_bw() + geom_hline(aes(yintercept = 0), linetype = 'dotdash') + theme(legend.position = 'top')

ggsave('QEstimatesALL.png', width = 9, height = 12)

ggplot(filter(Qs_all, Param %in% c('lambda1_k')), aes(x = Survey, y = Estimate)) +
	 facet_wrap(~spp, scale = 'free_y', nrow = 6) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	geom_point(data = filter(Qs_all, Survey =='CARLHELMAR', Param == 'lambda1_k', spp %in% c('cod_juv','bud_juv')), aes(y = Estimate), colour = 'red', shape = 2) +
	theme_bw() + ggtitle('Encounter Prob') + theme(axis.text.x = element_text(angle = -90, hjust = 0)) + geom_hline(aes(yintercept = 0), linetype = 'dotdash')
ggsave('QEstimatesGridEnc.png', width = 12, height = 15)

ggplot(filter(Qs_all, Param %in% c('lambda2_k')), aes(x = Survey, y = Estimate)) +
	 facet_wrap(~spp, scale = 'free_y', nrow = 6) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	geom_point(data = filter(Qs_all, Survey =='CARLHELMAR', Param == 'lambda2_k', spp %in% c('cod_juv','bud_juv')), aes(y = Estimate), colour = 'red', shape = 2) +
	theme_bw() + ggtitle('Positive catch rates') + theme(axis.text.x = element_text(angle = -90, hjust = 0))+ geom_hline(aes(yintercept = 0), linetype = 'dotdash')

ggsave('QEstimatesGridPos.png', width = 12, height = 15)

