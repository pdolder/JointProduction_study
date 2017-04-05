
library(TMB)
library(RColorBrewer)
library(ggthemes)

display.brewer.all(colorblindFriendly = T)
pal <- brewer.pal(6, 'Dark2')

load('QEstimates.RData')

QEst <- Qs

# List of survey and species
Survey <- c('CEXP','NWGFS','Q1SWBEAM','Q4SWIBTS','THA2', 'WCGFS')
spp <- paste(rep(c('cod','meg','bud','pisc','had','whg','hke','ple','sol'), each = 2),c('adu','juv'), sep = '_')

Qs_list <- lapply(Qs, summary, p.value = T)

Qs <- as.data.frame(do.call(rbind, Qs_list))

Qs$spp <- rep(spp, each = nrow(Qs_list[[1]]))
Qs$Survey <- Survey

Qs$Param <- rownames(Qs)

Qs <- dplyr::filter(Qs, Param %in% c('lambda1_k', 'lambda2_k'))

library(ggplot2)
library(dplyr)

colnames(Qs)[2] <- 'Std.Error'

ggplot(filter(Qs, Param %in% c('lambda1_k','lambda2_k')), aes(x = spp, y = Estimate)) +
	coord_flip() + facet_wrap(~Param) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	theme_bw()

ggsave('QEstimatesALL.png', width = 6, height = 12)


ggplot(filter(Qs, Param %in% c('lambda1_k','lambda2_k'), !spp %in% c('cod_juv','bud_juv')), aes(x = spp, y = Estimate)) +
	 coord_flip() + facet_wrap(~Param) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	theme_bw()

ggsave('QEstimates.png', width = 6, height = 12)

ggplot(filter(Qs, Param %in% c('lambda1_k')), aes(x = Survey, y = Estimate)) +
	 facet_wrap(~spp, scale = 'free', nrow = 6) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	theme_bw() + ggtitle('Encounter Prob') + theme(axis.text.x = element_text(angle = -90))
ggsave('QEstimatesGridEnc.png', width = 12, height = 15)

ggplot(filter(Qs, Param %in% c('lambda2_k')), aes(x = Survey, y = Estimate)) +
	 facet_wrap(~spp, scale = 'free', nrow = 6) + scale_color_manual(values = pal) +
	geom_pointrange(aes(y = Estimate, ymin = Estimate - 1.96 * Std.Error, ymax = Estimate + 1.96 * Std.Error, colour = Survey)) +
	theme_bw() + ggtitle('Positive catch rates') + theme(axis.text.x = element_text(angle = -90))

ggsave('QEstimatesGridPos.png', width = 12, height = 15)

