

In <- read.csv('Table_for_SS3.csv')

library(ggplot2)

ggplot(In, aes(x = Year, y = Estimate_metric_tons)) + geom_line() + facet_wrap(~Category, scale = 'free_y', ncol = 4) +
  geom_ribbon(aes(ymin = Estimate_metric_tons - 1.96 * SD_mt, ymax = Estimate_metric_tons  + 1.96 * SD_mt), alpha = 0.4) +
  theme_classic()
ggsave(file = 'IndexGGplot.png', width = 10, height = 12)

IndSc <- apply(In$Estimate_metric_tons, In$Category, scale)

library(dplyr)

In <- In %>% group_by(Category) %>% mutate(InSc = scale(Estimate_metric_tons))

scale(In[c(In$Estimate_metric_tons, In$SD_mt),])

ggplot(In, aes(x = Year, y = InSc)) + geom_line() + facet_wrap(~Category, scale = 'free_y', ncol = 4) +
  geom_ribbon(aes(ymin = InSc - LoSc, ymax = InSc + UpSc), alpha = 0.4) +   theme_classic()
ggsave(file = 'IndexScaledGGplot.png', width = 10, height = 12)


In2 <- transform(In, In.norm = ave(Estimate_metric_tons, Category, FUN = scale))
In2 <- transform(In2, SD.norm = ave(SD_mt, Category, FUN = scale))

ggplot(In2, aes(x = Year, y = In.norm)) + geom_line() + facet_wrap(~Category, scale = 'free_y', ncol = 4) +
  geom_errorbar(aes(ymin = In.norm - SD.norm, ymax = In.norm + SD.norm), alpha = 0.4) +   theme_classic()
