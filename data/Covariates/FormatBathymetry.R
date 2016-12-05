
BD <- read.csv('B3.DAT', sep = ';')

# restrict to 48,52 and -12, -2
library(dplyr)

colnames(BD) <- c('Lon','Lat','Depth','X')

BD2 <- filter(BD, Lon <= -2, Lon >= -12, Lat >= 48, Lat <= 52)

BD <- BD2[c('Lon', 'Lat', 'Depth')]

save(BD, file = 'Bathy.RData')
