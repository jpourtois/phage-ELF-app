
library(ggplot2)

setwd("~/Documents/Stanford/Research /DLS project") 

# Figure 2: PhagoBurn DLS vs time

phago.dls <- read.csv('phago_processed.csv')

phago.dls$diff_control <- NaN

for (j in unique(phago.dls$day)) {
    
  for (k in unique(phago.dls$phage[phago.dls$day == j])) {
      
      
    phago.dls$diff_control[phago.dls$day == j & phago.dls$phage == k] <- phago.dls$value[phago.dls$day == j & phago.dls$phage == k] - 
      phago.dls$value[phago.dls$day == 0 & phago.dls$phage == k]
      
    
  }
}
phago.dls$diff_control <- abs(phago.dls$diff_control)

area.log <- tapply(phago.dls$diff_control,phago.dls$setting,sum)
area.log <- data.frame(setting = names(area.log), diff_AUC = area.log, row.names=NULL)
area.log$day <- sapply(strsplit(area.log$setting, '[.]'), getElement, 3)
area.log$phage <- sapply(strsplit(area.log$setting, '[.]'), getElement, 1)

# TEST TITER

titer <- 1:56

area.log$titer <- titer

ggplot(area.log, aes(x = diff_AUC, y = titer, group = phage))+
  geom_line(aes(color = phage)) +
  geom_point(aes(color = phage)) + 
  scale_y_continuous(trans='log10') +
  geom_smooth(method='lm',formula= y ~ x)

ggplot(area.log, aes(x = diff_AUC, y = titer))+
  geom_point() + 
  scale_y_continuous(trans='log10') +
  geom_smooth(method='lm',formula= y ~ x)

# Don't forget to have the titer difference logged for this
summary(lm(titer ~ diff_AUC, area.log))


# Figures

ggplot(area.log, aes(x = day, y = diff_AUC, group = phage)) + 
  geom_line(aes(color = phage))

luz19 <- area.log[area.log$phage %in% c('LUZ19','LUZ24'),]

ggplot(luz19, aes(x = day, y = diff_AUC, group = phage)) + 
  geom_line(aes(color = phage))

ggplot(luz19, aes(x = day, y = diff_AUC, group = 1)) + 
  geom_line()

luz19 <- phago.dls[phago.dls$phage == 'LUZ19',]
ggplot(luz19, aes(x = size, y = value, color = setting)) +
  geom_line() + 
  scale_x_continuous(trans='log10') 
