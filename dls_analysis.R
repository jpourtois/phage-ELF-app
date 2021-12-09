

# Load packages
# Note: you will need to install these packages first if it is your first
# time using them
library('tidyr')
library(ggplot2)
library(ggpubr)
library(mclust)
library(sBIC)

# Set working directory --> change path to where you saved the csv files
setwd("~/Documents/Stanford/Research /DLS project") 

# Load data
dls <- read.csv('dls.csv')
titer <- read.csv('titer.csv')

dls$X <- NULL
titer$X <- NULL

dls$setting[dls$setting == "OMKO_control_per"] <- "OMKO_control_peroxide"

# Calculate area under the curve

dls$diff_control <- NaN

for (i in unique(dls$measure)) {
  
  for (j in unique(dls$experiment)) {
    
    for (k in unique(dls$phage[dls$experiment == j])) {
      
      
      dls$diff_control[dls$measure == i & dls$experiment == j & dls$phage == k] <- dls$value[dls$measure == i & dls$experiment == j & dls$phage == k] - 
        dls$value[dls$measure == i & dls$experiment == j & dls$phage == k & dls$group == "control"]
      
    }
  }
}


dls$diff_control <- abs(dls$diff_control)

dls.int <- dls[dls$measure == 'intensity',]
area.log <- tapply(dls.int$diff_control,dls.int$setting,sum)
area.log <- data.frame(setting = names(area.log), value = area.log, row.names=NULL)

dls.vol <- dls[dls$measure == 'volume',]
area.log.vol <- tapply(dls.vol$diff_control,dls.vol$setting,sum)
area.log.vol <- data.frame(setting = names(area.log.vol), area.vol = area.log.vol, row.names=NULL)

dls.num <- dls[dls$measure == 'number',]
area.log.num <- tapply(dls.num$diff_control,dls.num$setting,sum)
area.log.num <- data.frame(setting = names(area.log.num), area.num = area.log.num, row.names=NULL)

dls.titer <- merge(area.log,titer, by = "setting")
dls.titer <- merge(dls.titer, area.log.vol, by = 'setting')
dls.titer <- merge(dls.titer, area.log.num, by = 'setting')

# Calculate mode
max.size.dls <- dls.int$size[tapply(dls.int$value, dls.int$setting, which.max)]

# Calculate mean
dls.int$mult <- log10(dls.int$size)*dls.int$value
mean.size <- tapply(dls.int$mult, dls.int$setting, sum)/tapply(dls.int$value, dls.int$setting, sum)

# Calculate median
median.size <- rep(0, length(unique(dls.int$setting)))
row = 1

for (i in dls.titer$setting) {
  
  median.size[row] <- median(rep(log10(dls.int$size[dls.int$setting == i]), dls.int$value[dls.int$setting == i]*10))
  
  row = row + 1
}

# Calculate variance

variance.size <- rep(0, length(unique(dls.int$setting)))
row = 1

for (i in dls.titer$setting) {
  
  variance.size[row] <- var(rep(log10(dls.int$size[dls.int$setting == i]), dls.int$value[dls.int$setting == i]*10))
  
  row = row + 1
}


# Calculate titer difference

dls.titer$titer[dls.titer$titer == 0] <- 1

dls.titer$titer.diff <- NaN
dls.titer$titer.diff.relative <- NaN
dls.titer$mode.diff <- NaN

for (i in unique(dls.titer$experiment)){
  
  for (j in unique(dls.titer$phage[dls.titer$experiment == i])) {
  
    # Titer difference
    dls.titer$titer.diff[dls.titer$experiment == i & dls.titer$phage == j] <- log10(dls.titer$titer[dls.titer$experiment == i & dls.titer$group == 'control' & dls.titer$phage == j]/(dls.titer$titer[dls.titer$experiment == i & dls.titer$phage == j]))
    
    # Relative titer difference
    dls.titer$titer.diff.relative[dls.titer$experiment == i & dls.titer$phage == j] <- dls.titer$titer.diff[dls.titer$experiment == i & dls.titer$phage == j]/log10(dls.titer$titer[dls.titer$experiment == i & dls.titer$group == 'control' & dls.titer$phage == j])
    
    #Difference in mode
    dls.titer$mode.diff[dls.titer$experiment == i & dls.titer$phage == j] <-  abs(max.size.dls[dls.titer$experiment == i & dls.titer$phage == j] - max.size.dls[dls.titer$experiment == i & dls.titer$phage == j & dls.titer$group == 'control'])
    
    #Difference in mean
    dls.titer$mean.diff[dls.titer$experiment == i & dls.titer$phage == j] <- abs(mean.size[dls.titer$experiment == i & dls.titer$phage == j] - mean.size[dls.titer$experiment == i & dls.titer$phage == j & dls.titer$group == 'control'])
    
    #Difference in median 
    dls.titer$median.diff[dls.titer$experiment == i & dls.titer$phage == j] <- abs(median.size[dls.titer$experiment == i & dls.titer$phage == j] - median.size[dls.titer$experiment == i & dls.titer$phage == j & dls.titer$group == 'control'])
    
    #Difference in variance
    dls.titer$variance.diff[dls.titer$experiment == i & dls.titer$phage == j] <- abs(variance.size[dls.titer$experiment == i & dls.titer$phage == j] - variance.size[dls.titer$experiment == i & dls.titer$phage == j & dls.titer$group == 'control'])
    
    
  }
}

##### Cluster analysis

# Example for Tuesday meeting

# This function adds up the different gaussians, to plot with the data
create.gm <- function(size.v, gm.model) {
  
  n.gaussian <- gm.model$G
  
  g.sum <- 0
  
  for (i in 1:n.gaussian) {
    
    g.sum <- g.sum + dnorm(log10(size.v), gm.model$parameters$mean[i], sqrt(gm.model$parameters$variance$sigmasq[i]))*gm.model$parameters$pro[i]
    
  }
  
  g.sum <- g.sum/sum(g.sum)
  
  return(g.sum) 
}

# LPS control, peroxide
gm.data.1 <- rep(log10(dls.int$size[dls.int$setting == 'LPS_control_peroxide']), dls.int$value[dls.int$setting == 'LPS_control_peroxide']*10)

gm.1 <- Mclust(gm.data.1)
create.gm(dls.int$size[dls.int$setting == 'LPS_control_peroxide'], gm.1)

gm.df.1 <- data.frame(size = log10(dls.int$size[dls.int$setting == 'LPS_control_peroxide']),
                      intensity = dls.int$value[dls.int$setting == 'LPS_control_peroxide']/sum(dls.int$value[dls.int$setting == 'LPS_control_peroxide']), 
                      gaussian = create.gm(dls.int$size[dls.int$setting == 'LPS_control_peroxide'], gm.1), row.names=NULL)

ggplot(gm.df.1, aes(x = size, y = intensity)) +
  geom_line()+ 
  geom_line(aes(x = size, y = gaussian), color = 'red')

# LPS, peroxide concentration = 2
gm.data.2 <- rep(log10(dls.int$size[dls.int$setting == 'LPS_2']), dls.int$value[dls.int$setting == 'LPS_2']*10)
gm.2 <- Mclust(gm.data.2)

create.gm(dls.int$size[dls.int$setting == 'LPS_2'], gm.2)

gm.df.2 <- data.frame(size = log10(dls.int$size[dls.int$setting == 'LPS_2']),
                      intensity = dls.int$value[dls.int$setting == 'LPS_2']/sum(dls.int$value[dls.int$setting == 'LPS_2']), 
                      gaussian = create.gm(dls.int$size[dls.int$setting == 'LPS_2'], gm.2), row.names=NULL)

ggplot(gm.df.2, aes(x = size, y = intensity)) +
  geom_line()+ 
  geom_line(aes(x = size, y = gaussian), color = 'red')

row <- 1
means <- list()

for (i in dls.titer$setting) {
  
  test <- rep(log10(dls.int$size[dls.int$setting == i]), dls.int$value[dls.int$setting == i]*10)
  gm <- Mclust(test)
  means[[row]] <- gm$parameters$mean
  row <- row + 1
  
}

test <- rep(log10(dls.int$size[dls.int$setting == 'LPS_2']), dls.int$value[dls.int$setting == 'LPS_2']*10)
plot(density(test))
plot(mclustBIC(test))

gm <- Mclust(test)
summary(gm, parameters = TRUE)

plot(gm, what = 'density')



##### PLot by experiment #####

# Peroxide experiment

dls.titer.peroxide <- dls.titer[dls.titer$experiment == 'peroxide',]

area.log.peroxide <- ggplot(dls.titer.peroxide, aes(x = value, y = -titer.diff, color = phage)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Peroxide experiment", color = "Phage") + 
  theme_minimal()

# Salt experiment

dls.titer.salt <- dls.titer[dls.titer$experiment == 'salt',]

area.log.salt  <- ggplot(dls.titer.salt, aes(x = value, y = -titer.diff, color = phage)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Salt experiment", color = "Phage") + 
  theme_minimal()

# Physical experiment

dls.titer.physical <- dls.titer[dls.titer$experiment == 'physical',]

area.log.physical  <- ggplot(dls.titer.physical, aes(x = value, y = -titer.diff, color = phage)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Sonication/vortex experiment", color = "Phage") + 
  theme_minimal()

# Lyo experiment

dls.titer.lyo <- dls.titer[dls.titer$experiment == 'lyo',]

area.log.lyo <- ggplot(dls.titer.lyo, aes(x = value, y = -titer.diff, color = phage)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Lyophilization experiment", color = "Phage") + 
  theme_minimal()

# Reduction experiment

dls.titer.red <- dls.titer[dls.titer$experiment == 'reducing',]

area.log.red <- ggplot(dls.titer.red, aes(x = value, y = -titer.diff, color = phage)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Reduction experiment", color = "Phage") + 
  theme_minimal()


ggarrange(area.log.peroxide, area.log.physical, area.log.salt, area.log.lyo, area.log.red)


##### Plot by phage ######
### OMKO phage 

dls.titer.OMKO <- dls.titer[dls.titer$phage == "OMKO",]

area.log.OMKO <- ggplot(dls.titer.OMKO, aes(x = value, y = -titer.diff, color = experiment)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "OMKO phage", color = "Experiment") + 
  theme_minimal()

### LPS phage 

dls.titer.LPS <- dls.titer[dls.titer$phage == "LPS",]

area.log.LPS <- ggplot(dls.titer.LPS, aes(x = value, y = -titer.diff, color = experiment)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "LPS phage", color = "Experiment") + 
  theme_minimal()

### T4P 

dls.titer.T4P <- dls.titer[dls.titer$phage == "T4P",]

area.log.T4P <- ggplot(dls.titer.T4P, aes(x = value, y = -titer.diff, color = experiment)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "T4P phage", color = "Experiment") + 
  theme_minimal()

### LUZ 19

dls.titer.luz19 <- dls.titer[dls.titer$phage == "Luz19" | dls.titer$phage == "LUZ19",]

area.log.luz19 <- ggplot(dls.titer.luz19, aes(x = value, y = -titer.diff, color = experiment)) +
  geom_point() +
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Luz19 phage", color = "Experiment") + 
  theme_minimal()


### LUZ24

dls.titer.luz24 <- dls.titer[dls.titer$phage == "Luz24" | dls.titer$phage == "LUZ24",]

area.log.luz24 <- ggplot(dls.titer.luz24, aes(x = value, y = -titer.diff, color = experiment)) +
  geom_point() + 
  ylim(c(-11,0.5)) +
  labs(x = "Area under the curve", y = "Difference in titer (log10)", title = "Luz24 phage", color = "Experiment") + 
  theme_minimal()

ggarrange(area.log.OMKO, area.log.LPS, area.log.luz19, area.log.luz24, area.log.T4P)


#### Compare different metrics ####

# Area under the curve, intensity

area.log.phage <- ggplot(dls.titer, aes(x = value, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

area.log.exp <- ggplot(dls.titer, aes(x = value, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(area.log.phage, area.log.exp)

summary(lm(titer.diff ~ value - 1, dls.titer))
# Area under the curve,volume

area.vol.phage <- ggplot(dls.titer, aes(x = area.vol, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

area.vol.exp <- ggplot(dls.titer, aes(x = area.vol, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(area.vol.phage, area.vol.exp)

summary(lm(titer.diff ~ area.vol - 1, dls.titer))

# Area under the curve,number

area.num.phage <- ggplot(dls.titer, aes(x = area.num, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

area.num.exp <- ggplot(dls.titer, aes(x = area.num, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(area.num.phage, area.num.exp)

summary(lm(titer.diff ~ area.num - 1, dls.titer))

# Area under the curve, no physical

dls.titer.nophy <- dls.titer[!dls.titer$experiment == 'physical',]

area.log.phage <- ggplot(dls.titer.nophy, aes(x = value, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

area.log.exp <- ggplot(dls.titer.nophy, aes(x = value, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(area.log.phage, area.log.exp)

summary(lm(titer.diff ~ value, dls.titer.nophy))

# Difference in mode

mode.phage <- ggplot(dls.titer, aes(x = mode.diff + 1, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Difference in mode", y = "Difference in titer (log10)") + 
  scale_x_continuous(trans='log10') +
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

mode.experiment <- ggplot(dls.titer, aes(x = mode.diff + 1, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Difference in mode", y = "Difference in titer (log10)") + 
  scale_x_continuous(trans='log10') +
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(mode.phage, mode.experiment)

summary(lm(titer.diff ~ log10(mode.diff+1), dls.titer))

# Difference in mean

mean.phage <- ggplot(dls.titer, aes(x = mean.diff, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Difference in mean", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

mean.experiment <- ggplot(dls.titer, aes(x = mean.diff, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Difference in mean", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(mean.phage, mean.experiment)

summary(lm(titer.diff ~ mean.diff, dls.titer))

# Difference in median

median.phage <- ggplot(dls.titer, aes(x = median.diff, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Difference in median", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

median.exp <- ggplot(dls.titer, aes(x = median.diff, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Difference in median", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(median.phage, median.exp)

summary(lm(titer.diff ~ median.diff, dls.titer))


# Difference in variance
variance.phage <- ggplot(dls.titer, aes(x = variance.diff, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Difference in variance", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

variance.exp <- ggplot(dls.titer, aes(x = variance.diff, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Difference in variance", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(variance.phage, variance.exp)

summary(lm(titer.diff ~ variance.diff, dls.titer))

# Combination of mean and variance

combi.phage <- ggplot(dls.titer, aes(x = variance.diff + mean.diff*2, y = -titer.diff)) +
  geom_point(aes(color = phage)) + 
  labs(x = "Difference in variance", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

combi.exp <- ggplot(dls.titer, aes(x = variance.diff + mean.diff*2, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(x = "Difference in variance", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

ggarrange(combi.phage, combi.exp)

summary(lm(titer.diff ~ variance.diff + I(mean.diff*2), dls.titer))


# case study for LPS and OMKO phage

dls.titer.nophy.LPS <- dls.titer[(!dls.titer$experiment == 'physical') & (dls.titer$phage == 'LPS'),]

area.log.exp.LPS <- ggplot(dls.titer.nophy.LPS, aes(x = value, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(title = "LPS", x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

dls.titer.nophy.OMKO <- dls.titer[(!dls.titer$experiment == 'physical') & (dls.titer$phage == 'OMKO'),]

area.log.exp.OMKO <- ggplot(dls.titer.nophy.OMKO, aes(x = value, y = -titer.diff)) +
  geom_point(aes(color = experiment)) + 
  labs(title = 'OMKO', x = "Area under the curve", y = "Difference in titer (log10)") + 
  theme_minimal() +
  geom_smooth(method='lm',formula= y ~ x -1)

plot.t <- ggarrange(area.log.exp.LPS, area.log.exp.OMKO)

ggsave('plot_tejas.tiff', device='tiff', dpi=400)

summary(lm(titer.diff ~ value, dls.titer.nophy.LPS))

summary(lm(titer.diff ~ value, dls.titer.nophy.OMKO))

## Plot phage control

#LPS
dls.control.lps <- dls[dls$phage == 'Luz19' & dls$measure == 'intensity',]

ggplot(dls.control.lps, aes(x = size, y = value, color = setting)) +
  geom_line() + 
  scale_x_continuous(trans='log10') 



#T4P




### test OMKO stuff

dls.physical.omko <- dls[dls$experiment == 'physical' & dls$phage == 'OMKO' & dls$measure == 'intensity',]

ggplot(dls.physical.omko, aes(x = size, y = value, color = setting)) +
  geom_line() + 
  scale_x_continuous(trans='log10') 
  

dls.control.omko <- dls[dls$phage == 'T4P' & dls$measure == 'intensity',]

ggplot(dls.control.omko, aes(x = size, y = value, color = setting)) +
  geom_line() + 
  scale_x_continuous(trans='log10') 

tapply(dls.control.omko$value, dls.control.omko$setting, sum)
