# Load packages
# Note: you will need to install these packages first if it is your first
# time using them
library('tidyr') 

# Set working directory --> change path to where you saved the csv files
setwd("~/Documents/Stanford/Research /DLS project") 

# Load data
dls_wide <- read.csv('DLS_raw_data.csv')
titer <- read.csv('titer_raw_data.csv')
phago_wide <- read.csv('PhagoBurn_DLS.csv')

##### Data processing #####

# In general, it is easier to work with long, rather than wide data. 
# Instead of having a lot of columns, we will have more rows. This makes 
# some operations easier, and also allow you to group your phages and treatments
# in different categories. We will use a R function to do this. 

dls_wide <- dls_wide[complete.cases(dls_wide),]
dls <- gather(dls_wide, key = 'setting', value = 'value', T2:OMKO_20_BME)

phago_wide <- phago_wide[complete.cases(phago_wide),]
phago <- gather(phago_wide, key = 'setting', value = 'value', LUZ14.Day.0:PEV2.Day.68)

### Add different columns that can be used to subset the data later

## DLS ##

# Add a 'experiment' column
dls$experiment <- NaN

dls$experiment[dls$setting %in% c('T2','T3','T4','T6','T7')] <- 'time'

dls$experiment[dls$setting %in% c("OMKO_control_lyo","LPS_control_lyo",
                                  "T4P_control_lyo","OMKO_lyo","LPS_lyo","T4P_lyo")] <- 'lyo'

dls$experiment[dls$setting %in% c("LUZ19_control_salt", "LUZ19_highsalt","LUZ19_lowsalt",
                                  "LUZ24_control_salt", "LUZ24_highsalt","LUZ24_lowsalt")] <- 'salt'

dls$experiment[dls$setting %in% c("LPS_control_peroxide","LPS_0.0002","LPS_0.002","LPS_0.02","LPS_0.2",
                                  "LPS_2","LPS_20", "OMKO_control_per","OMKO_0.002","OMKO_0.02",
                                  "OMKO_0.2","OMKO_2","OMKO_20")] <- 'peroxide' 

dls$experiment[dls$setting %in% c("OMKO_control_sonic","LPS_control_sonic","T4P_control_sonic",
                                  "LUZ19_control_sonic","LUZ24_control_sonic","OMKO_sonic",
                                  "LPS_sonic","T4P_sonic","LUZ19_sonic","LUZ24_sonic",
                                  "OMKO_vortex","LPS_vortex", "T4P_vortex","LUZ19_vortex","LUZ24_vortex")] <- 'physical'

dls$experiment[dls$setting %in% c("LPS_control_red","LPS_20_ascorbate","LPS_20_cysteine","LPS_20_glutathione",
                                  "LPS_20_DTT","LPS_20_TCEP","LPS_20_BME","OMKO_control_red","OMKO_20_ascorbate",
                                  "OMKO_20_cysteine","OMKO_20_glutathione","OMKO_20_DTT","OMKO_20_TCEP",
                                  "OMKO_20_BME")] <- 'reducing'


sum(dls$experiment == "NaN")

phago$day <- sapply(strsplit(phago$setting,'[.]'), getElement, 3)

# Add a 'phage' column
dls$phage <- sapply(strsplit(dls$setting, "_"), getElement, 1)

sum(dls$phage == "NaN")

phago$phage <- sapply(strsplit(phago$setting,'[.]'), getElement, 1)

# Add a 'group' column
dls$group <- NaN
dls$group[dls$setting %in% c("OMKO_control_lyo", "LPS_control_lyo","T4P_control_lyo","LUZ24_control_salt","LUZ19_control_salt",
                             "LPS_control_peroxide","OMKO_control_per","OMKO_control_sonic","LPS_control_sonic",
                             "T4P_control_sonic", "LUZ19_control_sonic","LUZ24_control_sonic","LPS_control_red",
                             "OMKO_control_red")] <- 'control'

dls$group[dls$setting %in% c("T2","T3","T4","T6","T7","OMKO_lyo","LPS_lyo","T4P_lyo","LUZ19_highsalt","LUZ19_lowsalt",
                             "LUZ24_highsalt","LUZ24_lowsalt","LPS_0.0002","LPS_0.002","LPS_0.02","LPS_0.2","LPS_2","LPS_20",
                             "OMKO_0.002","OMKO_0.02","OMKO_0.2","OMKO_2","OMKO_20","OMKO_sonic","LPS_sonic","T4P_sonic", 
                             "LUZ19_sonic", "LUZ24_sonic","OMKO_vortex","LPS_vortex","T4P_vortex","LUZ19_vortex",
                             "LUZ24_vortex","LPS_20_ascorbate","LPS_20_cysteine","LPS_20_glutathione",
                             "LPS_20_DTT","LPS_20_TCEP","LPS_20_BME","OMKO_20_ascorbate","OMKO_20_cysteine","OMKO_20_glutathione","OMKO_20_DTT","OMKO_20_TCEP",
                             "OMKO_20_BME")] <- 'treatment'

sum(dls$group == "NaN")

## Titer

titer$experiment <- NaN
titer$experiment[titer$setting %in% c('T2','T3','T4','T6','T7')] <- 'time'
titer$experiment[titer$setting %in% c("OMKO_control_lyo","LPS_control_lyo",
                                      "T4P_control_lyo","OMKO_lyo","LPS_lyo","T4P_lyo")] <- 'lyo'

titer$experiment[titer$setting %in% c("LUZ19_control_salt", "LUZ19_highsalt","LUZ19_lowsalt",
                                      "LUZ24_control_salt", "LUZ24_highsalt","LUZ24_lowsalt")] <- 'salt'

titer$experiment[titer$setting %in% c("LPS_control_peroxide","LPS_0.0002","LPS_0.002","LPS_0.02","LPS_0.2",
                                      "LPS_2","LPS_20", "OMKO_control_peroxide","OMKO_0.002","OMKO_0.02",
                                      "OMKO_0.2","OMKO_2","OMKO_20")] <- 'peroxide' 

titer$experiment[titer$setting %in% c("OMKO_control_sonic","LPS_control_sonic","T4P_control_sonic",
                                      "LUZ19_control_sonic","LUZ24_control_sonic","OMKO_sonic",
                                      "LPS_sonic","T4P_sonic","LUZ19_sonic","LUZ24_sonic",
                                      "OMKO_vortex","LPS_vortex", "T4P_vortex","LUZ19_vortex","LUZ24_vortex")] <- 'physical'

titer$experiment[titer$setting %in% c("LPS_control_red","LPS_20_ascorbate","LPS_20_cysteine","LPS_20_glutathione",
                                  "LPS_20_DTT","LPS_20_TCEP","LPS_20_BME","OMKO_control_red","OMKO_20_ascorbate",
                                  "OMKO_20_cysteine","OMKO_20_glutathione","OMKO_20_DTT","OMKO_20_TCEP",
                                  "OMKO_20_BME")] <- 'reducing'

sum(titer$experiment == "NaN")

# Add a 'phage' column
titer$phage <- sapply(strsplit(titer$setting, "_"), getElement, 1)

sum(titer$phage == "NaN")

# Add a 'group' column
titer$group <- NaN
titer$group[titer$setting %in% c("OMKO_control_lyo", "LPS_control_lyo","T4P_control_lyo","LUZ24_control_salt","LUZ19_control_salt",
                                 "LPS_control_peroxide","OMKO_control_peroxide","OMKO_control_sonic","LPS_control_sonic",
                                 "T4P_control_sonic", "LUZ19_control_sonic","LUZ24_control_sonic","LPS_control_red",
                                 "OMKO_control_red")] <- 'control'

titer$group[titer$setting %in% c("T2","T3","T4","T6","T7","OMKO_lyo","LPS_lyo","T4P_lyo","LUZ19_highsalt","LUZ19_lowsalt",
                             "LUZ24_highsalt","LUZ24_lowsalt","LPS_0.0002","LPS_0.002","LPS_0.02","LPS_0.2","LPS_2","LPS_20",
                             "OMKO_0.002","OMKO_0.02","OMKO_0.2","OMKO_2","OMKO_20","OMKO_sonic","LPS_sonic","T4P_sonic", 
                             "LUZ19_sonic", "LUZ24_sonic","OMKO_vortex","LPS_vortex","T4P_vortex","LUZ19_vortex",
                             "LUZ24_vortex","LPS_20_ascorbate","LPS_20_cysteine","LPS_20_glutathione",
                             "LPS_20_DTT","LPS_20_TCEP","LPS_20_BME","OMKO_20_ascorbate","OMKO_20_cysteine","OMKO_20_glutathione","OMKO_20_DTT","OMKO_20_TCEP",
                             "OMKO_20_BME")] <- 'treatment'

sum(titer$group == "NaN")

### Save processed dataset

dls.time <- dls[dls$experiment == "time",]
dls <- dls[dls$experiment != "time",]

titer.time <- titer[titer$experiment == "time",]
titer <- titer[titer$experiment != 'time',]

write.csv(dls.time, 'dls_time.csv')
write.csv(dls, 'dls.csv')

write.csv(titer.time, 'titer_time.csv')
write.csv(titer, 'titer.csv')

write.csv(phago, 'phago_processed.csv')

