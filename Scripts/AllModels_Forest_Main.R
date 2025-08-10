## Master Script ##

##### VC01 #####

# Install and load packages
#install.packages("pacman")
library(pacman)
p_load("devtools", "esc", "metafor", "brms", "brmstools", "ggplot2", 
       "dplyr", "bayesplot", "robvis", "RobustBayesianCopas", "stringr", "stringi", "tidybayes", "forcats", "ggridges")

# Read in effect size estimates for Variable Category x
Relatives <- read.csv("VC01_RelativesInGroup_All.csv")

# Remove punctuation in study names
Relatives$study <- str_replace(Relatives$study, ",", "")   
Relatives$study <- str_replace(Relatives$study, "_", " ")
Relatives$study <- stri_replace_all_fixed(Relatives$study, "(", "")
Relatives$study <- stri_replace_all_fixed(Relatives$study, ")", "")

# Inspect data
glimpse(Relatives)

### Calculate summary effect measures ###
# Assemble log odds estimates
Relatives$yi <- as.numeric(ifelse(Relatives$effect_measure == 1 | 
                                    Relatives$effect_measure == 8,
                                  Relatives$effect_size, "NA"))

# Convert odds ratios to log odds
Relatives$yi <- ifelse(Relatives$effect_measure == 8, 
                       log(Relatives$yi),
                       Relatives$yi)

# Recode log odds to align estimates
Relatives$yi <- ifelse(Relatives$recode == 1, 
                       Relatives$yi * -1, 
                       Relatives$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
Relatives <- Relatives[!is.na(Relatives$yi), ]

#Relabel SE column in dataset as sei
names(Relatives)[names(Relatives) == "SE"] <- "sei"

#Convert to numeric
Relatives$sei <- as.numeric(Relatives$sei)

#Impute missing standard errors...
#From t values
Relatives$sei <- ifelse(is.na(Relatives$sei) & !is.na(Relatives$t_value),
                        Relatives$effect_size/Relatives$t_value,
                        Relatives$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
Relatives$p_value <- as.numeric(dplyr::recode(Relatives$p_value, "<0.01" = "0.0075", 
                                              "<0.05" = "0.025", ">0.05" = "0.075", 
                                              "<0.1" = "0.075", ">0.1" = "0.125",
                                              "<0.15" = "0.125"))

#Impute z values
Relatives$z_value <- ifelse(is.na(Relatives$sei) & 
                              is.na(Relatives$t_value) &
                              is.na(Relatives$z_value),
                            abs(qnorm(Relatives$p_value/2)),
                            Relatives$z_value)

#Impute standard errors from z values
Relatives$sei <- ifelse(is.na(Relatives$sei) & !is.na(Relatives$z_value),
                        Relatives$effect_size/Relatives$z_value,
                        Relatives$sei)

# Log the standard errors of odds ratios
Relatives$sei <- ifelse(Relatives$effect_measure == 8, 
                        log(Relatives$sei),
                        Relatives$sei)

#Ensure standard errors are non-negative
Relatives$sei <- abs(Relatives$sei)

### Meta-analytic models ###

# Bayesian random effects model (continuous predictors with low and medium risk of bias)
rem.brms.Relatives <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                          prior = c(prior(normal(0, 1), class = Intercept),
                                    prior(cauchy(0, 1), class = sd)),
                          data = subset(Relatives, risk_of_bias == "low" | 
                                          risk_of_bias == "medium" &
                                          predictor_measure == 2),
                          iter = 50000,
                          cores = 4,
                          seed = 42,
                          control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC02 #####

# Read in effect size estimates for Variable Category x
PrAcq <- read.csv("VC02_PriorAcquaintance_All.csv")

# Remove punctuation in study names
PrAcq$study <- str_replace(PrAcq$study, ",", "")   
PrAcq$study <- str_replace(PrAcq$study, "_", " ")
PrAcq$study <- stri_replace_all_fixed(PrAcq$study, "(", "")
PrAcq$study <- stri_replace_all_fixed(PrAcq$study, ")", "")

# Inspect data
glimpse(PrAcq)

### Calculate summary effect measures ###
# B: Assemble log odds estimates
PrAcq$yi <- as.numeric(ifelse(PrAcq$effect_measure == 1 | 
                                PrAcq$effect_measure == 8,
                              PrAcq$effect_size, "NA"))

# Convert odds ratios to log odds
PrAcq$yi <- ifelse(PrAcq$effect_measure == 8, 
                   log(PrAcq$yi),
                   PrAcq$yi)

# Recode log odds to align estimates
PrAcq$yi <- ifelse(PrAcq$recode == 1, 
                   PrAcq$yi * -1, 
                   PrAcq$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
PrAcq <- PrAcq[!is.na(PrAcq$yi), ]

#Relabel SE column in dataset as sei
names(PrAcq)[names(PrAcq) == "SE"] <- "sei"

#Convert to numeric
PrAcq$sei <- as.numeric(PrAcq$sei)

#Impute missing standard errors...
#From t values
PrAcq$sei <- ifelse(is.na(PrAcq$sei) & !is.na(PrAcq$t_value),
                    PrAcq$effect_size/PrAcq$t_value,
                    PrAcq$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
PrAcq$p_value <- as.numeric(dplyr::recode(PrAcq$p_value, "<0.01" = "0.0075", 
                                          "<0.05" = "0.025", ">0.05" = "0.075", 
                                          "<0.1" = "0.075", ">0.1" = "0.125",
                                          "<0.15" = "0.125"))

#Impute z values
PrAcq$z_value <- ifelse(is.na(PrAcq$sei) & 
                          is.na(PrAcq$t_value) &
                          is.na(PrAcq$z_value),
                        abs(qnorm(PrAcq$p_value/2)),
                        PrAcq$z_value)

#Impute standard errors from z values
PrAcq$sei <- ifelse(is.na(PrAcq$sei) & !is.na(PrAcq$z_value),
                    PrAcq$effect_size/PrAcq$z_value,
                    PrAcq$sei)

#Ensure standard errors are non-negative
PrAcq$sei <- abs(PrAcq$sei)

# Log the standard errors of odds ratios
PrAcq$sei <- ifelse(PrAcq$effect_measure == 8, 
                    log(PrAcq$sei),
                    PrAcq$sei)

#Restrict sample to studies with exchangeable predictors
PrAcq <- subset(PrAcq, predictor_measure == 1)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.PrAcq <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                      prior = c(prior(normal(0, 1), class = Intercept),
                                prior(cauchy(0, 1), class = sd)),
                      data = subset(PrAcq, risk_of_bias == "low" | 
                                      risk_of_bias == "medium"),
                      iter = 50000,
                      cores = 4,
                      seed = 42,
                      control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC03 #####
# Read in effect size estimates for Variable Category x
GroupTenure <- read.csv("VC03_GroupTenure_All.csv")

# Remove punctuation in study names
GroupTenure$study <- str_replace(GroupTenure$study, ",", "")   
GroupTenure$study <- str_replace(GroupTenure$study, "_", " ")
GroupTenure$study <- stri_replace_all_fixed(GroupTenure$study, "(", "")
GroupTenure$study <- stri_replace_all_fixed(GroupTenure$study, ")", "")

# Inspect data
glimpse(GroupTenure)

### Calculate summary effect measures ###
# B: Assemble log odds estimates
GroupTenure$yi <- as.numeric(ifelse(GroupTenure$effect_measure == 1 | 
                                      GroupTenure$effect_measure == 8,
                                    GroupTenure$effect_size, "NA"))

# Convert odds ratios to log odds
GroupTenure$yi <- ifelse(GroupTenure$effect_measure == 8, 
                         log(GroupTenure$yi),
                         GroupTenure$yi)

# Recode log odds to align estimates
GroupTenure$yi <- ifelse(GroupTenure$recode == 1, 
                         GroupTenure$yi * -1, 
                         GroupTenure$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
GroupTenure <- GroupTenure[!is.na(GroupTenure$yi), ]

#Relabel SE column in dataset as sei
names(GroupTenure)[names(GroupTenure) == "SE"] <- "sei"

#Convert to numeric
GroupTenure$sei <- as.numeric(GroupTenure$sei)

#Impute missing standard errors...
#From t values
GroupTenure$sei <- ifelse(is.na(GroupTenure$sei) & !is.na(GroupTenure$t_value),
                          GroupTenure$effect_size/GroupTenure$t_value,
                          GroupTenure$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
GroupTenure$p_value <- as.numeric(dplyr::recode(GroupTenure$p_value, "<0.01" = "0.0075", 
                                                "<0.05" = "0.025", ">0.05" = "0.075", 
                                                "<0.1" = "0.075", ">0.1" = "0.125",
                                                "<0.15" = "0.125"))

#Impute z values
GroupTenure$z_value <- ifelse(is.na(GroupTenure$sei) & 
                                is.na(GroupTenure$t_value) &
                                is.na(GroupTenure$z_value),
                              abs(qnorm(GroupTenure$p_value/2)),
                              GroupTenure$z_value)

#Impute standard errors from z values
GroupTenure$sei <- ifelse(is.na(GroupTenure$sei) & !is.na(GroupTenure$z_value),
                          GroupTenure$effect_size/GroupTenure$z_value,
                          GroupTenure$sei)

#Ensure standard errors are non-negative
GroupTenure$sei <- abs(GroupTenure$sei)

#If sei = 0, change to 0.001
GroupTenure$sei <- ifelse(GroupTenure$sei == 0, 
                          0.001, 
                          GroupTenure$sei)

# Log the standard errors of odds ratios
GroupTenure$sei <- ifelse(GroupTenure$effect_measure == 8, 
                          log(GroupTenure$sei),
                          GroupTenure$sei)

#Convert estimates from effect of group age (months) to group age (years)
GroupTenure$yi <- ifelse(GroupTenure$year_or_month == "m", 
                         GroupTenure$yi * 12, 
                         GroupTenure$yi)

GroupTenure$sei <- ifelse(GroupTenure$year_or_month == "m", 
                          GroupTenure$sei * 12, 
                          GroupTenure$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.GroupTenure <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                            prior = c(prior(normal(0, 1), class = Intercept),
                                      prior(cauchy(0, 1), class = sd)),
                            data = subset(GroupTenure, risk_of_bias == "low" | 
                                            risk_of_bias == "medium"),
                            iter = 50000,
                            cores = 4,
                            seed = 42,
                            control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC04 #####
#Insufficent data#

##### VC05 #####

# Read in effect size estimates for Variable Category x
GeogProx <- read.csv("VC05_GeographicProximity_All.csv")

# Remove punctuation in study names
GeogProx$study <- str_replace(GeogProx$study, ",", "")   
GeogProx$study <- str_replace(GeogProx$study, "_", " ")
GeogProx$study <- stri_replace_all_fixed(GeogProx$study, "(", "")
GeogProx$study <- stri_replace_all_fixed(GeogProx$study, ")", "")

# Inspect data
glimpse(GeogProx)

### Calculate summary effect measures ###
# B: Assemble log odds estimates
GeogProx$yi <- as.numeric(ifelse(GeogProx$effect_measure == 1 | 
                                   GeogProx$effect_measure == 8,
                                 GeogProx$effect_size, "NA"))

# Convert odds ratios to log odds
GeogProx$yi <- ifelse(GeogProx$effect_measure == 8, 
                      log(GeogProx$yi),
                      GeogProx$yi)

# Recode log odds to align estimates
GeogProx$yi <- ifelse(GeogProx$recode == 1, 
                      GeogProx$yi * -1, 
                      GeogProx$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
GeogProx <- GeogProx[!is.na(GeogProx$yi), ]

#Relabel SE column in dataset as sei
names(GeogProx)[names(GeogProx) == "SE"] <- "sei"

#Convert to numeric
GeogProx$sei <- as.numeric(GeogProx$sei)

#Impute missing standard errors...
#From t values
GeogProx$sei <- ifelse(is.na(GeogProx$sei) & !is.na(GeogProx$t_value),
                       GeogProx$effect_size/GeogProx$t_value,
                       GeogProx$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
GeogProx$p_value <- as.numeric(dplyr::recode(GeogProx$p_value, "<0.01" = "0.0075", 
                                             "<0.05" = "0.025", ">0.05" = "0.075", 
                                             "<0.1" = "0.075", ">0.1" = "0.125",
                                             "<0.15" = "0.125"))

#Impute z values
GeogProx$z_value <- ifelse(is.na(GeogProx$sei) & 
                             is.na(GeogProx$t_value) &
                             is.na(GeogProx$z_value),
                           abs(qnorm(GeogProx$p_value/2)),
                           GeogProx$z_value)

#Impute standard errors from z values
GeogProx$sei <- ifelse(is.na(GeogProx$sei) & !is.na(GeogProx$z_value),
                       GeogProx$effect_size/GeogProx$z_value,
                       GeogProx$sei)

# Log the standard errors of odds ratios
GeogProx$sei <- ifelse(GeogProx$effect_measure == 8, 
                       log(GeogProx$sei),
                       GeogProx$sei)

#Ensure standard errors are non-negative
GeogProx$sei <- abs(GeogProx$sei)

#Convert effect estimates from distance (m) to distance (km)
GeogProx$yi <- ifelse(GeogProx$m_or_km == "m", 
                      GeogProx$yi * 1000, 
                      GeogProx$yi)

GeogProx$sei <- ifelse(GeogProx$m_or_km == "m", 
                       GeogProx$sei * 1000, 
                       GeogProx$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.GeogProx <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                         prior = c(prior(normal(0, 1), class = Intercept),
                                   prior(cauchy(0, 1), class = sd)),
                         data = subset(GeogProx, risk_of_bias == "low" | 
                                         risk_of_bias == "medium"),
                         iter = 50000,
                         cores = 4,
                         seed = 42,
                         control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC06 #####
# Read in effect size estimates for Variable Category x
MemManage <- read.csv("VC06_MemberManagement_All.csv")

# Remove punctuation in study names
MemManage$study <- str_replace(MemManage$study, ",", "")   
MemManage$study <- str_replace(MemManage$study, "_", " ")
MemManage$study <- stri_replace_all_fixed(MemManage$study, "(", "")
MemManage$study <- stri_replace_all_fixed(MemManage$study, ")", "")

# Inspect data
glimpse(MemManage)

### Calculate summary effect measures ###
# B: Assemble log odds estimates
MemManage$yi <- as.numeric(ifelse(MemManage$effect_measure == 1 | 
                                    MemManage$effect_measure == 8,
                                  MemManage$effect_size, "NA"))

# Convert odds ratios to log odds
MemManage$yi <- ifelse(MemManage$effect_measure == 8, 
                       log(MemManage$yi),
                       MemManage$yi)

# Recode log odds to align estimates
MemManage$yi <- ifelse(MemManage$recode == 1, 
                       MemManage$yi * -1, 
                       MemManage$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
MemManage <- MemManage[!is.na(MemManage$yi), ]

#Relabel SE column in dataset as sei
names(MemManage)[names(MemManage) == "SE"] <- "sei"

#Convert to numeric
MemManage$sei <- as.numeric(MemManage$sei)

#Impute missing standard errors...
#From t values
MemManage$sei <- ifelse(is.na(MemManage$sei) & !is.na(MemManage$t_value),
                        MemManage$effect_size/MemManage$t_value,
                        MemManage$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
MemManage$p_value <- as.numeric(dplyr::recode(MemManage$p_value, "<0.01" = "0.0075", 
                                              "<0.05" = "0.025", ">0.05" = "0.075", 
                                              "<0.1" = "0.075", ">0.1" = "0.125",
                                              "<0.15" = "0.125"))

#Impute z values
MemManage$z_value <- ifelse(is.na(MemManage$sei) & 
                              is.na(MemManage$t_value) &
                              is.na(MemManage$z_value),
                            abs(qnorm(MemManage$p_value/2)),
                            MemManage$z_value)

#Impute standard errors from z values
MemManage$sei <- ifelse(is.na(MemManage$sei) & !is.na(MemManage$z_value),
                        MemManage$effect_size/MemManage$z_value,
                        MemManage$sei)

# Log the standard errors of odds ratios
MemManage$sei <- ifelse(MemManage$effect_measure == 8, 
                        log(MemManage$sei),
                        MemManage$sei)

#Ensure standard errors are non-negative
MemManage$sei <- abs(MemManage$sei)

#Round Kiros_2014 yi and sei (both 0) to 0.001
MemManage$yi <- ifelse(MemManage$model_ID == "Kiros_2014_M1", 
                       0.001, 
                       MemManage$yi)

MemManage$sei <- ifelse(MemManage$model_ID == "Kiros_2014_M1", 
                        0.001, 
                        MemManage$sei)

#Log the Noglo estimates
MemManage$yi <- ifelse(MemManage$model_ID == "NogloA_2015_01", 
                       -1 * log(16.473), 
                       MemManage$yi)

MemManage$sei <- ifelse(MemManage$model_ID == "NogloA_2015_01", 
                        log(MemManage$sei), 
                        MemManage$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.MemManage <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                          prior = c(prior(normal(0, 1), class = Intercept),
                                    prior(cauchy(0, 1), class = sd)),
                          data = subset(MemManage, risk_of_bias == "low" | 
                                          risk_of_bias == "medium"),
                          iter = 50000,
                          cores = 4,
                          seed = 42,
                          control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC07 #####

# Read in effect size estimates for Variable Category x
GroupSanctions <- read.csv("VC07_GroupSanctions_All.csv")

# Remove punctuation in study names
GroupSanctions$study <- str_replace(GroupSanctions$study, ",", "")   
GroupSanctions$study <- str_replace(GroupSanctions$study, "_", " ")
GroupSanctions$study <- stri_replace_all_fixed(GroupSanctions$study, "(", "")
GroupSanctions$study <- stri_replace_all_fixed(GroupSanctions$study, ")", "")

# Inspect data
glimpse(GroupSanctions)

### Calculate summary effect measures ###
# Assemble log odds estimates
GroupSanctions$yi <- as.numeric(ifelse(GroupSanctions$effect_measure == 1 | 
                                         GroupSanctions$effect_measure == 8,
                                       GroupSanctions$effect_size, "NA"))

# Convert odds ratios to log odds
GroupSanctions$yi <- ifelse(GroupSanctions$effect_measure == 8, 
                            log(GroupSanctions$yi),
                            GroupSanctions$yi)

# Recode log odds to align estimates
GroupSanctions$yi <- ifelse(GroupSanctions$recode == 1, 
                            GroupSanctions$yi * -1, 
                            GroupSanctions$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
GroupSanctions <- GroupSanctions[!is.na(GroupSanctions$yi), ]

#Relabel SE column in dataset as sei
names(GroupSanctions)[names(GroupSanctions) == "SE"] <- "sei"

#Convert to numeric
GroupSanctions$sei <- as.numeric(GroupSanctions$sei)

#Impute missing standard errors...
#From t values
GroupSanctions$sei <- ifelse(is.na(GroupSanctions$sei) & !is.na(GroupSanctions$t_value),
                             GroupSanctions$effect_size/GroupSanctions$t_value,
                             GroupSanctions$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
GroupSanctions$p_value <- as.numeric(dplyr::recode(GroupSanctions$p_value, "<0.01" = "0.0075", 
                                                   "<0.05" = "0.025", ">0.05" = "0.075", 
                                                   "<0.1" = "0.075", ">0.1" = "0.125",
                                                   "<0.15" = "0.125"))

#Impute z values
GroupSanctions$z_value <- ifelse(is.na(GroupSanctions$sei) & 
                                   is.na(GroupSanctions$t_value) &
                                   is.na(GroupSanctions$z_value),
                                 abs(qnorm(GroupSanctions$p_value/2)),
                                 GroupSanctions$z_value)

#Impute standard errors from z values
GroupSanctions$sei <- ifelse(is.na(GroupSanctions$sei) & !is.na(GroupSanctions$z_value),
                             GroupSanctions$effect_size/GroupSanctions$z_value,
                             GroupSanctions$sei)

# Log the standard errors of odds ratios
GroupSanctions$sei <- ifelse(GroupSanctions$effect_measure == 8, 
                             log(GroupSanctions$sei),
                             GroupSanctions$sei)

#Log the Noglo estimates
GroupSanctions$yi <- ifelse(GroupSanctions$study == "Noglo and Androuais 2015", 
                            log(GroupSanctions$yi), 
                            GroupSanctions$yi)

GroupSanctions$sei <- ifelse(GroupSanctions$study == "Noglo and Androuais 2015", 
                             -1*(log(GroupSanctions$sei)), 
                             GroupSanctions$sei)

#Ensure standard errors are non-negative
GroupSanctions$sei <- abs(GroupSanctions$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.GroupSanctions <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                               prior = c(prior(normal(0, 1), class = Intercept),
                                         prior(cauchy(0, 1), class = sd)),
                               data = subset(GroupSanctions, risk_of_bias == "low" | 
                                               risk_of_bias == "medium"),
                               iter = 50000,
                               cores = 4,
                               seed = 42,
                               control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC08 #####

# Read in effect size estimates for Variable Category x
PeerMon <- read.csv("VC08_PeerMonitoring_All.csv")

# Remove punctuation in study names
PeerMon$study <- str_replace(PeerMon$study, ",", "")   
PeerMon$study <- str_replace(PeerMon$study, "_", " ")
PeerMon$study <- stri_replace_all_fixed(PeerMon$study, "(", "")
PeerMon$study <- stri_replace_all_fixed(PeerMon$study, ")", "")

# Inspect data
glimpse(PeerMon)

### Calculate summary effect measures ###
# Select comparable studies
PeerMon$effect_size <- as.numeric(ifelse(PeerMon$predictor_measure == 1,
                                         PeerMon$effect_size, "NA"))

# Assemble log odds estimates
PeerMon$yi <- as.numeric(ifelse(PeerMon$effect_measure == 1 | 
                                  PeerMon$effect_measure == 8,
                                PeerMon$effect_size, "NA"))

# Convert odds ratios to log odds
PeerMon$yi <- ifelse(PeerMon$effect_measure == 8, 
                     log(PeerMon$yi),
                     PeerMon$yi)

# Recode log odds to align estimates
PeerMon$yi <- ifelse(PeerMon$recode == 1, 
                     PeerMon$yi * -1, 
                     PeerMon$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
PeerMon <- PeerMon[!is.na(PeerMon$yi), ]

#Relabel SE column in dataset as sei
names(PeerMon)[names(PeerMon) == "SE"] <- "sei"

#Convert to numeric
PeerMon$sei <- as.numeric(PeerMon$sei)

#Impute missing standard errors...
#From t values
PeerMon$sei <- ifelse(is.na(PeerMon$sei) & !is.na(PeerMon$t_value),
                      PeerMon$effect_size/PeerMon$t_value,
                      PeerMon$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
PeerMon$p_value <- as.numeric(dplyr::recode(PeerMon$p_value, "<0.01" = "0.0075", 
                                            "<0.05" = "0.025", ">0.05" = "0.075", 
                                            "<0.1" = "0.075", ">0.1" = "0.125",
                                            "<0.15" = "0.125"))

#Impute z values
PeerMon$z_value <- ifelse(is.na(PeerMon$sei) & 
                            is.na(PeerMon$t_value) &
                            is.na(PeerMon$z_value),
                          abs(qnorm(PeerMon$p_value/2)),
                          PeerMon$z_value)

#Impute standard errors from z values
PeerMon$sei <- ifelse(is.na(PeerMon$sei) & !is.na(PeerMon$z_value),
                      PeerMon$effect_size/PeerMon$z_value,
                      PeerMon$sei)

# Log the standard errors of odds ratios
PeerMon$sei <- ifelse(PeerMon$effect_measure == 8, 
                      log(PeerMon$sei),
                      PeerMon$sei)

#Log the Noglo estimates
PeerMon$yi <- ifelse(PeerMon$study == "Noglo and Androuais, 2015", 
                     log(PeerMon$yi), 
                     PeerMon$yi)

PeerMon$sei <- ifelse(PeerMon$study == "Noglo and Androuais, 2015", 
                      -1*(log(PeerMon$sei)), 
                      PeerMon$sei)

#Ensure standard errors are non-negative
PeerMon$sei <- abs(PeerMon$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.PeerMon <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(cauchy(0, 1), class = sd)),
                        data = subset(PeerMon, risk_of_bias == "low" | 
                                        risk_of_bias == "medium"),
                        iter = 50000,
                        cores = 4,
                        seed = 42,
                        control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC09 #####

# Read in effect size estimates for Variable Category x
ExtMonitoring <- read.csv("VC09_ExternalMonitoring_All.csv")

# Remove punctuation in study names
ExtMonitoring$study <- str_replace(ExtMonitoring$study, ",", "")   
ExtMonitoring$study <- str_replace(ExtMonitoring$study, "_", " ")
ExtMonitoring$study <- stri_replace_all_fixed(ExtMonitoring$study, "(", "")
ExtMonitoring$study <- stri_replace_all_fixed(ExtMonitoring$study, ")", "")

# Inspect data
glimpse(ExtMonitoring)

### Calculate summary effect measures ###
# Assemble log odds estimates
ExtMonitoring$yi <- as.numeric(ifelse(ExtMonitoring$effect_measure == 1 | 
                                        ExtMonitoring$effect_measure == 8,
                                      ExtMonitoring$effect_size, "NA"))

# Convert odds ratios to log odds
ExtMonitoring$yi <- ifelse(ExtMonitoring$effect_measure == 8, 
                           log(ExtMonitoring$yi),
                           ExtMonitoring$yi)

# Recode log odds to align estimates
ExtMonitoring$yi <- ifelse(ExtMonitoring$recode == 1, 
                           ExtMonitoring$yi * -1, 
                           ExtMonitoring$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
ExtMonitoring <- ExtMonitoring[!is.na(ExtMonitoring$yi), ]

#Relabel SE column in dataset as sei
names(ExtMonitoring)[names(ExtMonitoring) == "SE"] <- "sei"

#Convert to numeric
ExtMonitoring$sei <- as.numeric(ExtMonitoring$sei)

#Impute missing standard errors...
#From t values
ExtMonitoring$sei <- ifelse(is.na(ExtMonitoring$sei) & !is.na(ExtMonitoring$t_value),
                            ExtMonitoring$effect_size/ExtMonitoring$t_value,
                            ExtMonitoring$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
ExtMonitoring$p_value <- as.numeric(dplyr::recode(ExtMonitoring$p_value, "<0.01" = "0.0075", 
                                                  "<0.05" = "0.025", ">0.05" = "0.075", 
                                                  "<0.1" = "0.075", ">0.1" = "0.125",
                                                  "<0.15" = "0.125"))

#Impute z values
ExtMonitoring$z_value <- ifelse(is.na(ExtMonitoring$sei) & 
                                  is.na(ExtMonitoring$t_value) &
                                  is.na(ExtMonitoring$z_value),
                                abs(qnorm(ExtMonitoring$p_value/2)),
                                ExtMonitoring$z_value)

#Impute standard errors from z values
ExtMonitoring$sei <- ifelse(is.na(ExtMonitoring$sei) & !is.na(ExtMonitoring$z_value),
                            ExtMonitoring$effect_size/ExtMonitoring$z_value,
                            ExtMonitoring$sei)

# Log the standard errors of odds ratios
ExtMonitoring$sei <- ifelse(ExtMonitoring$effect_measure == 8, 
                            log(ExtMonitoring$sei),
                            ExtMonitoring$sei)

#Ensure standard errors are non-negative
ExtMonitoring$sei <- abs(ExtMonitoring$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.ExtMonitoring <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                              prior = c(prior(normal(0, 1), class = Intercept),
                                        prior(cauchy(0, 1), class = sd)),
                              data = subset(ExtMonitoring, risk_of_bias == "low" | 
                                              risk_of_bias == "medium"),
                              iter = 50000,
                              cores = 4,
                              seed = 42,
                              control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC010 #####
# Read in effect size estimates for Variable Category x
GroupSize <- read.csv("VC10_GroupSize_All.csv")

# Remove punctuation in study names
GroupSize$study <- str_replace(GroupSize$study, ",", "")   
GroupSize$study <- str_replace(GroupSize$study, "_", " ")
GroupSize$study <- stri_replace_all_fixed(GroupSize$study, "(", "")
GroupSize$study <- stri_replace_all_fixed(GroupSize$study, ")", "")

# Inspect data
glimpse(GroupSize)

### Calculate summary effect measures ###
# Assemble log odds estimates
GroupSize$yi <- as.numeric(ifelse(GroupSize$effect_measure == 1 | 
                                    GroupSize$effect_measure == 8,
                                  GroupSize$effect_size, "NA"))

# Convert odds ratios to log odds
GroupSize$yi <- ifelse(GroupSize$effect_measure == 8, 
                       log(GroupSize$yi),
                       GroupSize$yi)

# Recode log odds to align estimates
GroupSize$yi <- ifelse(GroupSize$recode == 1, 
                       GroupSize$yi * -1, 
                       GroupSize$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
GroupSize <- GroupSize[!is.na(GroupSize$yi), ]

#Relabel SE column in dataset as sei
names(GroupSize)[names(GroupSize) == "SE"] <- "sei"

#Convert to numeric
GroupSize$sei <- as.numeric(GroupSize$sei)

#Impute missing standard errors...
#From t values
GroupSize$sei <- ifelse(is.na(GroupSize$sei) & !is.na(GroupSize$t_value),
                        GroupSize$effect_size/GroupSize$t_value,
                        GroupSize$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
GroupSize$p_value <- as.numeric(dplyr::recode(GroupSize$p_value, "<0.01" = "0.0075", 
                                              "<0.05" = "0.025", ">0.05" = "0.075", 
                                              "<0.1" = "0.075", ">0.1" = "0.125",
                                              "<0.15" = "0.125"))

#Impute z values
GroupSize$z_value <- ifelse(is.na(GroupSize$sei) & 
                              is.na(GroupSize$t_value) &
                              is.na(GroupSize$z_value),
                            abs(qnorm(GroupSize$p_value/2)),
                            GroupSize$z_value)

#Impute standard errors from z values
GroupSize$sei <- ifelse(is.na(GroupSize$sei) & !is.na(GroupSize$z_value),
                        GroupSize$effect_size/GroupSize$z_value,
                        GroupSize$sei)

# Log the standard errors of odds ratios
GroupSize$sei <- ifelse(GroupSize$effect_measure == 8, 
                        log(GroupSize$sei),
                        GroupSize$sei)

#Log the Noglo estimates
#Make estimate positive before logging
GroupSize$yi <- ifelse(GroupSize$study == "Noglo and Androuais 2015", 
                       abs(GroupSize$yi), 
                       GroupSize$yi)

#Log estimate
GroupSize$yi <- ifelse(GroupSize$study == "Noglo and Androuais 2015", 
                       log(GroupSize$yi), 
                       GroupSize$yi)

#Return estimate to negative
GroupSize$yi <- ifelse(GroupSize$study == "Noglo and Androuais 2015", 
                       -1*(GroupSize$yi), 
                       GroupSize$yi)

#Log standard error
GroupSize$sei <- ifelse(GroupSize$study == "Noglo and Androuais 2015", 
                        -1*(log(GroupSize$sei)), 
                        GroupSize$sei)

#Ensure standard errors are non-negative
GroupSize$sei <- abs(GroupSize$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.GroupSize <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                          prior = c(prior(normal(0, 1), class = Intercept),
                                    prior(cauchy(0, 1), class = sd)),
                          data = subset(GroupSize, risk_of_bias == "low" | 
                                          risk_of_bias == "medium"),
                          iter = 50000,
                          cores = 4,
                          seed = 42,
                          control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC011 #####
# Read in effect size estimates for Variable Category x
BorrowerAge <- read.csv("VC11_BorrowerAge_All.csv")

# Remove punctuation in study names
BorrowerAge$study <- str_replace(BorrowerAge$study, ",", "")   
BorrowerAge$study <- str_replace(BorrowerAge$study, "_", " ")
BorrowerAge$study <- stri_replace_all_fixed(BorrowerAge$study, "(", "")
BorrowerAge$study <- stri_replace_all_fixed(BorrowerAge$study, ")", "")

# Inspect data
glimpse(BorrowerAge)

### Calculate summary effect measures ###
# Assemble log odds estimates
BorrowerAge$yi <- as.numeric(ifelse(BorrowerAge$effect_measure == 1 | 
                                      BorrowerAge$effect_measure == 8,
                                    BorrowerAge$effect_size, "NA"))

# Convert odds ratios to log odds
BorrowerAge$yi <- ifelse(BorrowerAge$effect_measure == 8, 
                         log(BorrowerAge$yi),
                         BorrowerAge$yi)

# Recode log odds to align estimates
BorrowerAge$yi <- ifelse(BorrowerAge$recode == 1, 
                         BorrowerAge$yi * -1, 
                         BorrowerAge$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
BorrowerAge <- BorrowerAge[!is.na(BorrowerAge$yi), ]

#Relabel SE column in dataset as sei
names(BorrowerAge)[names(BorrowerAge) == "SE"] <- "sei"

#Convert to numeric
BorrowerAge$sei <- as.numeric(BorrowerAge$sei)

#Impute missing standard errors...
#From t values
BorrowerAge$sei <- ifelse(is.na(BorrowerAge$sei) & !is.na(BorrowerAge$t_value),
                          BorrowerAge$effect_size/BorrowerAge$t_value,
                          BorrowerAge$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
BorrowerAge$p_value <- as.numeric(dplyr::recode(BorrowerAge$p_value, "<0.01" = "0.0075", 
                                                "<0.05" = "0.025", ">0.05" = "0.075", 
                                                "<0.1" = "0.075", ">0.1" = "0.125",
                                                "<0.15" = "0.125"))

#Impute z values
BorrowerAge$z_value <- ifelse(is.na(BorrowerAge$sei) & 
                                is.na(BorrowerAge$t_value) &
                                is.na(BorrowerAge$z_value),
                              abs(qnorm(BorrowerAge$p_value/2)),
                              BorrowerAge$z_value)

#Impute standard errors from z values
BorrowerAge$sei <- ifelse(is.na(BorrowerAge$sei) & !is.na(BorrowerAge$z_value),
                          BorrowerAge$effect_size/BorrowerAge$z_value,
                          BorrowerAge$sei)

# Log the standard errors of odds ratios
BorrowerAge$sei <- ifelse(BorrowerAge$effect_measure == 8, 
                          log(BorrowerAge$sei),
                          BorrowerAge$sei)

#Ensure standard errors are non-negative
BorrowerAge$sei <- abs(BorrowerAge$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.BorrowerAge <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                            prior = c(prior(normal(0, 1), class = Intercept),
                                      prior(cauchy(0, 1), class = sd)),
                            data = subset(BorrowerAge, risk_of_bias == "low" | 
                                            risk_of_bias == "medium"),
                            iter = 50000,
                            cores = 4,
                            seed = 42,
                            control = list(adapt_delta = 0.999, max_treedepth = 12))

##### VC012 #####
# Read in effect size estimates for Variable Category x
BorrowerSex <- read.csv("VC12_BorrowerSex_All.csv")

# Remove punctuation in study names
BorrowerSex$study <- str_replace(BorrowerSex$study, ",", "")   
BorrowerSex$study <- str_replace(BorrowerSex$study, "_", " ")
BorrowerSex$study <- stri_replace_all_fixed(BorrowerSex$study, "(", "")
BorrowerSex$study <- stri_replace_all_fixed(BorrowerSex$study, ")", "")

# Inspect data
glimpse(BorrowerSex)

### Calculate summary effect measures ###
# Assemble log odds estimates
BorrowerSex$yi <- as.numeric(ifelse(BorrowerSex$effect_measure == 1 | 
                                      BorrowerSex$effect_measure == 8,
                                    BorrowerSex$effect_size, "NA"))

# Convert odds ratios to log odds
BorrowerSex$yi <- ifelse(BorrowerSex$effect_measure == 8, 
                         log(BorrowerSex$yi),
                         BorrowerSex$yi)

# Recode log odds to align estimates
BorrowerSex$yi <- ifelse(BorrowerSex$recode == 1, 
                         BorrowerSex$yi * -1, 
                         BorrowerSex$yi)

### Impute missing standard errors ###
# Remove studies without relevant estimates
BorrowerSex <- BorrowerSex[!is.na(BorrowerSex$yi), ]

#Relabel SE column in dataset as sei
names(BorrowerSex)[names(BorrowerSex) == "SE"] <- "sei"

#Convert to numeric
BorrowerSex$sei <- as.numeric(BorrowerSex$sei)

#Impute missing standard errors...
#From t values
BorrowerSex$sei <- ifelse(is.na(BorrowerSex$sei) & !is.na(BorrowerSex$t_value),
                          BorrowerSex$effect_size/BorrowerSex$t_value,
                          BorrowerSex$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
BorrowerSex$p_value <- as.numeric(dplyr::recode(BorrowerSex$p_value, "<0.01" = "0.0075", 
                                                "<0.05" = "0.025", ">0.05" = "0.075", 
                                                "<0.1" = "0.075", ">0.1" = "0.125",
                                                "<0.15" = "0.125"))

#Impute z values
BorrowerSex$z_value <- ifelse(is.na(BorrowerSex$sei) & 
                                is.na(BorrowerSex$t_value) &
                                is.na(BorrowerSex$z_value),
                              abs(qnorm(BorrowerSex$p_value/2)),
                              BorrowerSex$z_value)

#Impute standard errors from z values
BorrowerSex$sei <- ifelse(is.na(BorrowerSex$sei) & !is.na(BorrowerSex$z_value),
                          BorrowerSex$effect_size/BorrowerSex$z_value,
                          BorrowerSex$sei)

# Log the standard errors of odds ratios
BorrowerSex$sei <- ifelse(BorrowerSex$effect_measure == 8, 
                          log(BorrowerSex$sei),
                          BorrowerSex$sei)

#Ensure standard errors are non-negative
BorrowerSex$sei <- abs(BorrowerSex$sei)

### Meta-analytic models ###

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms.BorrowerSex <- brm(yi | se(sei) ~ 1 + (1 | model_ID),
                            prior = c(prior(normal(0, 1), class = Intercept),
                                      prior(cauchy(0, 1), class = sd)),
                            data = subset(BorrowerSex, risk_of_bias == "low" | 
                                            risk_of_bias == "medium"),
                            iter = 50000,
                            cores = 4,
                            seed = 42,
                            control = list(adapt_delta = 0.999, max_treedepth = 12))

##### Meta-analysis forest plots #####

##### VC01 Relatives in Group ##### 
# Study-specific effects are deviations + average
out_r_Relatives <- spread_draws(rem.brms.Relatives, r_model_ID[model_ID], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept)

# Average effect
out_f_Relatives <- spread_draws(rem.brms.Relatives, b_Intercept) %>% 
  mutate(model_ID = "Relatives in Group Effect") 

# Combine average and study-specific effects' data frames
out_all_Relatives <- bind_rows(out_r_Relatives, out_f_Relatives) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Relatives in Group Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept)) 

# Data frame of summary numbers
out_all_Relatives_sum <- group_by(out_all_Relatives, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_Relatives %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_Relatives_sum, size = 1) +
  geom_text(data = mutate_if(out_all_Relatives_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 6), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = Relatives %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC02 Prior Acquaintance ##### 
# Study-specific effects are deviations + average
out_r_PrAcq <- spread_draws(rem.brms.PrAcq, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_PrAcq <- spread_draws(rem.brms.PrAcq, b_Intercept) %>% 
  mutate(model_ID = "Prior Acquaintance Effect") 

# Combine average and study-specific effects' data frames
out_all_PrAcq <- bind_rows(out_r_PrAcq, out_f_PrAcq) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Prior Acquaintance Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_PrAcq_sum <- group_by(out_all_PrAcq, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_PrAcq %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_PrAcq_sum, size = 1) +
  geom_text(data = mutate_if(out_all_PrAcq_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 9), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = PrAcq %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC03 Group Tenure #####
# Study-specific effects are deviations + average
out_r_GroupTenure <- spread_draws(rem.brms.GroupTenure, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_GroupTenure <- spread_draws(rem.brms.GroupTenure, b_Intercept) %>% 
  mutate(model_ID = "Group Tenure Effect") 

# Combine average and study-specific effects' data frames
out_all_GroupTenure <- bind_rows(out_r_GroupTenure, out_f_GroupTenure) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Group Tenure Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_GroupTenure_sum <- group_by(out_all_GroupTenure, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_GroupTenure %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_GroupTenure_sum, size = 1) +
  geom_text(data = mutate_if(out_all_GroupTenure_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 5), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = GroupTenure %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC04 ##### 
# insufficient data

##### VC05 Geographic Proximity #####
# Study-specific effects are deviations + average
out_r_GeogProx <- spread_draws(rem.brms.GeogProx, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_GeogProx <- spread_draws(rem.brms.GeogProx, b_Intercept) %>% 
  mutate(model_ID = "Geographic Proximity Effect") 

# Combine average and study-specific effects' data frames
out_all_GeogProx <- bind_rows(out_r_GeogProx, out_f_GeogProx) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Geographic Proximity Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_GeogProx_sum <- group_by(out_all_GeogProx, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_GeogProx %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_GeogProx_sum, size = 1) +
  geom_text(data = mutate_if(out_all_GeogProx_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 2), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = GeogProx %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC06 Member Management #####
# Study-specific effects are deviations + average
out_r_MemManage <- spread_draws(rem.brms.MemManage, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_MemManage <- spread_draws(rem.brms.MemManage, b_Intercept) %>% 
  mutate(model_ID = "Member Management Effect") 

# Combine average and study-specific effects' data frames
out_all_MemManage <- bind_rows(out_r_MemManage, out_f_MemManage) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Member Management Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_MemManage_sum <- group_by(out_all_MemManage, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_MemManage %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_MemManage_sum, size = 1) +
  geom_text(data = mutate_if(out_all_MemManage_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 5), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = MemManage %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC07 Group Sanctions #####
# Study-specific effects are deviations + average
out_r_GroupSanctions <- spread_draws(rem.brms.GroupSanctions, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_GroupSanctions <- spread_draws(rem.brms.GroupSanctions, b_Intercept) %>% 
  mutate(model_ID = "Group Sanctions Effect") 

# Combine average and study-specific effects' data frames
out_all_GroupSanctions <- bind_rows(out_r_GroupSanctions, out_f_GroupSanctions) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Group Sanctions Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_GroupSanctions_sum <- group_by(out_all_GroupSanctions, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_GroupSanctions %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_GroupSanctions_sum, size = 1) +
  geom_text(data = mutate_if(out_all_GroupSanctions_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 9), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = GroupSanctions %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC08 Peer Monitoring #####
# Study-specific effects are deviations + average
out_r_PeerMon <- spread_draws(rem.brms.PeerMon, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_PeerMon <- spread_draws(rem.brms.PeerMon, b_Intercept) %>% 
  mutate(model_ID = "Peer Monitoring Effect") 

# Combine average and study-specific effects' data frames
out_all_PeerMon <- bind_rows(out_r_PeerMon, out_f_PeerMon) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Peer Monitoring Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_PeerMon_sum <- group_by(out_all_PeerMon, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_PeerMon %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_PeerMon_sum, size = 1) +
  geom_text(data = mutate_if(out_all_PeerMon_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 4), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = PeerMon %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC09 External Monitoring #####
# Study-specific effects are deviations + average
out_r_ExtMonitoring <- spread_draws(rem.brms.ExtMonitoring, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_ExtMonitoring <- spread_draws(rem.brms.ExtMonitoring, b_Intercept) %>% 
  mutate(model_ID = "External Monitoring Effect") 

# Combine average and study-specific effects' data frames
out_all_ExtMonitoring <- bind_rows(out_r_ExtMonitoring, out_f_ExtMonitoring) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "External Monitoring Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_ExtMonitoring_sum <- group_by(out_all_ExtMonitoring, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_ExtMonitoring %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_ExtMonitoring_sum, size = 1) +
  geom_text(data = mutate_if(out_all_ExtMonitoring_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 3), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = ExtMonitoring %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC10 Group Size #####
# Study-specific effects are deviations + average
out_r_GroupSize <- spread_draws(rem.brms.GroupSize, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_GroupSize <- spread_draws(rem.brms.GroupSize, b_Intercept) %>% 
  mutate(model_ID = "Group Size Effect") 

# Combine average and study-specific effects' data frames
out_all_GroupSize <- bind_rows(out_r_GroupSize, out_f_GroupSize) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Group Size Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_GroupSize_sum <- group_by(out_all_GroupSize, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_GroupSize %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_GroupSize_sum, size = 1) +
  geom_text(data = mutate_if(out_all_GroupSize_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0, 3), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = GroupSize %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC11 Borrower Age #####
# Study-specific effects are deviations + average
out_r_BorrowerAge <- spread_draws(rem.brms.BorrowerAge, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_BorrowerAge <- spread_draws(rem.brms.BorrowerAge, b_Intercept) %>% 
  mutate(model_ID = "Borrower Age Effect") 

# Combine average and study-specific effects' data frames
out_all_BorrowerAge <- bind_rows(out_r_BorrowerAge, out_f_BorrowerAge) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Borrower Age Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_BorrowerAge_sum <- group_by(out_all_BorrowerAge, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_BorrowerAge %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_BorrowerAge_sum, size = 1) +
  geom_text(data = mutate_if(out_all_BorrowerAge_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0.8, 1.5), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = BorrowerAge %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### VC12 Borrower Sex #####
# Study-specific effects are deviations + average
out_r_BorrowerSex <- spread_draws(rem.brms.BorrowerSex, r_model_ID[model_ID,], b_Intercept) %>% 
  mutate(b_Intercept = r_model_ID + b_Intercept) 

# Average effect
out_f_BorrowerSex <- spread_draws(rem.brms.BorrowerSex, b_Intercept) %>% 
  mutate(model_ID = "Borrower Sex Effect") 

# Combine average and study-specific effects' data frames
out_all_BorrowerSex <- bind_rows(out_r_BorrowerSex, out_f_BorrowerSex) %>% 
  ungroup() %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T)) %>% 
  # Ensure that Average effect is on the bottom of the forest plot
  mutate(model_ID = fct_relevel(model_ID, "Borrower Sex Effect")) %>% 
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept))

# Data frame of summary numbers
out_all_BorrowerSex_sum <- group_by(out_all_BorrowerSex, model_ID) %>% 
  mean_qi(b_Intercept)

# Draw plot
out_all_BorrowerSex %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_BorrowerSex_sum, size = 1) +
  geom_text(data = mutate_if(out_all_BorrowerSex_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0.25, 3.5), n.breaks = 8) +
  #Add points for original study estimates
  geom_point(
    data = BorrowerSex %>% filter(risk_of_bias == "low" | risk_of_bias == "medium"), 
    aes(x=exp(yi)), position = position_nudge(y = -.05), shape = 1)

##### Combine all meta-analytic estimates #####


# Combine average and study-specific effects' data frames
out_all_metas <- bind_rows(out_f_Relatives, out_f_PrAcq, out_f_GroupTenure, 
                           out_f_GeogProx, out_f_MemManage, 
                           out_f_GroupSanctions, out_f_PeerMon, 
                           out_f_ExtMonitoring, out_f_GroupSize, 
                           out_f_BorrowerAge, out_f_BorrowerSex) %>% 
  ungroup() %>%
  # Exponentiate coefficients to provide odds ratios
  mutate(b_Intercept = exp(b_Intercept)) %>%
  # Order studies by magnitude of estimate
  mutate(model_ID = fct_reorder(model_ID, b_Intercept, .desc = T))

# Data frame of summary numbers
out_all_metas_sum <- group_by(out_all_metas, model_ID) %>% 
  median_qi(b_Intercept)

# Draw plot
out_all_metas %>%   
  ggplot(aes(b_Intercept, model_ID)) +
  geom_density_ridges(rel_min_height = 0.01, 
                      col = NA,
                      scale = 1, 
                      fill = "steelblue1") +
  geom_pointintervalh(data = out_all_metas_sum, size = 1) +
  geom_text(data = mutate_if(out_all_metas_sum, is.numeric, round, 2),
            # Use glue package to combine strings
            aes(label = glue::glue("{b_Intercept} [{.lower}, {.upper}]"), x = Inf),
            hjust = "inward") +
  # Add dashed line to 0
  geom_vline(xintercept=1, color = "black", linetype = "dashed") +
  # X axis label
  xlab("Odds Ratio") +
  # Y axis label
  ylab("Meta-Analytic Estimate") +
  # Theme
  theme_minimal() +
  # Set x axis limits
  scale_x_continuous(limits = c(0.2, 5.5), n.breaks = 8)
