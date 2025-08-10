#####################################################################################
#
# R code to reproduce analysis for 'Testing Evolutionary Theories of Cooperation 
# via Meta-Analysis of Predictors of Repayment for Joint Liability Microfinance Loans'
#
# Written by Dugald Foster, contact dugaldwefoster"at"gmail.com
#
# OSF pre-registration available from https://osf.io/3y9fu
#
#####################################################################################

# Install and load packages
#install.packages("pacman")
library(pacman)
p_load("devtools", "esc", "metafor", "brms", "brmstools", "ggplot2", 
       "dplyr", "bayesplot", "robvis", "RobustBayesianCopas", "stringr", "stringi")

# Read in effect size estimates for Variable Category x
GroupSize <- read.csv("VC10_GroupSize_All.csv")

# Remove punctuation in study names
GroupSize$study <- str_replace(GroupSize$study, ",", "")   
GroupSize$study <- str_replace(GroupSize$study, "_", " ")
GroupSize$study <- stri_replace_all_fixed(GroupSize$study, "(", "")
GroupSize$study <- stri_replace_all_fixed(GroupSize$study, ")", "")

# Inspect data
glimpse(GroupSize)

##### Calculate summary effect measures #####
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

##### Impute missing standard errors #####
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

##### Meta-analytic models #####

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

#Reveal stan code behind model
stancode(rem.brms.GroupSize)

#Plot priors
#mu
norm_ggplot <- ggplot(data = tibble(x = c(-20, 20)), aes(x = x)) +
  stat_function(fun = dnorm, n = 500, args = list(mean = 0, sd = 1)) +
  labs(title = "normal(0, 1)")

norm_ggplot

#tau
cauchy_ggplot <- ggplot(data = tibble(x = c(0, 5)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 500, args = list(location = 0, scale = 1)) +
  labs(title = "HalfCauchy(0, 1)") 

cauchy_ggplot

# Trace plots
# Intercepts for mu and tau
plot(rem.brms.GroupSize, pars = c("^b_", "^sd_"))

# Assess convergence via posterior predictive checks
# Compare densities of sampled datasets to data
pp_check(rem.brms.GroupSize, nsamples = 100)

# Predictive intervals
pp_check(rem.brms.GroupSize, type = 'intervals', nsamples = 100)

# Compare distribution of median pp values in replicates against median pp value of the data
pp_check(rem.brms.GroupSize, type = "stat", stat = "median", nsamples = 1000) +
  legend_move("bottom")

# Check uniformity of distribution of LOO-PIT values
pp_check(rem.brms.GroupSize, type = "loo_pit")

# Model summary with 89% credible intervals
summary(rem.brms.GroupSize, prob = .89)

# Bayesian R2 model fit
bayes_R2(rem.brms.GroupSize, summary = T)

##### Extract odds ratio for mu with 95% credible interval #####
exp(fixef(rem.brms.GroupSize))

##### Assess I2 in addition to tau2? ###
# See https://academic.oup.com/jeea/article/18/6/3045/5908781?login=true for calculation

# Get prediction intervals
predictive_interval(rem.brms.GroupSize)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms.GroupSize)

### Create density plots for mu and tau
# Extract posterior samples 
post.samples <- posterior_samples(rem.brms.GroupSize, c("^b", "^sd"), add_chain = T)
names(post.samples) <- c("mu", "tau", "chain", "iter")

# Denisty plot for posterior distributions of mu
ggplot(aes(x = mu), data = post.samples) +
  geom_density(fill = "lightblue", color = "lightblue", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samples$mu)) +
  labs(x = expression(mu),
       y = element_blank()) +
  theme_minimal()

# Density plot showing draws by chain
mcmc_dens_overlay(post.samples, pars = c("mu")) +
  theme(panel.grid = element_blank())

# Denisty plot for posterior distributions of tau
ggplot(aes(x = tau), data = post.samples) +
  geom_density(fill = "lightgreen", color = "lightgreen", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samples$tau)) +
  labs(x = expression(tau),
       y = "Density") +
  theme_minimal()

# Check the probability of mu < 0.20 based on empirical cumulative distribution function
mu.ecdf <- ecdf(post.samples$mu) 
plot(mu.ecdf)
mu.ecdf(0.2)

#Alternative?
hypothesis(rem.brms.GroupSize, "Intercept < 0.01")

# Create forest plot of posterior distributions for each study's effect size
# Black dots show estimates and CIs of effect size based on the model
# White dots show observed effect sizes from studies
forest(rem.brms.GroupSize,
       show_data = T,
       sort = T,
       av_name = "Meta-Analytic Estimate",
       fill = "steelblue1") +
  geom_vline(xintercept = 0, size = .25, lty = 2)

##### Assess Heterogeneity #####

# Split studies by level of analysis (individual vs group outcomes)
rem.brms.GroupSize_borrowers <- brm(yi | se(sei) ~ 1 + (1 | study),
                                         prior = c(prior(normal(0, 1), class = Intercept),
                                                   prior(cauchy(0, 1), class = sd)),
                                         data = subset(GroupSize, risk_of_bias == "low" | 
                                                         risk_of_bias == "medium" &
                                                         level_of_analysis == 1),
                                         iter = 50000,
                                         cores = 4,
                                         seed = 42,
                                         control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.GroupSize_borrowers)

rem.brms.GroupSize_groups <- brm(yi | se(sei) ~ 1 + (1 | study),
                                      prior = c(prior(normal(0, 1), class = Intercept),
                                                prior(cauchy(0, 1), class = sd)),
                                      data = subset(GroupSize, risk_of_bias == "low" | 
                                                      risk_of_bias == "medium" &
                                                      level_of_analysis == 2),
                                      iter = 50000,
                                      cores = 4,
                                      seed = 42,
                                      control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.GroupSize_groups)

# Investigate moderators
# number of predictors in each model
rem.brms.GroupSize_nPredictors <- brm(yi | se(sei) ~ 1 + n_predictors + (1 | study),
                                           prior = c(prior(normal(0, 1), class = Intercept),
                                                     prior(cauchy(0, 1), class = sd)),
                                           data = subset(GroupSize, risk_of_bias == "low" | 
                                                           risk_of_bias == "medium"),
                                           iter = 50000,
                                           cores = 4,
                                           seed = 42,
                                           control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.GroupSize_nPredictors)

# proportion male in population 
rem.brms.GroupSize_sex <- brm(yi | se(sei) ~ 1 + prop_male + (1 | study),
                                   prior = c(prior(normal(0, 1), class = Intercept),
                                             prior(cauchy(0, 1), class = sd)),
                                   data = subset(GroupSize, risk_of_bias == "low" | 
                                                   risk_of_bias == "medium"),
                                   iter = 50000,
                                   cores = 4,
                                   seed = 42,
                                   control = list(adapt_delta = 0.999, max_treedepth = 12))

# average age of population 
rem.brms.GroupSize_age <- brm(yi | se(sei) ~ 1 + age + (1 | study),
                                   prior = c(prior(normal(0, 1), class = Intercept),
                                             prior(cauchy(0, 1), class = sd)),
                                   data = subset(GroupSize, risk_of_bias == "low" | 
                                                   risk_of_bias == "medium"),
                                   iter = 50000,
                                   cores = 4,
                                   seed = 42,
                                   control = list(adapt_delta = 0.999, max_treedepth = 12))

# Compare models using leave-one-out cross validation checks
loo(rem.brms.GroupSize, rem.brms.GroupSize_nPredictors)

# Plot Pareto k diagnostics
loo_rem.brms.GroupSize <- loo(rem.brms.GroupSize, save_psis = TRUE)
loo_rem.brms.GroupSize_nPredictors <- loo(rem.brms.GroupSize_nPredictors, save_psis = TRUE)

plot(loo_rem.brms.GroupSize)
plot(loo_rem.brms.GroupSize_nPredictors)

# Test sensitivity of results to inclusion of studies at higher risk of bias
rem.brms.GroupSize_highrob <- brm(yi | se(sei) ~ 1 + (1 | study),
                                       prior = c(prior(normal(0, 1), class = Intercept),
                                                 prior(cauchy(0, 1), class = sd)),
                                       data = subset(GroupSize, risk_of_bias == "low" | 
                                                       risk_of_bias == "medium" | 
                                                       risk_of_bias == "high"),
                                       iter = 50000,
                                       cores = 4,
                                       seed = 42,
                                       control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms.GroupSize_highrob)

# Model summary
summary(rem.brms.GroupSize_highrob)

# Get prediction intervals
predictive_interval(rem.brms.GroupSize_highrob)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms.GroupSize_highrob)

# Extract odds ratio for mu with 95% credible interval 
exp(fixef(rem.brms.GroupSize_highrob))

### Create density plots for mu and tau
# Extract posterior samples 
post.samplesHR <- posterior_samples(rem.brms.GroupSize_highrob, c("^b", "^sd"))
names(post.samplesHR) <- c("mu", "tau")

# Denisty plot for posterior distributions of mu
ggplot(aes(x = mu), data = post.samplesHR) +
  geom_density(fill = "lightblue", color = "lightblue", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samplesHR$mu)) +
  labs(x = expression(mu),
       y = element_blank()) +
  theme_minimal()

# Denisty plot for posterior distributions of tau
ggplot(aes(x = tau), data = post.samplesHR) +
  geom_density(fill = "lightgreen", color = "lightgreen", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samplesHR$tau)) +
  labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()

# Check the probability of mu < 0.20 based on empirical cumulative distribution function
mu.ecdf <- ecdf(post.samplesHR$mu) 
plot(mu.ecdf)
mu.ecdf(0.2)

# Create forest plot of posterior distributions for each study's effect size
forest(rem.brms.GroupSize_highrob,
       show_data = TRUE,
       av_name = "Meta-Analytic Estimate") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  ggtitle("Group Size, All RoB Studies")

##### Risk of bias plots #####
GroupSize_ROBOS <- read.csv("ROBOS_VC10.csv")

# See tool = "Generic" option: https://mcguinlu.shinyapps.io/robvis/
rob_summary(data = GroupSize_ROBOS,
            tool = "Generic",
            colour = "colourblind")

rob_traffic_light(GroupSize_ROBOS, tool = "Generic")

