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
BorrowerSex <- read.csv("VC12_BorrowerSex_All.csv")

# Remove punctuation in study names
BorrowerSex$study <- str_replace(BorrowerSex$study, ",", "")   
BorrowerSex$study <- str_replace(BorrowerSex$study, "_", " ")
BorrowerSex$study <- stri_replace_all_fixed(BorrowerSex$study, "(", "")
BorrowerSex$study <- stri_replace_all_fixed(BorrowerSex$study, ")", "")

# Inspect data
glimpse(BorrowerSex)

##### Calculate summary effect measures #####
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

##### Impute missing standard errors #####
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

##### Meta-analytic models #####

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

#Reveal stan code behind model
stancode(rem.brms.BorrowerSex)

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
plot(rem.brms.BorrowerSex, pars = c("^b_", "^sd_"))

# Assess convergence via posterior predictive checks
# Compare densities of sampled datasets to data
pp_check(rem.brms.BorrowerSex, nsamples = 100)

# Predictive intervals
pp_check(rem.brms.BorrowerSex, type = 'intervals', nsamples = 100)

# Compare distribution of median pp values in replicates against median pp value of the data
pp_check(rem.brms.BorrowerSex, type = "stat", stat = "median", nsamples = 1000) +
  legend_move("bottom")

# Check uniformity of distribution of LOO-PIT values
pp_check(rem.brms.BorrowerSex, type = "loo_pit")

# Model summary with 89% credible intervals
summary(rem.brms.BorrowerSex, prob = .89)

# Bayesian R2 model fit
bayes_R2(rem.brms.BorrowerSex, summary = T)

##### Extract odds ratio for mu with 95% credible interval #####
exp(fixef(rem.brms.BorrowerSex))

##### Assess I2 in addition to tau2? ###
# See https://academic.oup.com/jeea/article/18/6/3045/5908781?login=true for calculation

# Get prediction intervals
predictive_interval(rem.brms.BorrowerSex)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms.BorrowerSex)

### Create density plots for mu and tau
# Extract posterior samples 
post.samples <- posterior_samples(rem.brms.BorrowerSex, c("^b", "^sd"), add_chain = T)
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
hypothesis(rem.brms.BorrowerSex, "Intercept < 0.01")

# Create forest plot of posterior distributions for each study's effect size
# Black dots show estimates and CIs of effect size based on the model
# White dots show observed effect sizes from studies
forest(rem.brms.BorrowerSex,
       show_data = T,
       sort = T,
       av_name = "Meta-Analytic Estimate",
       fill = "steelblue1") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  scale_x_continuous(limits = c(-1, 3), n.breaks = 7)

##### Assess Heterogeneity #####

# Split studies by level of analysis (individual vs group outcomes)
rem.brms.BorrowerSex_borrowers <- brm(yi | se(sei) ~ 1 + (1 | study),
                                    prior = c(prior(normal(0, 1), class = Intercept),
                                              prior(cauchy(0, 1), class = sd)),
                                    data = subset(BorrowerSex, risk_of_bias == "low" | 
                                                    risk_of_bias == "medium" &
                                                    level_of_analysis == 1),
                                    iter = 50000,
                                    cores = 4,
                                    seed = 42,
                                    control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.BorrowerSex_borrowers)

rem.brms.BorrowerSex_groups <- brm(yi | se(sei) ~ 1 + (1 | study),
                                 prior = c(prior(normal(0, 1), class = Intercept),
                                           prior(cauchy(0, 1), class = sd)),
                                 data = subset(BorrowerSex, risk_of_bias == "low" | 
                                                 risk_of_bias == "medium" &
                                                 level_of_analysis == 2),
                                 iter = 50000,
                                 cores = 4,
                                 seed = 42,
                                 control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.BorrowerSex_groups)

# Investigate moderators
# number of predictors in each model
rem.brms.BorrowerSex_nPredictors <- brm(yi | se(sei) ~ 1 + n_predictors + (1 | study),
                                      prior = c(prior(normal(0, 1), class = Intercept),
                                                prior(cauchy(0, 1), class = sd)),
                                      data = subset(BorrowerSex, risk_of_bias == "low" | 
                                                      risk_of_bias == "medium"),
                                      iter = 50000,
                                      cores = 4,
                                      seed = 42,
                                      control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.BorrowerSex_nPredictors)

# proportion male in population 
rem.brms.BorrowerSex_sex <- brm(yi | se(sei) ~ 1 + prop_male + (1 | study),
                              prior = c(prior(normal(0, 1), class = Intercept),
                                        prior(cauchy(0, 1), class = sd)),
                              data = subset(BorrowerSex, risk_of_bias == "low" | 
                                              risk_of_bias == "medium"),
                              iter = 50000,
                              cores = 4,
                              seed = 42,
                              control = list(adapt_delta = 0.999, max_treedepth = 12))

# average age of population 
rem.brms.BorrowerSex_age <- brm(yi | se(sei) ~ 1 + age + (1 | study),
                              prior = c(prior(normal(0, 1), class = Intercept),
                                        prior(cauchy(0, 1), class = sd)),
                              data = subset(BorrowerSex, risk_of_bias == "low" | 
                                              risk_of_bias == "medium"),
                              iter = 50000,
                              cores = 4,
                              seed = 42,
                              control = list(adapt_delta = 0.999, max_treedepth = 12))

# Compare models using leave-one-out cross validation checks
loo(rem.brms.BorrowerSex, rem.brms.BorrowerSex_nPredictors)

# Plot Pareto k diagnostics
loo_rem.brms.BorrowerSex <- loo(rem.brms.BorrowerSex, save_psis = TRUE)
loo_rem.brms.BorrowerSex_nPredictors <- loo(rem.brms.BorrowerSex_nPredictors, save_psis = TRUE)

plot(loo_rem.brms.BorrowerSex)
plot(loo_rem.brms.BorrowerSex_nPredictors)

# Test sensitivity of results to inclusion of studies at higher risk of bias
rem.brms.BorrowerSex_highrob <- brm(yi | se(sei) ~ 1 + (1 | study),
                                  prior = c(prior(normal(0, 1), class = Intercept),
                                            prior(cauchy(0, 1), class = sd)),
                                  data = subset(BorrowerSex, risk_of_bias == "low" | 
                                                  risk_of_bias == "medium" | 
                                                  risk_of_bias == "high"),
                                  iter = 50000,
                                  cores = 4,
                                  seed = 42,
                                  control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms.BorrowerSex_highrob)

# Model summary
summary(rem.brms.BorrowerSex_highrob)

# Get prediction intervals
predictive_interval(rem.brms.BorrowerSex_highrob)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms.BorrowerSex_highrob)

# Extract odds ratio for mu with 95% credible interval 
exp(fixef(rem.brms.BorrowerSex_highrob))

### Create density plots for mu and tau
# Extract posterior samples 
post.samplesHR <- posterior_samples(rem.brms.BorrowerSex_highrob, c("^b", "^sd"))
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
forest(rem.brms.BorrowerSex_highrob,
       show_data = TRUE,
       av_name = "Meta-Analytic Estimate") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  ggtitle("Borrower Sex, All RoB Studies")

##### Risk of bias plots #####
BorrowerSex_ROBOS <- read.csv("ROBOS_VC12.csv")

# See tool = "Generic" option: https://mcguinlu.shinyapps.io/robvis/
rob_summary(data = BorrowerSex_ROBOS,
            tool = "Generic",
            colour = "colourblind")

rob_traffic_light(BorrowerSex_ROBOS, tool = "Generic")

