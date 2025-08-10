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
BorrowerAge <- read.csv("VC11_BorrowerAge_All.csv")

# Remove punctuation in study names
BorrowerAge$study <- str_replace(BorrowerAge$study, ",", "")   
BorrowerAge$study <- str_replace(BorrowerAge$study, "_", " ")
BorrowerAge$study <- stri_replace_all_fixed(BorrowerAge$study, "(", "")
BorrowerAge$study <- stri_replace_all_fixed(BorrowerAge$study, ")", "")

# Inspect data
glimpse(BorrowerAge)

##### Calculate summary effect measures #####
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

##### Impute missing standard errors #####
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

##### Meta-analytic models #####

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

#Reveal stan code behind model
stancode(rem.brms.BorrowerAge)

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
plot(rem.brms.BorrowerAge, pars = c("^b_", "^sd_"))

# Assess convergence via posterior predictive checks
# Compare densities of sampled datasets to data
pp_check(rem.brms.BorrowerAge, nsamples = 100)

# Predictive intervals
pp_check(rem.brms.BorrowerAge, type = 'intervals', nsamples = 100)

# Compare distribution of median pp values in replicates against median pp value of the data
pp_check(rem.brms.BorrowerAge, type = "stat", stat = "median", nsamples = 1000) +
  legend_move("bottom")

# Check uniformity of distribution of LOO-PIT values
pp_check(rem.brms.BorrowerAge, type = "loo_pit")

# Model summary with 89% credible intervals
summary(rem.brms.BorrowerAge, prob = .89)

# Bayesian R2 model fit
bayes_R2(rem.brms.BorrowerAge, summary = T)

##### Extract odds ratio for mu with 95% credible interval #####
exp(fixef(rem.brms.BorrowerAge))

##### Assess I2 in addition to tau2? ###
# See https://academic.oup.com/jeea/article/18/6/3045/5908781?login=true for calculation

# Get prediction intervals
predictive_interval(rem.brms.BorrowerAge)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms.BorrowerAge)

### Create density plots for mu and tau
# Extract posterior samples 
post.samples <- posterior_samples(rem.brms.BorrowerAge, c("^b", "^sd"), add_chain = T)
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
hypothesis(rem.brms.BorrowerAge, "Intercept < 0.01")

# Create forest plot of posterior distributions for each study's effect size
# Black dots show estimates and CIs of effect size based on the model
# White dots show observed effect sizes from studies
forest(rem.brms.BorrowerAge,
       show_data = T,
       sort = T,
       av_name = "Meta-Analytic Estimate",
       fill = "steelblue1") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  scale_x_continuous(limits = c(-0.2, 0.4), n.breaks = 7)

##### Assess Heterogeneity #####

# Split studies by level of analysis (individual vs group outcomes)
rem.brms.BorrowerAge_borrowers <- brm(yi | se(sei) ~ 1 + (1 | study),
                                      prior = c(prior(normal(0, 1), class = Intercept),
                                                prior(cauchy(0, 1), class = sd)),
                                      data = subset(BorrowerAge, risk_of_bias == "low" | 
                                                      risk_of_bias == "medium" &
                                                      level_of_analysis == 1),
                                      iter = 50000,
                                      cores = 4,
                                      seed = 42,
                                      control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.BorrowerAge_borrowers)

rem.brms.BorrowerAge_groups <- brm(yi | se(sei) ~ 1 + (1 | study),
                                   prior = c(prior(normal(0, 1), class = Intercept),
                                             prior(cauchy(0, 1), class = sd)),
                                   data = subset(BorrowerAge, risk_of_bias == "low" | 
                                                   risk_of_bias == "medium" &
                                                   level_of_analysis == 2),
                                   iter = 50000,
                                   cores = 4,
                                   seed = 42,
                                   control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.BorrowerAge_groups)

# Investigate moderators
# number of predictors in each model
rem.brms.BorrowerAge_nPredictors <- brm(yi | se(sei) ~ 1 + n_predictors + (1 | study),
                                        prior = c(prior(normal(0, 1), class = Intercept),
                                                  prior(cauchy(0, 1), class = sd)),
                                        data = subset(BorrowerAge, risk_of_bias == "low" | 
                                                        risk_of_bias == "medium"),
                                        iter = 50000,
                                        cores = 4,
                                        seed = 42,
                                        control = list(adapt_delta = 0.999, max_treedepth = 12))

summary(rem.brms.BorrowerAge_nPredictors)

# proportion male in population 
rem.brms.BorrowerAge_sex <- brm(yi | se(sei) ~ 1 + prop_male + (1 | study),
                                prior = c(prior(normal(0, 1), class = Intercept),
                                          prior(cauchy(0, 1), class = sd)),
                                data = subset(BorrowerAge, risk_of_bias == "low" | 
                                                risk_of_bias == "medium"),
                                iter = 50000,
                                cores = 4,
                                seed = 42,
                                control = list(adapt_delta = 0.999, max_treedepth = 12))

# average age of population 
rem.brms.BorrowerAge_age <- brm(yi | se(sei) ~ 1 + age + (1 | study),
                                prior = c(prior(normal(0, 1), class = Intercept),
                                          prior(cauchy(0, 1), class = sd)),
                                data = subset(BorrowerAge, risk_of_bias == "low" | 
                                                risk_of_bias == "medium"),
                                iter = 50000,
                                cores = 4,
                                seed = 42,
                                control = list(adapt_delta = 0.999, max_treedepth = 12))

# Compare models using leave-one-out cross validation checks
loo(rem.brms.BorrowerAge, rem.brms.BorrowerAge_nPredictors)

# Plot Pareto k diagnostics
loo_rem.brms.BorrowerAge <- loo(rem.brms.BorrowerAge, save_psis = TRUE)
loo_rem.brms.BorrowerAge_nPredictors <- loo(rem.brms.BorrowerAge_nPredictors, save_psis = TRUE)

plot(loo_rem.brms.BorrowerAge)
plot(loo_rem.brms.BorrowerAge_nPredictors)

# Test sensitivity of results to inclusion of studies at higher risk of bias
rem.brms.BorrowerAge_highrob <- brm(yi | se(sei) ~ 1 + (1 | study),
                                    prior = c(prior(normal(0, 1), class = Intercept),
                                              prior(cauchy(0, 1), class = sd)),
                                    data = subset(BorrowerAge, risk_of_bias == "low" | 
                                                    risk_of_bias == "medium" | 
                                                    risk_of_bias == "high"),
                                    iter = 50000,
                                    cores = 4,
                                    seed = 42,
                                    control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms.BorrowerAge_highrob)

# Model summary
summary(rem.brms.BorrowerAge_highrob)

# Get prediction intervals
predictive_interval(rem.brms.BorrowerAge_highrob)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms.BorrowerAge_highrob)

# Extract odds ratio for mu with 95% credible interval 
exp(fixef(rem.brms.BorrowerAge_highrob))

### Create density plots for mu and tau
# Extract posterior samples 
post.samplesHR <- posterior_samples(rem.brms.BorrowerAge_highrob, c("^b", "^sd"))
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
forest(rem.brms.BorrowerAge_highrob,
       show_data = TRUE,
       av_name = "Meta-Analytic Estimate") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  ggtitle("Borrower Age, All RoB Studies")

##### Risk of bias plots #####
BorrowerAge_ROBOS <- read.csv("ROBOS_VC11.csv")

# See tool = "Generic" option: https://mcguinlu.shinyapps.io/robvis/
rob_summary(data = BorrowerAge_ROBOS,
            tool = "Generic",
            colour = "colourblind")

rob_traffic_light(BorrowerAge_ROBOS, tool = "Generic")

