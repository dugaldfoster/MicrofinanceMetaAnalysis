#####################################################################################
#
# R code to reproduce analysis for 'Testing Evolutionary Theories of Cooperation 
# via Meta-Analysis of Predictors of Repayment for Joint Liability Microfinance Loans'
#
# Written by Dugald Foster, contact dugaldwefoster"at"gmail.com
#
# OSF pre-registration available from XXX
#
#####################################################################################

# Install and load packages
#install.packages("pacman")
library(pacman)
p_load("devtools", "esc", "metafor", "brms", "brmstools", "ggplot2", 
       "dplyr", "bayesplot", "robvis", "RobustBayesianCopas")

# Read in effect size estimates for Variable Category x
VariableCategory <- read.csv("VariableCategory.csv")

# Inspect data
glimpse(VariableCategory)

##### Calculate summary effect measures #####
# B: Assemble log odds estimates
VariableCategory$yi <- as.numeric(ifelse(VariableCategory$effect_measure == 1 | 
                                      VariableCategory$effect_measure == 8,
                                    VariableCategory$effect_size, "NA"))

# Convert odds ratios to log odds
VariableCategory$yi <- ifelse(VariableCategory$effect_measure == 8, 
                         log(VariableCategory$yi),
                         VariableCategory$yi)

# Recode log odds to align estimates
VariableCategory$yi <- ifelse(VariableCategory$recode == 1, 
                         VariableCategory$yi * -1, 
                         VariableCategory$yi)

##### Impute missing standard errors #####
# Remove studies without relevant estimates
VariableCategory <- VariableCategory[!is.na(VariableCategory$yi), ]

#Relabel SE column in dataset as sei
names(VariableCategory)[names(VariableCategory) == "SE"] <- "sei"

#Convert to numeric
VariableCategory$sei <- as.numeric(VariableCategory$sei)

#Impute missing standard errors...
#From t values
VariableCategory$sei <- ifelse(is.na(VariableCategory$sei) & 
                                 !is.na(VariableCategory$t_value),
                          VariableCategory$effect_size/VariableCategory$t_value,
                          VariableCategory$sei)

#From z values
#If both t and z values are missing, first impute z values from p values
#Recode all p-values using threshold midpoints
VariableCategory$p_value <- as.numeric(recode(VariableCategory$p_value, 
                                              "<0.01" = "0.0075", 
                                         "<0.05" = "0.025", ">0.05" = "0.075", 
                                         "<0.1" = "0.075", ">0.1" = "0.125",
                                         "<0.15" = "0.125"))

#Impute z values
VariableCategory$z_value <- ifelse(is.na(VariableCategory$sei) & 
                                is.na(VariableCategory$t_value) &
                                is.na(VariableCategory$z_value),
                              abs(qnorm(VariableCategory$p_value/2)),
                              VariableCategory$z_value)

#Impute standard errors from z values
VariableCategory$sei <- ifelse(is.na(VariableCategory$sei) & 
                                 !is.na(VariableCategory$z_value),
                          VariableCategory$effect_size/VariableCategory$z_value,
                          VariableCategory$sei)

#Ensure standard errors are non-negative
VariableCategory$sei <- abs(VariableCategory$sei)

# Log the standard errors of odds ratios
VariableCategory$sei <- ifelse(VariableCategory$effect_measure == 8, 
                          log(VariableCategory$sei),
                          VariableCategory$sei)

##### Meta-analytic models #####

# Bayesian random effects model (low and medium risk of bias studies only)
rem.brms <- brm(yi | se(sei) ~ 1 + (1 | study),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(cauchy(0, 1), class = sd)),
                data = subset(VariableCategory, risk_of_bias == "low" | 
                                risk_of_bias == "medium"),
                iter = 50000,
                cores = 4,
                seed = 42,
                control = list(adapt_delta = 0.999, max_treedepth = 12))

# Trace plots
# Intercepts for mu and tau
plot(rem.brms, pars = c("^b_", "^sd_"))

# Assess convergence via posterior predictive checks
# Compare densities of sampled datasets to data
pp_check(rem.brms, nsamples = 100)

# Predictive intervals
pp_check(rem.brms, type = 'intervals', nsamples = 100)

# Compare distribution of median pp values in replicates against median pp value of the data
pp_check(rem.brms, type = "stat", stat = "median", nsamples = 1000) +
  legend_move("bottom")

# Check uniformity of distribution of LOO-PIT values
pp_check(rem.brms, type = "loo_pit")

# Model summary 
summary(rem.brms)

# Bayesian R2 model fit
bayes_R2(rem.brms, summary = T)

# Extract odds ratio for mu
exp(fixef(rem.brms))

# Get prediction intervals
predictive_interval(rem.brms)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms)

### Create density plots for mu and tau
# Extract posterior samples 
post.samples <- posterior_samples(rem.brms, c("^b", "^sd"), add_chain = T)
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
       y = element_blank()) +
  theme_minimal()

# Check the probability of mu < 0.20 based on empirical cumulative distribution function
mu.ecdf <- ecdf(post.samples$mu) 
plot(mu.ecdf)
mu.ecdf(0.2)

# Create forest plot of posterior distributions for each study's effect size
# Black dots show estimates and CIs of effect size based on the model
# White dots show observed effect sizes from studies
forest(rem.brms,
       show_data = T,
       sort = T,
       av_name = "Meta-Analytic Estimate") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  ggtitle("Group Tenure, Low/Moderate RoB Studies") 

##### Sensitivity analyses #####
# Alternative prior (a) for mu
rem.brms_mu_a <- brm(yi | se(sei) ~ 1 + (1 | study),
                prior = c(prior(normal(0, 10), class = Intercept),
                          prior(cauchy(0, 1), class = sd)),
                data = subset(VariableCategory, risk_of_bias == "low" | 
                                risk_of_bias == "medium"),
                iter = 50000,
                cores = 4,
                seed = 42,
                control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_mu_a)

# Model Summary
summary(rem.brms_mu_a)

# Alternative prior (b) for mu
rem.brms_mu_b <- brm(yi | se(sei) ~ 1 + (1 | study),
                   prior = c(prior(normal(0, 0.5), class = Intercept),
                             prior(cauchy(0, 1), class = sd)),
                   data = subset(VariableCategory, risk_of_bias == "low" | 
                                   risk_of_bias == "medium"),
                   iter = 50000,
                   cores = 4,
                   seed = 42,
                   control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_mu_b)

# Model Summary
summary(rem.brms_mu_b)

# Alternative prior (a) for tau
rem.brms_tau_a <- brm(yi | se(sei) ~ 1 + (1 | study),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(cauchy(0, .5), class = sd)),
                data = subset(VariableCategory, risk_of_bias == "low" | 
                                risk_of_bias == "medium"),
                iter = 50000,
                cores = 4,
                seed = 42,
                control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_tau_a)

# Model Summary
summary(rem.brms_tau_a)

# Alternative prior (b) for tau
rem.brms_tau_b <- brm(yi | se(sei) ~ 1 + (1 | study),
                      prior = c(prior(normal(0, 1), class = Intercept),
                                prior(cauchy(0, .3), class = sd)),
                      data = subset(VariableCategory, risk_of_bias == "low" | 
                                      risk_of_bias == "medium"),
                      iter = 50000,
                      cores = 4,
                      seed = 42,
                      control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_tau_b)

# Model Summary
summary(rem.brms_tau_b)
        
# Split studies by level of analysis (individual vs group outcomes)
# Repayment studies at the level of the borrower
rem.brms_borrowers <- brm(yi | se(sei) ~ 1 + (1 | study),
                          prior = c(prior(normal(0, 1), class = Intercept),
                                    prior(cauchy(0, 1), class = sd)),
                          data = subset(VariableCategory, risk_of_bias == "low" | 
                                          risk_of_bias == "medium" |
                                          level_of_analysis == 1),
                          iter = 50000,
                          cores = 4,
                          seed = 42,
                          control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_borrowers)

# Model Summary
summary(rem.brms_borrowers)

# Repayment studies at the level of the loan group
rem.brms_groups <- brm(yi | se(sei) ~ 1 + (1 | study),
                       prior = c(prior(normal(0, 1), class = Intercept),
                                 prior(cauchy(0, 1), class = sd)),
                       data = subset(VariableCategory, risk_of_bias == "low" | 
                                       risk_of_bias == "medium"|
                                       level_of_analysis == 2),
                       iter = 50000,
                       cores = 4,
                       seed = 42,
                       control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_groups)

# Model Summary
summary(rem.brms_groups)

# Test sensitivity of results to inclusion of studies at higher risk of bias
rem.brms_highrob <- brm(yi | se(sei) ~ 1 + (1 | study),
                        prior = c(prior(normal(0, 1), class = Intercept),
                                  prior(cauchy(0, 1), class = sd)),
                        data = subset(VariableCategory, risk_of_bias == "low" | 
                                        risk_of_bias == "medium" | 
                                        risk_of_bias == "high"),
                        iter = 50000,
                        cores = 4,
                        seed = 42,
                        control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_highrob)

# Model summary
summary(rem.brms_highrob)

# Get prediction intervals
predictive_interval(rem.brms_highrob)

# Extract random effects for each study, and their deviation from pooled effect
ranef(rem.brms_highrob)

### Create density plots for mu and tau
# Extract posterior samples 
post.samples <- posterior_samples(rem.brms_highrob, c("^b", "^sd"))
names(post.samples) <- c("mu", "tau")

# Denisty plot for posterior distributions of mu
ggplot(aes(x = mu), data = post.samples) +
  geom_density(fill = "lightblue", color = "lightblue", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samples$mu)) +
  labs(x = expression(mu),
       y = element_blank()) +
  theme_minimal()

# Denisty plot for posterior distributions of tau
ggplot(aes(x = tau), data = post.samples) +
  geom_density(fill = "lightgreen", color = "lightgreen", alpha = 0.7) +
  geom_point(y = 0, x = mean(post.samples$tau)) +
  labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()

# Check the probability that mu < 0.20 based on empirical cumulative distribution function
mu.ecdf <- ecdf(post.samples$mu) 
plot(mu.ecdf)
mu.ecdf(0.2)

# Create forest plot of posterior distributions for each study's effect size
forest(rem.brms_highrob,
       show_data = TRUE,
       av_name = "Meta-Analytic Estimate") +
  geom_vline(xintercept = 0, size = .25, lty = 2) +
  ggtitle("Group Tenure, All RoB Studies") 

##### Investigate heterogeneity via meta-regression #####
# Meta-regression 1: number of predictors in each model
rem.brms_nPredictors <- update(rem.brms,  ~ . + n_predictors,
                               prior = c(prior(normal(0, 1), class = Intercept),
                                         prior(cauchy(0, 1), class = sd)),
                               newdata = subset(VariableCategory, risk_of_bias == "low" | 
                                                  risk_of_bias == "medium"),
                               iter = 50000,
                               cores = 4,
                               seed = 42,
                               control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_nPredictors)

# Model summary
summary(rem.brms_nPredictors)

# Compare models using leave-one-out cross validation checks
loo(rem.brms, rem.brms_nPredictors)

# Plot Pareto k diagnostics
loo_rem.brms <- loo(rem.brms, save_psis = TRUE)
loo_rem.brms_nPredictors <- loo(rem.brms_nPredictors, save_psis = TRUE)

plot(loo_rem.brms)
plot(loo_rem.brms_nPredictors)

# Meta-regression 2: group size
rem.brms_grpsize <- update(rem.brms,  ~ . + group_size,
                               prior = c(prior(normal(0, 1), class = Intercept),
                                         prior(cauchy(0, 1), class = sd)),
                               newdata = subset(VariableCategory, risk_of_bias == "low" | 
                                                  risk_of_bias == "medium"),
                               iter = 50000,
                               cores = 4,
                               seed = 42,
                               control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_grpsize)

# Model summary
summary(rem.brms_grpsize)

# Compare models using leave-one-out cross validation checks
loo(rem.brms, rem.brms_grpsize)

# Meta-regression 3: borrower age
rem.brms_age <- update(rem.brms,  ~ . + age,
                               prior = c(prior(normal(0, 1), class = Intercept),
                                         prior(cauchy(0, 1), class = sd)),
                               newdata = subset(VariableCategory, risk_of_bias == "low" | 
                                                  risk_of_bias == "medium"),
                               iter = 50000,
                               cores = 4,
                               seed = 42,
                               control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_age)

# Model summary
summary(rem.brms_age)

# Compare models using leave-one-out cross validation checks
loo(rem.brms, rem.brms_age)

# Meta-regression 4: borrower sex
rem.brms_propmale <- update(rem.brms,  ~ . + prop_male,
                               prior = c(prior(normal(0, 1), class = Intercept),
                                         prior(cauchy(0, 1), class = sd)),
                               newdata = subset(VariableCategory, risk_of_bias == "low" | 
                                                  risk_of_bias == "medium"),
                               iter = 50000,
                               cores = 4,
                               seed = 42,
                               control = list(adapt_delta = 0.999, max_treedepth = 12))

# Assess convergence
pp_check(rem.brms_propmale)

# Model summary
summary(rem.brms_propmale)

# Compare models using leave-one-out cross validation checks
loo(rem.brms, rem.brms_propmale)

##### Risk of bias plots #####
VariableCategory_ROBOS <- read.csv("ROBOS_VariableCategory.csv")

rob_summary(data = VariableCategory_ROBOS,
            tool = "Generic",
            colour = "colourblind")

rob_traffic_light(VariableCategory_ROBOS, tool = "Generic")

#Alternatively, create these plots using the web app: https://mcguinlu.shinyapps.io/robvis/

##### Test for publication bias #####

### Tests which take frequentist model as input ###
# Contour-enhanced funnel plot
rem_freq <- rma(yi = yi, vi = ((sei^2)*n), data = VariableCategory) 

funnel(rem_freq, refline=0, level=c(90, 95, 99), 
       shade=c("white", "gray", "darkgray")) 

# Rank correlation test
ranktest(rem_freq)

### Bayesian Copas Selection Model ###

#Fit the RBC model with different random effects distributions 
#and compare them using the DIC
pubbias.some <- RobustBayesianCopas(y = VariableCategory$yi, s = VariableCategory$sei, 
                                    init=NULL, re.dist=c("normal"),
                                    burn=10000, nmc=10000)

pubbias.someST <- RobustBayesianCopas(y = VariableCategory$yi, s = VariableCategory$sei, 
                                      init=NULL, re.dist=c("StudentsT"),
                                      burn=10000, nmc=10000)

pubbias.someL <- RobustBayesianCopas(y = VariableCategory$yi, s = VariableCategory$sei, 
                                     init=NULL, re.dist=c("Laplace"),
                                     burn=10000, nmc=10000)

pubbias.someSL <- RobustBayesianCopas(y = VariableCategory$yi, s = VariableCategory$sei, 
                                      init=NULL, re.dist=c("slash"),
                                      burn=10000, nmc=10000)

pubbias.some$DIC
pubbias.someST$DIC
pubbias.someL$DIC
pubbias.someSL$DIC

# Plot posterior for rho
hist(pubbias.some$rho.samples)

# Point estimate for theta
pubbias.some$theta.hat 

# Standard error for theta
sd(pubbias.some$theta.samples)

# 95% posterior credible interval for theta
quantile(pubbias.some$theta.samples, probs=c(0.025,0.975))

# Obtain odds ratio estimates
OR.samples.RBC = exp(pubbias.some$theta.samples)

# Posterior mean OR
mean(OR.samples.RBC) 

# 95% posterior credible interval for OR
quantile(OR.samples.RBC, probs=c(0.025,0.975))

#Calculate D measure
#Model without publication bias
pubbias.none <- BayesNonBiasCorrected(y = VariableCategory$yi, s = VariableCategory$sei,
                                      init = NULL, re.dist = ("normal"),
                                      burn=10000, nmc=10000)

# D measure
D.measure(pubbias.some$theta.samples, pubbias.none$theta.samples)

### Repeat for Variable Categories x:n ###