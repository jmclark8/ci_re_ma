### As described in the README, our simulation study was run using computing
### resources available at our institution. This following script defines the full parameter
### grid used in our simulation study and includes code necessary to run a single
### iteration for a given set of parameters, i.e., one row from that grid. In 
### practice, we used the `foreach` package with %dopar% to run many such iterations in 
### parallel in our HPC environment.

### We can provide the exact seeds necessary to reproduce our simulation results
### upon request.

library(here)

source(here("main_sim_study", "helper_fcns", "outcome_model_preds.R"))
source(here("main_sim_study", "helper_fcns", "gamma_squared_interval.R"))
source(here("main_sim_study", "helper_fcns","simple_bootstrap_interval.R"))
source(here("main_sim_study", "helper_fcns","wild_bootstrap_interval.R"))

library(MASS)
library(data.table)
library(magrittr)
library(EnvStats)

# Setting up parameter grid -----------------------------------------------

param_grid <- expand.grid(n_studies  = c(5, 10, 15, 30, 50), 
                          var_dist = c("Uniform", "Normal", "Exponential", 
                                       "Pareto"))

# Code necessary for one iteration ----------------------------------------

### Generating data for the studies and the target population. 

### Saving the potential outcome shifts. 

if (var_dist == "Uniform"){
  
  deltas <- runif(n = n_studies + 1, min = -2, max = 2)
  
}

if (var_dist == "Normal"){
  
  deltas <- rnorm(n = n_studies + 1, mean = 0, sd = 1)
  
}

if (var_dist == "Exponential"){
  
  deltas <- rexp(n = n_studies + 1, rate = 1)-1
  
}

if (var_dist == "Pareto"){
  
  deltas <- rpareto(n = n_studies + 1, location = 1, shape = 3)
  
}


delta_0 <- deltas[1]
study_deltas <- deltas[2:length(deltas)]

### Defining the covariate sigma matrix that will be common across settings
vcov_matrix <- diag(x = 1, nrow = 3, ncol = 3)
off_diag_entries <- row(vcov_matrix)-col(vcov_matrix) != 0
vcov_matrix[off_diag_entries] <- 0.5

### The mean vectors for the covariates in the studies will be determined
### via an equally spaced partition of the interval [0, 1.5], with number 
### of partitions equal to the number of studies.

mu_0 <- 1*rep(1, 3)
study_mus <- seq(from = 0, to = 1.5, length.out = n_studies) %>%
  ### We need 3 identical mu values for each study, since the three covariates
  ### for each person in a given study are assumed to arise from a multivariate
  ### normal distribution with a mean vector of the form (c, c, c), for scalar
  ### c.
  lapply(., function(x){x*rep(1, 3)})

### Generating the data

### Covariate information

target_pop_data <- setnames(data.table(mvrnorm(n = 1000, mu = mu_0, 
                                               Sigma = vcov_matrix)), 
                            c("X_1", "X_2", "X_3"))

study_data <- lapply(1:n_studies, 
                     function(x){setnames(data.table(mvrnorm(n = 100, 
                                                             mu = study_mus[[x]], 
                                                             Sigma = vcov_matrix), 
                                                     ### Adding the delta value
                                                     ### will make things easier
                                                     ### when generating the outcome
                                                     delta_value = study_deltas[x]), 
                                          c("X_1", "X_2", "X_3", "delta_value"))})

### Assigning participants in each study to treatment or control. Note
### that we focus on outcomes under treatment, and the delta deviations are
### only relevant to such outcomes. 

invisible(lapply(study_data, 
                 function(x){x[, assignment := rbinom(n = .N, 
                                                      size = 1, 
                                                      prob = 0.5)]}))

### Adding potential outcome data to each dataset. 

target_pop_data[, y_a := 0.5+0.5*X_1+0.5*X_2+0.5*X_3+delta_0+rnorm(.N)]

true_target_outcome <- 2+delta_0

### Note that `data.table` allows us to modify the tables in-place, 
### so we can add the outcome to our `study_data`object directly.

### Again, these are potential outcomes under treatment; the key thing, 
### of course, is that we actually observe such outcomes in the treatment
### group for those in a given study!

invisible(lapply(study_data, 
                 function(x){x[, y_a := 0.5+0.5*X_1+0.5*X_2+0.5*X_3+delta_value+rnorm(.N)]}))

### Constructing a prediction interval using the bootstrapped values
### Passing copies to these functions so that we can make changes to the argument within
### each function w/o changing the underlying data. These functions
### are included in the relevant helper function scripts.

pred_interval_gamma <- get_interval_gamma(target_data = target_pop_data,
                                          trial_data = study_data, 
                                          number_studies = n_studies)

pred_interval_simple <- get_interval_simple(target_data = copy(target_pop_data), 
                                            trial_data = lapply(study_data, copy),
                                            number_studies = n_studies)

pred_interval_wild <- get_interval_wild(target_data = copy(target_pop_data), 
                                        trial_data = lapply(study_data, copy),
                                        number_studies = n_studies)

########################################
########################################
##
##        aggregate all performance 
##        metrics and setting info
##
########################################
########################################

res_tmp <- data.frame(seed               = seed, 
                      sim_rep            = s,
                      n_studies          = n_studies,
                      var_dist           = var_dist,
                      pred_interval_low_wild_quant  = pred_interval_wild$`Wild Bootstrap Interval Quantiles.lwr`, 
                      pred_interval_high_wild_quant = pred_interval_wild$`Wild Bootstrap Interval Quantiles.upper`,
                      pred_interval_low_wild_norm  = pred_interval_wild$`Wild Bootstrap Interval Normal Approx.lwr`, 
                      pred_interval_high_wild_norm = pred_interval_wild$`Wild Bootstrap Interval Normal Approx.upper`,
                      pred_interval_low_simple_quant = pred_interval_simple$`Simple Bootstrap Interval Quantiles.lwr`, 
                      pred_interval_high_simple_quant = pred_interval_simple$`Simple Bootstrap Interval Quantiles.upper`,
                      pred_interval_low_simple_norm = pred_interval_simple$`Simple Bootstrap Interval Normal Approx.lwr`, 
                      pred_interval_high_simple_norm = pred_interval_simple$`Simple Bootstrap Interval Normal Approx.upper`,
                      pred_interval_low_gamma = pred_interval_gamma$`Gamma Squared Interval.lwr`, 
                      pred_interval_high_gamma = pred_interval_gamma$`Gamma Squared Interval.upper`,
                      truth = true_target_outcome)

