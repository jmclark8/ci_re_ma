### This script implements a simple comparison between our random effects CI-MA
### approach and that of Dahabreh et al. (2023). This is a very
### simple simulation setting designed to show at a high level how our approaches
### differ and when ideas from RE-CI-MA might be more appropriate.

### Loading packages

library(MASS) ### For mvrnorm
library(data.table) ### For organizing results into a table
library(here) ### For interacting with the working directory
library(ggplot2)
library(pracma)
library(magrittr)


# Defining the simulation -------------------------------------------------

sim_rep <- function(deltas, ### Size of potential outcome shift
                    sample_sizes, ### Number of individuals in each setting
                    ### (including both the trials and target pop)
                    combine_data ### Specify whether predictions should come 
                    ### from combined data or not
)
{
  
  ##### Saving important values #####
  
  ### Potential outcome shifts:
  delta_1 <- deltas[1]
  delta_2 <- deltas[2]
  delta_3 <- deltas[3]
  
  ### Setting-specific sample sizes:
  n_0 <- sample_sizes[1]
  n_1 <- sample_sizes[2]
  n_2 <- sample_sizes[3]
  n_3 <- sample_sizes[4]
  
  ### Sigma for the covariate draws
  vcov_matrix <- diag(x = 1, nrow = 3, ncol = 3)
  off_diag_entries <- row(vcov_matrix)-col(vcov_matrix) != 0
  vcov_matrix[off_diag_entries] <- 0.5
  
  ### Mean vectors for the covariate draws
  mu_0 <- 1*rep(1, 3)
  mu_1 <- 1.5*rep(1, 3)
  mu_2 <- 0.5*rep(1, 3)
  mu_3 <- 0*rep(1, 3)
  
  ##### Generating the data #####
  
  ### Covariate information
  
  target_pop_data <- setnames(data.table(mvrnorm(n = n_0, mu = mu_0, 
                                                 Sigma = vcov_matrix)), 
                              c("X_1", "X_2", "X_3"))
  study_1_data <- setnames(data.table(mvrnorm(n = n_1, mu = mu_1, 
                                              Sigma = vcov_matrix)), 
                           c("X_1", "X_2", "X_3"))
  study_2_data <- setnames(data.table(mvrnorm(n = n_2, mu = mu_2, 
                                              Sigma = vcov_matrix)), 
                           c("X_1", "X_2", "X_3"))
  study_3_data <- setnames(data.table(mvrnorm(n = n_3, mu = mu_3, 
                                              Sigma = vcov_matrix)), 
                           c("X_1", "X_2", "X_3"))
  
  ### Outcome information
  
  study_1_data[, y_a := 0.5+0.5*X_1+0.5*X_2+0.5*X_3+delta_1+rnorm(n_1)]
  study_2_data[, y_a := 0.5+0.5*X_1+0.5*X_2+0.5*X_3+delta_2+rnorm(n_2)]
  study_3_data[, y_a := 0.5+0.5*X_1+0.5*X_2+0.5*X_3+delta_3+rnorm(n_3)]
  
  ### Creating a combined dataset for the Dahabreh et al. method and adding a trial
  ### participation indicator:
  
  combined_trial_data <- rbindlist(list(study_1_data, study_2_data, study_3_data))
  
  ### Fitting the outcome models. (Note that we focus here on average
  ### potential outcomes alone, not average differences of potential
  ### outcomes.)
  
  study_1_model <- lm(y_a ~ X_1+X_2+X_3, data = study_1_data)
  study_2_model <- lm(y_a ~ X_1+X_2+X_3, data = study_2_data)
  study_3_model <- lm(y_a ~ X_1+X_2+X_3, data = study_3_data)
  combined_model <- lm(y_a ~ X_1+X_2+X_3, data = combined_trial_data)
  
  ### Adding the predictions to the target population data
  
  target_pop_data[, study_1_pred := predict(study_1_model, 
                                            newdata = .SD[, .(X_1, X_2, X_3)])]
  
  target_pop_data[, study_2_pred := predict(study_2_model, 
                                            newdata = .SD[, .(X_1, X_2, X_3)])]
  
  target_pop_data[, study_3_pred := predict(study_3_model, 
                                            newdata = .SD[, .(X_1, X_2, X_3)])]
  
  
  ### This prediction is a very simple analogue to Dahabreh's method
  target_pop_data[, combined_data_pred := predict(combined_model, 
                                                  newdata = .SD[, .(X_1, X_2, X_3)])]
  
  ### Now applying our method, which averages over the predictions from the
  ### study-specific models
  
  target_pop_data[, re_pred := rowMeans(.SD), .SDcols = c("study_1_pred", 
                                                          "study_2_pred", 
                                                          "study_3_pred")]
  
  if (combine_data == FALSE){
    
    y_a_pred <- target_pop_data[, mean(re_pred)]
    return(y_a_pred)
    
  }
  
  if (combine_data == TRUE){
    
    y_a_pred <- target_pop_data[, mean(combined_data_pred)]
    return(y_a_pred)
    
  }
}


# Running the simulation and compiling results ----------------------------

### Function to produce sims and results

sim_results_fcn <- function(deltas, sample_sizes, combine_data){
  
  set.seed(1)
  
  sims <- replicate(10000, sim_rep(deltas = deltas, 
                                   sample_sizes = sample_sizes, 
                                   combine_data = combine_data))
  
  est_mean <- mean(sims)
  
  return(est_mean)
}

### Generating the table

sim_scenarios <- data.table("delta_1" = c(1.5, 1.5, 1.5, 1.5, 0, 0), 
                            "delta_2" = c(0.5, 0.5, 0.75, 0.75, 0, 0), 
                            "delta_3" = c(5, 5, 1.25, 1.25, 0, 0), 
                            "n_0" = c(100, 100, 100, 100, 100, 100), 
                            "n_1" = c(100, 100, 100, 100, 100, 100), 
                            "n_2" = c(100, 100, 100, 100, 100, 100), 
                            "n_3" = c(400, 400, 100, 100, 100, 100), 
                            combine_data = c(TRUE, FALSE, TRUE, FALSE, 
                                             TRUE, FALSE))

sim_scenarios[, scenario := paste0("Scenario ", .I)]

### Computing the mean bias and rmse

sim_results <- sim_scenarios[, sim_results_fcn(deltas = c(delta_1, 
                                                          delta_2, delta_3), 
                                               sample_sizes = c(n_0, n_1, 
                                                                n_2, n_3), 
                                               combine_data = combine_data), 
                             by = scenario]

sim_table <- merge(sim_scenarios, sim_results, 
                   by = "scenario")

### Need to do some re-shaping:

sim_table_wide <- dcast(sim_table[, !c("scenario")], 
                        delta_1 + delta_2 + delta_3 +
                          n_0 + n_1 + n_2 + n_3 ~ combine_data, 
                        value.var = c("V1"))

sim_table_rounded <- copy(sim_table_wide) %>%
  .[, c("FALSE", "TRUE") := lapply(.SD, 
                                   function(x){format(round(x, digits = 3), nsmall = 3)}),
                                   .SDcols = c("FALSE", "TRUE")]

### Saving the results for later:

# saveRDS(object = sim_table_wide,
#          file = here("estimator_comparison_res_raw.RDS"))

### Also computing the true values:

true_values_combine <- copy(sim_scenarios) %>%
  .[, total_trial_size := n_1+n_2+n_3] %>%
  .[combine_data == TRUE, true_value := (n_1/total_trial_size)*(1.5+.5+delta_1)+
      (n_2/total_trial_size)*(1.5+.5+delta_2)+
      (n_3/total_trial_size)*(1.5+.5+delta_3)]

### Doing so for the Dahabreh method is more complex:

get_true_dahabreh <- function(delta_1, delta_2, delta_3, 
                              n_1, n_2, n_3){
  
  
  ### First, saving the true coefficients for each study's outcome model
  coef_1 <- matrix(c(0.5+delta_1, .5, .5, .5), nrow = 4, ncol = 1)
  coef_2 <- matrix(c(0.5+delta_2, .5, .5, .5), nrow = 4, ncol = 1)
  coef_3 <- matrix(c(0.5+delta_3, .5, .5, .5), nrow = 4, ncol = 1)
  
  ### Saving the sample size proportions
  
  total_trial_n <- n_1+n_2+n_3
  p_1 <- n_1/total_trial_n
  p_2 <- n_2/total_trial_n
  p_3 <- n_3/total_trial_n
  ### Now, saving the true weights for computing the true coefficients of the 
  ### combined outcome model
  
  X_1_prod <- matrix(c(1, 1.5, 1.5, 1.5, 
                       1.5, 3.25, 2.75, 2.75,
                       1.5, 2.75, 3.25, 2.75,
                       1.5, 2.75, 2.75, 3.25),
                     ncol = 4, nrow = 4, byrow = TRUE)
  
  X_2_prod <- matrix(c(1.00, 0.50, 0.50, 0.50,
                       0.50, 1.25, 0.75, 0.75,
                       0.50, 0.75, 1.25, 0.75,
                       0.50, 0.75, 0.75, 1.25),
                     ncol = 4, nrow = 4, byrow = TRUE)
  
  X_3_prod <- matrix(c(1.00, 0.00, 0.00, 0.00, 
                       0.00, 1.00, 0.50, 0.50,
                       0.00, 0.50, 1.00, 0.50,
                       0.00, 0.50, 0.50, 1.00),
                     ncol = 4, nrow = 4, byrow = TRUE)
  
  X_prod <- p_1*X_1_prod+p_2*X_2_prod+p_3*X_3_prod
  X_prod_inv <- inv(X_prod)
  
  ### Computing weights for the overall coefficient:
  W_1 <- p_1*X_prod_inv%*%X_1_prod
  W_2 <- p_2*X_prod_inv%*%X_2_prod
  W_3 <- p_3*X_prod_inv%*%X_3_prod
  overall_coeff <- W_1%*%coef_1+W_2%*%coef_2+W_3%*%coef_3
  
  ### Plugging in the mean covariate values for the target population
  true_value <- matrix(c(1, 1, 1, 1), nrow = 1, ncol = 4)%*%overall_coeff %>%
    as.numeric()
  true_value_floor <- floor(true_value)
  true_value_resid <- true_value-true_value_floor
  true_value_frac <- true_value_resid %>%
    fractions() %>%
    as.character()
  
  true_value_char <- paste0(true_value_floor, " + ", true_value_frac)
  
  return(true_value_char)
  
  
  
}

true_values_dahabreh <- copy(sim_scenarios) %>%
  .[combine_data == TRUE] %>%
  .[, .(delta_1, delta_2, delta_3, n_1, n_2, n_3, scenario)] %>%
  .[, get_true_dahabreh(delta_1 = delta_1, delta_2 = delta_2, 
                        delta_3 = delta_3, 
                        n_1 = n_1, n_2 = n_2, n_3 = n_3), 
    by = "scenario"]



































